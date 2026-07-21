/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2021 Queen Mary University of London.
 * Copyright (c) 2012-2016 Charité Universitätsmedizin Berlin and Genome Research Ltd.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


package org.monarchinitiative.exomiser.core.genome;

import jakarta.annotation.Nonnull;
import org.monarchinitiative.exomiser.core.model.GeneStatistics;
import org.monarchinitiative.exomiser.core.genome.dao.*;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencyData;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityScore;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource;
import org.monarchinitiative.svart.GenomicInterval;
import org.monarchinitiative.svart.GenomicVariant;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.concurrent.CompletableFuture;

import static org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource.*;

/**
 * Default implementation of the VariantDataService.
 *
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public class VariantDataServiceImpl implements VariantDataService {

    private static final Logger logger = LoggerFactory.getLogger(VariantDataServiceImpl.class);

    // Sources which only score non-coding variants and so should not be consulted for coding variants. This mirrors the
    // per-source Applicability declared on the tabix DAO configs, but is needed here too because callers may register a
    // mock/stub DAO in place of the real one.
    private static final Set<PathogenicitySource> NON_CODING_ONLY_SOURCES = EnumSet.of(REMM, NCBOOST);

    private final VariantWhiteList whiteList;
    // Default data sources
    private final FrequencyDao defaultFrequencyDao;
    private final PathogenicityDao defaultPathogenicityDao;
    private final ClinVarDao clinVarDao;

    // Optional data sources
    private final FrequencyDao localFrequencyDao;
    // Registry of tabix-backed pathogenicity DAOs keyed by the source they provide. Adding a new tabix source is a
    // matter of registering it on the builder - no code in this class needs to change (bar NON_CODING_ONLY_SOURCES if
    // the source is non-coding-only).
    private final Map<PathogenicitySource, PathogenicityDao> tabixPathogenicityDaos;

    // Structural variant data sources
    private final FrequencyDao svFrequencyDao;
    private final PathogenicityDao svPathogenicityDao;

    private VariantDataServiceImpl(Builder builder) {

        this.whiteList = Objects.requireNonNull(builder.variantWhiteList);

        this.defaultFrequencyDao = Objects.requireNonNull(builder.defaultFrequencyDao, "defaultFrequencyDao required!");
        this.defaultPathogenicityDao = Objects.requireNonNull(builder.defaultPathogenicityDao, "defaultPathogenicityDao required!");
        this.clinVarDao = Objects.requireNonNull(builder.clinVarDao, "clinVarDao required!");

        this.localFrequencyDao = builder.localFrequencyDao;
        this.tabixPathogenicityDaos = new LinkedHashMap<>(builder.tabixPathogenicityDaos);

        this.svFrequencyDao = builder.svFrequencyDao;
        this.svPathogenicityDao = builder.svPathogenicityDao;
    }

    @Override
    public boolean variantIsWhiteListed(Variant variant) {
        return whiteList.contains(variant);
    }

    @Override
    public FrequencyData getVariantFrequencyData(Variant variant, Set<FrequencySource> frequencySources) {

        if (isStructural(variant)) {
            return svFrequencyDao.getFrequencyData(variant);
        }
        // This could be run alongside the pathogenicities as they are all stored in the same datastore
        FrequencyData defaultFrequencyData = defaultFrequencyDao.getFrequencyData(variant);

        FrequencyData.Builder frequencyDataBuilder = defaultFrequencyData.toBuilder();
        frequencyDataBuilder.filterSources(frequencySources);

        if (frequencySources.contains(FrequencySource.LOCAL)) {
            FrequencyData localFrequencyData = localFrequencyDao.getFrequencyData(variant);
            frequencyDataBuilder.mergeFrequencyData(localFrequencyData);
        }

        return frequencyDataBuilder.build();
    }

    // PacBio data contains lots of longer non-symbolic variants with an SVTYPE
    // so our working definition of 'structural' is any symbolic allele or allele over 50 bp
    private boolean isStructural(Variant variant) {
        return variant.isSymbolic() || variant.length() >= 50;
    }

    @Override
    public PathogenicityData getVariantPathogenicityData(Variant variant, Set<PathogenicitySource> pathogenicitySources) {

        if (isStructural(variant)) {
            return svPathogenicityDao.getPathogenicityData(variant);
        }

        ClinVarData clinVarData = clinVarDao.getClinVarData(variant);

        if (pathogenicitySources.isEmpty()) {
            // Fast-path for the unlikely case when no sources are defined - we'll just return the ClinVar data
            return PathogenicityData.of(clinVarData);
        }

        PathogenicityData defaultPathogenicityData;
        List<PathogenicityScore> allPathScores = new ArrayList<>();
        if (containsTabixSource(pathogenicitySources)) {
            CompletableFuture<PathogenicityData> futureDefaultData = CompletableFuture.supplyAsync(() -> defaultPathogenicityDao
                    .getPathogenicityData(variant));
            // run async - tabix sources are slow compared to MVStore
            List<CompletableFuture<PathogenicityData>> futurePathData = new ArrayList<>();
            // Fire each requested tabix source that applies to this variant. Non-coding-only sources (REMM, NCBOOST) are
            // skipped for coding variants - they are trained on the non-coding genome, where they outperform CADD.
            for (Map.Entry<PathogenicitySource, PathogenicityDao> entry : tabixPathogenicityDaos.entrySet()) {
                PathogenicitySource source = entry.getKey();
                if (pathogenicitySources.contains(source) && sourceAppliesToVariant(source, variant)) {
                    PathogenicityDao pathogenicityDao = entry.getValue();
                    futurePathData.add(CompletableFuture.supplyAsync(() -> pathogenicityDao.getPathogenicityData(variant)));
                }
            }
            for (CompletableFuture<PathogenicityData> pathogenicityDataCompletableFuture : futurePathData) {
                PathogenicityData pathogenicityData = pathogenicityDataCompletableFuture.join();
                allPathScores.addAll(pathogenicityData.pathogenicityScores());
            }
            defaultPathogenicityData = futureDefaultData.join();
        } else {
            defaultPathogenicityData = defaultPathogenicityDao.getPathogenicityData(variant);
        }

        // Previously we deliberately ignored synonymous variants from dbNSFP as these shouldn't be there
        // e.g. ?assembly=hg37&chr=1&start=158581087&ref=G&alt=A has a MutationTaster score of 1
        // However, cryptic splice variants can be synonymous variants, so these should be retained.
        addAllWantedScores(pathogenicitySources, defaultPathogenicityData, allPathScores);

        return PathogenicityData.of(clinVarData, allPathScores);
    }

    private boolean containsTabixSource(Set<PathogenicitySource> pathogenicitySources) {
        for (PathogenicitySource source : tabixPathogenicityDaos.keySet()) {
            if (pathogenicitySources.contains(source)) {
                return true;
            }
        }
        return false;
    }

    private static boolean sourceAppliesToVariant(PathogenicitySource source, Variant variant) {
        return !NON_CODING_ONLY_SOURCES.contains(source) || variant.isNonCodingVariant();
    }

    private void addAllWantedScores(Set<PathogenicitySource> pathogenicitySources, PathogenicityData defaultPathogenicityData, List<PathogenicityScore> allPathScores) {
        for (PathogenicityScore score : defaultPathogenicityData.pathogenicityScores()) {
            if (pathogenicitySources.contains(score.source())) {
                allPathScores.add(score);
            }
        }
    }

    @Override
    public ClinVarData getClinVarData(@Nonnull Variant variant) {
        return clinVarDao.getClinVarData(variant);
    }

    @Override
    public ClinVarData getClinVarData(@Nonnull GenomicVariant genomicVariant) {
        return clinVarDao.getClinVarData(genomicVariant);
    }

    @Override
    public Map<GenomicVariant, ClinVarData> findClinVarRecordsOverlappingInterval(@Nonnull GenomicInterval genomicInterval) {
        return clinVarDao.findClinVarRecordsOverlappingInterval(genomicInterval);
    }

    @Override
    public GeneStatistics getGeneStatistics(@Nonnull String geneSymbol) {
        return clinVarDao.getGeneStatistics(geneSymbol);
    }

    public static Builder builder() {
        return new Builder();
    }

    public static class Builder {

        private VariantWhiteList variantWhiteList = InMemoryVariantWhiteList.empty();

        private ClinVarDao clinVarDao;

        private FrequencyDao defaultFrequencyDao;
        private PathogenicityDao defaultPathogenicityDao;

        private FrequencyDao localFrequencyDao;

        private final Map<PathogenicitySource, PathogenicityDao> tabixPathogenicityDaos = new LinkedHashMap<>();

        private FrequencyDao svFrequencyDao = new StubFrequencyDao();
        private PathogenicityDao svPathogenicityDao = new StubPathogenicityDao();

        public Builder variantWhiteList(VariantWhiteList variantWhiteList) {
            this.variantWhiteList = variantWhiteList;
            return this;
        }

        public Builder clinVarDao(ClinVarDao clinVarDao) {
            this.clinVarDao = clinVarDao;
            return this;
        }

        public Builder defaultFrequencyDao(FrequencyDao defaultFrequencyDao) {
            this.defaultFrequencyDao = defaultFrequencyDao;
            return this;
        }

        public Builder defaultPathogenicityDao(PathogenicityDao defaultPathogenicityDao) {
            this.defaultPathogenicityDao = defaultPathogenicityDao;
            return this;
        }

        public Builder localFrequencyDao(FrequencyDao localFrequencyDao) {
            this.localFrequencyDao = localFrequencyDao;
            return this;
        }

        /**
         * Registers a tabix-backed pathogenicity DAO under the given source. This is the generic entry point; adding a
         * new tabix pathogenicity source needs nothing more than a call to this. Passing a {@code null} DAO leaves the
         * source unregistered (so it simply contributes no scores) rather than failing at query time.
         */
        public Builder pathogenicityDao(PathogenicitySource source, PathogenicityDao pathogenicityDao) {
            if (pathogenicityDao != null) {
                this.tabixPathogenicityDaos.put(source, pathogenicityDao);
            }
            return this;
        }

        public Builder caddDao(PathogenicityDao caddDao) {
            return pathogenicityDao(PathogenicitySource.CADD, caddDao);
        }

        public Builder remmDao(PathogenicityDao remmDao) {
            return pathogenicityDao(PathogenicitySource.REMM, remmDao);
        }

        public Builder dannDao(PathogenicityDao dannDao) {
            return pathogenicityDao(PathogenicitySource.DANN, dannDao);
        }

        public Builder ncBoostDao(PathogenicityDao ncBoostDao) {
            return pathogenicityDao(PathogenicitySource.NCBOOST, ncBoostDao);
        }

        public Builder testPathScoreDao(PathogenicityDao testPathScoreDao) {
            return pathogenicityDao(PathogenicitySource.TEST, testPathScoreDao);
        }

        public Builder svFrequencyDao(FrequencyDao svFrequencyDao) {
            this.svFrequencyDao = svFrequencyDao;
            return this;
        }

        public Builder svPathogenicityDao(PathogenicityDao svPathogenicityDao) {
            this.svPathogenicityDao = svPathogenicityDao;
            return this;
        }

        public VariantDataServiceImpl build() {
            return new VariantDataServiceImpl(this);
        }
    }

    private static class StubFrequencyDao implements FrequencyDao {
        @Override
        public FrequencyData getFrequencyData(Variant variant) {
            return FrequencyData.empty();
        }
    }

    private static class StubPathogenicityDao implements PathogenicityDao {
        @Override
        public PathogenicityData getPathogenicityData(Variant variant) {
            return PathogenicityData.empty();
        }
    }
}
