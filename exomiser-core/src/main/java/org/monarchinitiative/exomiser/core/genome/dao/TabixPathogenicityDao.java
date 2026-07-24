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

package org.monarchinitiative.exomiser.core.genome.dao;

import htsjdk.tribble.readers.TabixReader;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityScore;
import org.monarchinitiative.svart.VariantType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.List;
import java.util.Objects;

/**
 * Generic tabix-backed {@link PathogenicityDao} driven by a {@link TabixPathogenicitySource} config. This is the single
 * code path behind all tabix pathogenicity sources (CADD, REMM, DANN, NCBOOST, ...); the per-source subclasses that
 * exist alongside it ({@link CaddDao}, {@link RemmDao}, {@link DannDao}, {@link NCBoostDao}) add nothing but a config
 * preset and their own {@code @Cacheable} cache name.
 * <p>
 * The raw score column is parsed and passed straight to {@link PathogenicityScore#of}, which applies any per-source
 * scaling. Because every such scaling is monotonic, {@link TabixPathogenicitySource.ReadMode#RANGE_MAX} can safely take
 * the max of the raw column values before scaling, and {@link TabixPathogenicitySource.ReadMode#EXACT_MATCH} can take
 * the max across multiple tabix files (e.g. official CADD plus sample-specific novel scores under one CADD identity).
 *
 * @since 15.1.0
 */
public class TabixPathogenicityDao implements PathogenicityDao {

    private final Logger logger = LoggerFactory.getLogger(TabixPathogenicityDao.class);

    private final TabixPathogenicitySource config;
    // For EXACT_MATCH sources SNVs and indels are held in separate file lists; RANGE_MAX sources typically have a
    // single file, in which case both lists are the same one-element list. Multiple files are queried and the max
    // matching raw score is kept (one PathogenicitySource identity).
    private final List<TabixDataSource> snvTabixDataSources;
    private final List<TabixDataSource> indelTabixDataSources;

    /**
     * Constructor for {@link TabixPathogenicitySource.ReadMode#EXACT_MATCH} sources with separate SNV and indel file
     * lists. Each list is queried in order; when more than one file hits the same allele the highest raw score wins.
     */
    public TabixPathogenicityDao(TabixPathogenicitySource config, List<TabixDataSource> snvTabixDataSources, List<TabixDataSource> indelTabixDataSources) {
        this.config = Objects.requireNonNull(config);
        this.snvTabixDataSources = List.copyOf(Objects.requireNonNull(snvTabixDataSources));
        this.indelTabixDataSources = List.copyOf(Objects.requireNonNull(indelTabixDataSources));
        if (this.snvTabixDataSources.isEmpty() || this.indelTabixDataSources.isEmpty()) {
            throw new IllegalArgumentException("Tabix pathogenicity sources for " + config.source() + " must not be empty lists");
        }
    }

    /**
     * Convenience constructor for a single SNV file and a single indel file.
     */
    public TabixPathogenicityDao(TabixPathogenicitySource config, TabixDataSource snvTabixDataSource, TabixDataSource indelTabixDataSource) {
        this(config, List.of(snvTabixDataSource), List.of(indelTabixDataSource));
    }

    /**
     * Constructor for single-file sources (typically {@link TabixPathogenicitySource.ReadMode#RANGE_MAX}).
     */
    public TabixPathogenicityDao(TabixPathogenicitySource config, TabixDataSource tabixDataSource) {
        this(config, List.of(tabixDataSource), List.of(tabixDataSource));
    }

    @Override
    public PathogenicityData getPathogenicityData(Variant variant) {
        logger.debug("Getting {} data for {}", config.source(), variant);
        if (!config.appliesTo(variant)) {
            return PathogenicityData.empty();
        }
        return switch (config.readMode()) {
            case EXACT_MATCH -> exactMatch(variant);
            case RANGE_MAX -> rangeMax(variant);
        };
    }

    private synchronized PathogenicityData exactMatch(Variant variant) {
        List<TabixDataSource> sources = variant.variantType() == VariantType.SNV ? snvTabixDataSources : indelTabixDataSources;
        String ref = variant.ref();
        String alt = variant.alt();
        int start = variant.start();
        String contig = variant.contigName();
        float best = Float.NaN;
        for (TabixDataSource tabixDataSource : sources) {
            try {
                TabixReader.Iterator results = tabixDataSource.query(contig, start, start);
                String line;
                // there can be 0 - N results at a position, one per ref/alt allele
                while ((line = results.next()) != null) {
                    String[] elements = line.split("\t");
                    if (elements[config.refColumn()].equals(ref) && elements[config.altColumn()].equals(alt)) {
                        float value = parseScore(elements);
                        best = Float.isNaN(best) ? value : Math.max(best, value);
                        // one matching allele per file is enough; other files may still raise the max
                        break;
                    }
                }
            } catch (IOException e) {
                logger.error("Unable to read from {} tabix file {}", config.source(), tabixDataSource.getSource(), e);
            }
        }
        if (!Float.isNaN(best)) {
            return toPathogenicityData(best);
        }
        return PathogenicityData.empty();
    }

    private synchronized PathogenicityData rangeMax(Variant variant) {
        String chromosome = variant.contigName();
        int start = variant.start();
        // test bases either side of an insertion, or all bases of a SNV/MNV/DEL
        int end = variant.variantType() == VariantType.INS ? variant.end() + 1 : variant.end();
        float score = Float.NaN;
        // RANGE_MAX historically used a single file; if multiple are configured, max across all of them too
        for (TabixDataSource tabixDataSource : snvTabixDataSources) {
            try {
                String line;
                TabixReader.Iterator results = tabixDataSource.query(chromosome, start, end);
                while ((line = results.next()) != null) {
                    String[] elements = line.split("\t");
                    float value = parseScore(elements);
                    score = Float.isNaN(score) ? value : Math.max(score, value);
                }
            } catch (IOException e) {
                logger.error("Unable to read from {} tabix file {}", config.source(), tabixDataSource.getSource(), e);
            }
        }
        if (!Float.isNaN(score)) {
            return toPathogenicityData(score);
        }
        return PathogenicityData.empty();
    }

    private float parseScore(String[] elements) {
        int column = config.scoreColumn() == TabixPathogenicitySource.LAST_COLUMN ? elements.length - 1 : config.scoreColumn();
        return Float.parseFloat(elements[column]);
    }

    private PathogenicityData toPathogenicityData(float rawScore) {
        return PathogenicityData.of(PathogenicityScore.of(config.source(), rawScore));
    }
}
