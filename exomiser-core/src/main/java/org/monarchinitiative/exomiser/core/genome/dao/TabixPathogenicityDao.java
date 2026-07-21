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

/**
 * Generic tabix-backed {@link PathogenicityDao} driven by a {@link TabixPathogenicitySource} config. This is the single
 * code path behind all tabix pathogenicity sources (CADD, REMM, DANN, NCBOOST, ...); the per-source subclasses that
 * exist alongside it ({@link CaddDao}, {@link RemmDao}, {@link DannDao}, {@link NCBoostDao}) add nothing but a config
 * preset and their own {@code @Cacheable} cache name.
 * <p>
 * The raw score column is parsed and passed straight to {@link PathogenicityScore#of}, which applies any per-source
 * scaling. Because every such scaling is monotonic, {@link TabixPathogenicitySource.ReadMode#RANGE_MAX} can safely take
 * the max of the raw column values before scaling.
 *
 * @since 15.1.0
 */
public class TabixPathogenicityDao implements PathogenicityDao {

    private final Logger logger = LoggerFactory.getLogger(TabixPathogenicityDao.class);

    private final TabixPathogenicitySource config;
    // For EXACT_MATCH sources SNVs and indels are held in separate files; RANGE_MAX sources have a single file, in which
    // case both references point at the same source.
    private final TabixDataSource snvTabixDataSource;
    private final TabixDataSource indelTabixDataSource;

    /**
     * Constructor for {@link TabixPathogenicitySource.ReadMode#EXACT_MATCH} sources with separate SNV and indel files.
     */
    public TabixPathogenicityDao(TabixPathogenicitySource config, TabixDataSource snvTabixDataSource, TabixDataSource indelTabixDataSource) {
        this.config = config;
        this.snvTabixDataSource = snvTabixDataSource;
        this.indelTabixDataSource = indelTabixDataSource;
    }

    /**
     * Constructor for single-file sources (typically {@link TabixPathogenicitySource.ReadMode#RANGE_MAX}).
     */
    public TabixPathogenicityDao(TabixPathogenicitySource config, TabixDataSource tabixDataSource) {
        this(config, tabixDataSource, tabixDataSource);
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
        TabixDataSource tabixDataSource = variant.variantType() == VariantType.SNV ? snvTabixDataSource : indelTabixDataSource;
        String ref = variant.ref();
        String alt = variant.alt();
        int start = variant.start();
        try {
            TabixReader.Iterator results = tabixDataSource.query(variant.contigName(), start, start);
            String line;
            //there can be 0 - N results at a position, one per ref/alt allele
            while ((line = results.next()) != null) {
                String[] elements = line.split("\t");
                if (elements[config.refColumn()].equals(ref) && elements[config.altColumn()].equals(alt)) {
                    return toPathogenicityData(parseScore(elements));
                }
            }
        } catch (IOException e) {
            logger.error("Unable to read from {} tabix file {}", config.source(), tabixDataSource.getSource(), e);
        }
        return PathogenicityData.empty();
    }

    private synchronized PathogenicityData rangeMax(Variant variant) {
        String chromosome = variant.contigName();
        int start = variant.start();
        // test bases either side of an insertion, or all bases of a SNV/MNV/DEL
        int end = variant.variantType() == VariantType.INS ? variant.end() + 1 : variant.end();
        try {
            float score = Float.NaN;
            String line;
            TabixReader.Iterator results = snvTabixDataSource.query(chromosome, start, end);
            while ((line = results.next()) != null) {
                String[] elements = line.split("\t");
                float value = parseScore(elements);
                score = Float.isNaN(score) ? value : Math.max(score, value);
            }
            if (!Float.isNaN(score)) {
                return toPathogenicityData(score);
            }
        } catch (IOException e) {
            logger.error("Unable to read from {} tabix file {}", config.source(), snvTabixDataSource.getSource(), e);
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
