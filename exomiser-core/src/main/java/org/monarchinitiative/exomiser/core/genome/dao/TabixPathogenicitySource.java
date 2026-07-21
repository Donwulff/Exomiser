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

import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource;

/**
 * Declarative description of a tabix-backed pathogenicity source, consumed by {@link TabixPathogenicityDao}.
 * <p>
 * A single generic DAO parameterised by one of these replaces the previous one-class-per-source approach, so that
 * adding a new tabix pathogenicity source (e.g. DANN, NCBOOST) is a matter of declaring a config preset rather than
 * duplicating tabix-parsing logic. The score column is read as the <em>raw</em> value and handed to
 * {@link PathogenicitySource#name()}'s {@code PathogenicityScore.of(source, raw)} factory, which applies whatever
 * per-source scaling is required - so no scaling logic lives here.
 *
 * @see #cadd()
 * @see #remm()
 * @see #dann()
 * @see #ncboost()
 * @since 15.1.0
 */
public record TabixPathogenicitySource(
        PathogenicitySource source,
        ReadMode readMode,
        Applicability applicability,
        int refColumn,
        int altColumn,
        int scoreColumn
) {

    /**
     * Sentinel value for {@link #scoreColumn()} meaning "use the last column of each tabix row". This handles data
     * files whose score is the final column regardless of how many annotation columns precede it - e.g. CADD is
     * distributed both as a trimmed 6-column file and as a full annotation TSV, and in both the PHRED score is last.
     * Solving that mismatch without editing the DAO was the original motivation for making the source configurable.
     */
    public static final int LAST_COLUMN = -1;

    public enum ReadMode {
        /** Match {@code ref}/{@code alt} at a single position (CADD, DANN). */
        EXACT_MATCH,
        /** Take the maximum score over every base spanned by the variant (REMM, NCBOOST). */
        RANGE_MAX
    }

    public enum Applicability {
        /** Score every variant (CADD, DANN). */
        ALL,
        /** Score only non-coding variants (REMM, NCBOOST); these predictors are not trained on coding variants. */
        NON_CODING;

        boolean test(Variant variant) {
            return this == ALL || variant.isNonCodingVariant();
        }
    }

    /**
     * @return {@code true} if this source should be consulted for the given variant.
     */
    public boolean appliesTo(Variant variant) {
        return applicability.test(variant);
    }

    // Presets for the currently-shipped tabix sources. Column indices are 0-based.

    /**
     * CADD: {@code #Chrom Pos Ref Alt RawScore PHRED}. The PHRED score is the last column in both the trimmed and the
     * full-annotation distributions, so {@link #LAST_COLUMN} is used deliberately.
     */
    public static TabixPathogenicitySource cadd() {
        return new TabixPathogenicitySource(PathogenicitySource.CADD, ReadMode.EXACT_MATCH, Applicability.ALL, 2, 3, LAST_COLUMN);
    }

    /**
     * DANN (CADD-style, SNV + indel): {@code Chrom Pos Ref Alt score}, score already scaled 0-1 in the final column.
     */
    public static TabixPathogenicitySource dann() {
        return new TabixPathogenicitySource(PathogenicitySource.DANN, ReadMode.EXACT_MATCH, Applicability.ALL, 2, 3, LAST_COLUMN);
    }

    /**
     * REMM: {@code Chrom Pos Score}; per-base non-coding score taken as the max over the variant's range.
     */
    public static TabixPathogenicitySource remm() {
        return new TabixPathogenicitySource(PathogenicitySource.REMM, ReadMode.RANGE_MAX, Applicability.NON_CODING, -1, -1, 2);
    }

    /**
     * NCBOOST: {@code Chrom Pos ... Score}; per-base non-coding score taken as the max over the variant's range.
     */
    public static TabixPathogenicitySource ncboost() {
        return new TabixPathogenicitySource(PathogenicitySource.NCBOOST, ReadMode.RANGE_MAX, Applicability.NON_CODING, -1, -1, LAST_COLUMN);
    }
}
