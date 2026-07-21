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

import de.charite.compbio.jannovar.annotation.VariantEffect;
import htsjdk.tribble.readers.TabixReader;
import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.genome.TestFactory;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.CaddScore;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityScore;

import java.util.Arrays;
import java.util.List;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource.DANN;
import static org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource.NCBOOST;

/**
 * Tests for the generic tabix pathogenicity engine. Deliberately Mockito-free (a plain fake {@link TabixDataSource}
 * feeds canned rows) so it exercises the parsing logic without mocking the concrete {@code TabixReader}.
 */
class TabixPathogenicityDaoTest {

    /** A {@link TabixDataSource} that returns the same canned tab-separated rows for any query. */
    private static TabixDataSource tabixReturning(String... rows) {
        List<String> lines = Arrays.asList(rows);
        return new TabixDataSource() {
            @Override
            public TabixReader.Iterator query(String query) {
                return MockTabixIterator.of(lines);
            }

            @Override
            public TabixReader.Iterator query(String chromosome, int start, int end) {
                return MockTabixIterator.of(lines);
            }

            @Override
            public String getSource() {
                return "test";
            }

            @Override
            public void close() {
                // no-op
            }
        };
    }

    private static VariantEvaluation snv(String ref, String alt, VariantEffect effect) {
        return TestFactory.variantBuilder(1, 2, ref, alt).variantEffect(effect).build();
    }

    @Test
    void exactMatchReadsLastColumnForFullCaddAnnotationRow() {
        // A full CADD annotation TSV has many columns before the trailing PHRED score; LAST_COLUMN must pick it out,
        // which is the whole reason the DAO was made configurable.
        TabixDataSource tabix = tabixReturning("1\t2\tA\tT\tANNO\t-0.234\t25.3");
        TabixPathogenicityDao instance = new TabixPathogenicityDao(TabixPathogenicitySource.cadd(), tabix);

        PathogenicityData result = instance.getPathogenicityData(snv("A", "T", VariantEffect.MISSENSE_VARIANT));

        assertThat(result, equalTo(PathogenicityData.of(CaddScore.of(25.3f))));
    }

    @Test
    void exactMatchTrimmedSixColumnCaddRowStillWorks() {
        // The trimmed distribution ends in PHRED too, so LAST_COLUMN is identical to the old fixed column 5.
        TabixDataSource tabix = tabixReturning("1\t2\tA\tT\t-0.234\t3.45");
        TabixPathogenicityDao instance = new TabixPathogenicityDao(TabixPathogenicitySource.cadd(), tabix);

        PathogenicityData result = instance.getPathogenicityData(snv("A", "T", VariantEffect.MISSENSE_VARIANT));

        assertThat(result, equalTo(PathogenicityData.of(CaddScore.of(3.45f))));
    }

    @Test
    void exactMatchReturnsEmptyWhenNoAlleleMatches() {
        TabixDataSource tabix = tabixReturning("1\t2\tA\tC\tANNO\t-0.234\t25.3");
        TabixPathogenicityDao instance = new TabixPathogenicityDao(TabixPathogenicitySource.dann(), tabix);

        PathogenicityData result = instance.getPathogenicityData(snv("A", "T", VariantEffect.MISSENSE_VARIANT));

        assertThat(result, equalTo(PathogenicityData.empty()));
    }

    @Test
    void dannExactMatchUsesAlreadyScaledScore() {
        TabixDataSource tabix = tabixReturning("1\t2\tA\tT\t0.98");
        TabixPathogenicityDao instance = new TabixPathogenicityDao(TabixPathogenicitySource.dann(), tabix);

        PathogenicityData result = instance.getPathogenicityData(snv("A", "T", VariantEffect.MISSENSE_VARIANT));

        assertThat(result, equalTo(PathogenicityData.of(PathogenicityScore.of(DANN, 0.98f))));
    }

    @Test
    void ncBoostRangeMaxTakesHighestSpannedScore() {
        TabixDataSource tabix = tabixReturning("1\t1\t0.1", "1\t2\t0.9", "1\t3\t0.3");
        TabixPathogenicityDao instance = new TabixPathogenicityDao(TabixPathogenicitySource.ncboost(), tabix);

        PathogenicityData result = instance.getPathogenicityData(snv("A", "T", VariantEffect.REGULATORY_REGION_VARIANT));

        assertThat(result, equalTo(PathogenicityData.of(PathogenicityScore.of(NCBOOST, 0.9f))));
    }

    @Test
    void nonCodingSourceReturnsEmptyForCodingVariant() {
        TabixDataSource tabix = tabixReturning("1\t2\t0.9");
        TabixPathogenicityDao instance = new TabixPathogenicityDao(TabixPathogenicitySource.ncboost(), tabix);

        // a coding (missense) variant must never be scored by a non-coding-only source, regardless of the data present
        PathogenicityData result = instance.getPathogenicityData(snv("A", "T", VariantEffect.MISSENSE_VARIANT));

        assertThat(result, equalTo(PathogenicityData.empty()));
    }
}
