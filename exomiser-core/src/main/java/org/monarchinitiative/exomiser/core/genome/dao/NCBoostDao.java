/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2018 Queen Mary University of London.
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

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.monarchinitiative.exomiser.core.genome.dao;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import htsjdk.tribble.readers.TabixReader;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.NCBoostScore;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.cache.annotation.Cacheable;
import org.springframework.cache.annotation.Caching;

import java.io.IOException;

/**
 *
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public class NCBoostDao implements PathogenicityDao {

    private final Logger logger = LoggerFactory.getLogger(NCBoostDao.class);

    private final TabixDataSource ncboostTabixDataSource;

    public NCBoostDao(TabixDataSource ncboostTabixDataSource) {
        this.ncboostTabixDataSource = ncboostTabixDataSource;
    }

    @Caching(cacheable = {
            @Cacheable(cacheNames = "hg19.ncboost", keyGenerator = "variantKeyGenerator", condition = "#variant.genomeAssembly == T(org.monarchinitiative.exomiser.core.genome.GenomeAssembly).HG19"),
            @Cacheable(cacheNames = "hg38.ncboost", keyGenerator = "variantKeyGenerator", condition = "#variant.genomeAssembly == T(org.monarchinitiative.exomiser.core.genome.GenomeAssembly).HG38"),
    })
    @Override
    public PathogenicityData getPathogenicityData(Variant variant) {
        logger.debug("Getting NCBoost data for {}", variant);
        // NCBoost has not been trained on missense variants so skip these
        if (variant.getVariantEffect() == VariantEffect.MISSENSE_VARIANT) {
            return PathogenicityData.empty();
        }
        return processResults(variant);
    }

    private PathogenicityData processResults(Variant variant) {
        String chromosome = variant.getChromosomeName();
        int start = variant.getPosition();
        int end = calculateEndPosition(variant);
        return getNCBoostData(chromosome, start, end);
    }

    private int calculateEndPosition(Variant variant) {
        int pos = variant.getPosition();

        //we're doing this here in order not to have to count all this each time we need the value
        int refLength = variant.getRef().length();
        int altLength = variant.getAlt().length();
        //What about MNV?
        if (refLength == altLength) {
            return pos;
        }
        //these end positions are calculated according to recommendation by Max and Peter who produced the NCBoost score
        //don't change this unless they say.
        if (isDeletion(refLength, altLength)) {
            // test all deleted bases (being 1-based we need to correct the length)
            return pos + refLength - 1;
        } else if (isInsertion(refLength, altLength)) {
            // test bases either side of insertion
            return pos + 1;
        }
        return pos;
    }

    private static boolean isDeletion(int refLength, int altLength) {
        return refLength > altLength;
    }

    private static boolean isInsertion(int refLength, int altLength) {
        return refLength < altLength;
    }

    private synchronized PathogenicityData getNCBoostData(String chromosome, int start, int end) {
        try {
            float score = Float.NaN;
            String line;
            TabixReader.Iterator results = ncboostTabixDataSource.query(chromosome + ":" + start + "-" + end);
            while ((line = results.next()) != null) {
                String[] elements = line.split("\t");
                if (Float.isNaN(score)) {
                    score = Float.parseFloat(elements[5]);
                } else {
                    score = Math.max(score, Float.parseFloat(elements[5]));
                }
            }
            if (!Float.isNaN(score)) {
                return PathogenicityData.of(NCBoostScore.of(score));
            }
        } catch (IOException e) {
            logger.error("Unable to read from NCBoost tabix file {}", ncboostTabixDataSource.getSource(), e);
        }
        return PathogenicityData.empty();
    }

}
