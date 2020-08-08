/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2019 Queen Mary University of London.
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

import htsjdk.tribble.readers.TabixReader;
import org.monarchinitiative.exomiser.core.model.AllelePosition;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.pathogenicity.DannScore;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.cache.annotation.Cacheable;
import org.springframework.cache.annotation.Caching;

import java.io.IOException;
 
/**
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class DannDao implements PathogenicityDao {
 
    private final Logger logger = LoggerFactory.getLogger(DannDao.class);

    private final TabixDataSource dannInDelTabixDataSource;
    private final TabixDataSource dannSnvTabixDataSource;

    public DannDao(TabixDataSource dannInDelTabixDataSource, TabixDataSource dannSnvTabixDataSource) {
        this.dannInDelTabixDataSource = dannInDelTabixDataSource;
        this.dannSnvTabixDataSource = dannSnvTabixDataSource;
    }

    @Caching(cacheable = {
            @Cacheable(cacheNames = "hg19.dann", keyGenerator = "variantKeyGenerator", condition = "#variant.genomeAssembly == T(org.monarchinitiative.exomiser.core.genome.GenomeAssembly).HG19"),
            @Cacheable(cacheNames = "hg38.dann", keyGenerator = "variantKeyGenerator", condition = "#variant.genomeAssembly == T(org.monarchinitiative.exomiser.core.genome.GenomeAssembly).HG38"),
    })
    @Override
    public PathogenicityData getPathogenicityData(Variant variant) {
        logger.debug("Getting DANN data for {}", variant);
        return processResults(variant);
    }

    private PathogenicityData processResults(Variant variant) {
        String chromosome = variant.getChromosomeName();
        String ref = variant.getRef();
        String alt = variant.getAlt();
        int start = variant.getPosition();
        if (AllelePosition.isSnv(ref, alt)) {
            return getDannPathogenicityData(dannSnvTabixDataSource, chromosome, start, ref, alt);
        }
        return getDannPathogenicityData(dannInDelTabixDataSource, chromosome, start, ref, alt);
    }

    private PathogenicityData getDannPathogenicityData(TabixDataSource tabixDataSource, String chromosome, int start, String ref, String alt) {
        try {
            TabixReader.Iterator results = tabixDataSource.query(chromosome + ":" + start + "-" + start);
            String line;
            //there can be 0 - N results in this format:
            //1       10001   T       A       0.16461391399220135
            //1       10001   T       C       0.4396994049749739
            //1       10001   T       G       0.38108629377072734
            while ((line = results.next()) != null) {
                String[] elements = line.split("\t");
                String dannRef = elements[2];
                String dannAlt = elements[3];
                if (dannRef.equals(ref) && dannAlt.equals(alt)) {
                    return makeDannPathData(elements[4]);
                }
            }
        } catch (IOException e) {
            logger.error("Unable to read from DANN tabix file {}", tabixDataSource.getSource(), e);
        }
        return PathogenicityData.empty();
    }
 
    private PathogenicityData makeDannPathData(String phredScaledDannScore) {
        float score = Float.parseFloat(phredScaledDannScore);
        DannScore dannScore = DannScore.of(score);
        return PathogenicityData.of(dannScore);
    }
}
