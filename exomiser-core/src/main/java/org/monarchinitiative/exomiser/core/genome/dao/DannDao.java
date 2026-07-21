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
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.springframework.cache.annotation.Cacheable;
import org.springframework.cache.annotation.Caching;

/**
 * Reads DANN scores (CADD-style, SNV + indel, already scaled 0-1) from tabix files. All parsing is handled by
 * {@link TabixPathogenicityDao}; this class only supplies the {@link TabixPathogenicitySource#dann()} preset and its own
 * cache name.
 */
public class DannDao extends TabixPathogenicityDao {

    public DannDao(TabixDataSource dannInDelTabixDataSource, TabixDataSource dannSnvTabixDataSource) {
        super(TabixPathogenicitySource.dann(), dannSnvTabixDataSource, dannInDelTabixDataSource);
    }

    @Caching(cacheable = {
            @Cacheable(cacheNames = "hg19.dann", keyGenerator = "variantKeyGenerator", condition = "#variant.genomeAssembly == T(org.monarchinitiative.exomiser.core.genome.GenomeAssembly).HG19"),
            @Cacheable(cacheNames = "hg38.dann", keyGenerator = "variantKeyGenerator", condition = "#variant.genomeAssembly == T(org.monarchinitiative.exomiser.core.genome.GenomeAssembly).HG38"),
    })
    @Override
    public PathogenicityData getPathogenicityData(Variant variant) {
        return super.getPathogenicityData(variant);
    }
}
