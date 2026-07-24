/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2022 Queen Mary University of London.
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

import java.util.List;

/**
 * Reads CADD PHRED scores from tabix file(s). All parsing is handled by {@link TabixPathogenicityDao}; this class only
 * supplies the {@link TabixPathogenicitySource#cadd()} preset and its own cache name.
 * <p>
 * Multiple SNV and/or indel tabix files may be supplied; for each variant the highest matching PHRED score across the
 * applicable files is returned as a single {@link org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource#CADD}
 * score (official + novel sample scores under one identity).
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class CaddDao extends TabixPathogenicityDao {

    public CaddDao(TabixDataSource caddInDelTabixDataSource, TabixDataSource caddSnvTabixDataSource) {
        this(List.of(caddInDelTabixDataSource), List.of(caddSnvTabixDataSource));
    }

    /**
     * @param caddInDelTabixDataSources indel tabix sources (order is not significant; max score wins)
     * @param caddSnvTabixDataSources   SNV tabix sources (order is not significant; max score wins)
     */
    public CaddDao(List<TabixDataSource> caddInDelTabixDataSources, List<TabixDataSource> caddSnvTabixDataSources) {
        super(TabixPathogenicitySource.cadd(), caddSnvTabixDataSources, caddInDelTabixDataSources);
    }

    @Caching(cacheable = {
            @Cacheable(cacheNames = "hg19.cadd", keyGenerator = "variantKeyGenerator", condition = "#variant.genomeAssembly == T(org.monarchinitiative.exomiser.core.genome.GenomeAssembly).HG19"),
            @Cacheable(cacheNames = "hg38.cadd", keyGenerator = "variantKeyGenerator", condition = "#variant.genomeAssembly == T(org.monarchinitiative.exomiser.core.genome.GenomeAssembly).HG38"),
    })
    @Override
    public PathogenicityData getPathogenicityData(Variant variant) {
        return super.getPathogenicityData(variant);
    }
}
