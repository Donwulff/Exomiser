/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2017 Queen Mary University of London.
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

package org.monarchinitiative.exomiser.core.genome.dao.serialisers;

import org.h2.mvstore.MVMap;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.proto.AlleleProto.AlleleKey;
import org.monarchinitiative.exomiser.core.proto.AlleleProto.AlleleProperties;

/**
 * Utility class for helping read and write Alleles to the {@link org.h2.mvstore.MVStore}
 *
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public class MvStoreUtil {

    private MvStoreUtil() {
        //static utility class - not instantiable
    }

    public static MVMap.Builder<AlleleKey, AlleleProperties> alleleMapBuilder() {
        return new MVMap.Builder<AlleleKey, AlleleProperties>()
                .keyType(AlleleKeyDataType.INSTANCE)
                .valueType(AllelePropertiesDataType.INSTANCE);
    }

    public static AlleleKey generateAlleleKey(Variant variant) {
        return AlleleKey.newBuilder()
                .setChr(variant.getChromosome())
                .setPosition(variant.getPosition())
                .setRef(variant.getRef())
                .setAlt(variant.getAlt())
                .build();
    }
}
