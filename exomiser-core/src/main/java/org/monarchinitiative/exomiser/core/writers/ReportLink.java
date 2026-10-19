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

package org.monarchinitiative.exomiser.core.writers;

import com.fasterxml.jackson.annotation.JsonProperty;

/**
 * Definition of a single external-information link shown next to a variant in the HTML report. The
 * {@code url} is a template that may contain placeholders resolved per-variant by {@link ReportLinks}.
 * <p>
 * {@code key} marks a link as important enough for space-constrained views (the web UI report renders
 * only key links; the full HTML report renders all of them).
 *
 * @since 15.1.0
 */
public record ReportLink(@JsonProperty("label") String label, @JsonProperty("url") String url,
                         @JsonProperty("key") boolean key) {

    public ReportLink(String label, String url) {
        this(label, url, false);
    }
}
