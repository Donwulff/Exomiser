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
import com.fasterxml.jackson.dataformat.yaml.YAMLMapper;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A configurable set of external-information links for the HTML report. Links are defined as a label
 * plus a URL template (see {@code report-links.yml} on the classpath); the template placeholders are
 * filled in per-variant by {@link #resolve(VariantEvaluation)}.
 * <p>
 * Supported placeholders: {@code {assembly}} (hg19/hg38), {@code {chr}}, {@code {pos}}, {@code {ref}},
 * {@code {alt}}, {@code {gene}}, {@code {rsId}} (the numeric part of the dbSNP id). A link is only
 * emitted for a variant when every placeholder it references has a value for that variant, so e.g. an
 * rsId-based link is silently omitted for variants with no rsId.
 * <p>
 * Note: outbound links disclose which variants are being investigated to third-party sites.
 *
 * @since 15.1.0
 */
public class ReportLinks {

    private static final Logger logger = LoggerFactory.getLogger(ReportLinks.class);

    private static final String DEFAULT_RESOURCE = "/report-links.yml";
    private static final YAMLMapper YAML_MAPPER = new YAMLMapper();
    private static final Pattern PLACEHOLDER = Pattern.compile("\\{([a-zA-Z][a-zA-Z0-9]*)}");

    private final List<ReportLink> links;

    public ReportLinks(List<ReportLink> links) {
        this.links = links == null ? List.of() : List.copyOf(links);
    }

    /**
     * Loads the default link definitions bundled on the classpath ({@code /report-links.yml}). Any
     * problem loading the resource is logged and results in an empty (no-op) set of links rather than
     * a failure, so a malformed override can never break report generation.
     */
    public static ReportLinks loadDefault() {
        return load(DEFAULT_RESOURCE);
    }

    static ReportLinks load(String resourcePath) {
        try (InputStream inputStream = ReportLinks.class.getResourceAsStream(resourcePath)) {
            if (inputStream == null) {
                logger.warn("No report links resource found at {} - external report links disabled", resourcePath);
                return new ReportLinks(List.of());
            }
            Config config = YAML_MAPPER.readValue(inputStream, Config.class);
            return new ReportLinks(config.links());
        } catch (IOException e) {
            logger.error("Error reading report links resource {} - external report links disabled", resourcePath, e);
            return new ReportLinks(List.of());
        }
    }

    public List<ReportLink> getLinks() {
        return links;
    }

    /**
     * @return the configured links resolved for the given variant, omitting any link whose template
     * references a placeholder that has no value for this variant.
     */
    public List<ResolvedLink> resolve(VariantEvaluation variant) {
        if (links.isEmpty()) {
            return List.of();
        }
        Map<String, String> values = placeholderValues(variant);
        List<ResolvedLink> resolved = new ArrayList<>(links.size());
        for (ReportLink link : links) {
            String url = interpolate(link.url(), values);
            if (url != null) {
                resolved.add(new ResolvedLink(link.label(), url));
            }
        }
        return resolved;
    }

    private static Map<String, String> placeholderValues(VariantEvaluation variant) {
        Map<String, String> values = new HashMap<>();
        values.put("assembly", variant.genomeAssembly().toString());
        values.put("chr", variant.contigName());
        values.put("pos", Integer.toString(variant.start()));
        values.put("ref", variant.ref());
        values.put("alt", variant.alt());
        values.put("gene", blankIfPlaceholder(variant.geneSymbol()));
        values.put("rsId", rsIdNumber(variant.frequencyData().getRsId()));
        return values;
    }

    /**
     * Substitutes {@code {placeholder}} tokens in the template. Returns {@code null} if the template
     * references a placeholder with no value, signalling that the link should be skipped for this
     * variant.
     */
    static String interpolate(String template, Map<String, String> values) {
        Matcher matcher = PLACEHOLDER.matcher(template);
        StringBuilder result = new StringBuilder();
        while (matcher.find()) {
            String value = values.get(matcher.group(1));
            if (value == null || value.isEmpty()) {
                return null;
            }
            matcher.appendReplacement(result, Matcher.quoteReplacement(value));
        }
        matcher.appendTail(result);
        return result.toString();
    }

    /**
     * @return the numeric part of a dbSNP id (a leading, case-insensitive "rs" is stripped) so that
     * templates can freely prefix it, e.g. {@code Rs{rsId}} for SNPedia or {@code rs{rsId}} for dbSNP.
     */
    static String rsIdNumber(String rsId) {
        if (rsId == null) {
            return "";
        }
        String trimmed = rsId.trim();
        if (trimmed.regionMatches(true, 0, "rs", 0, 2)) {
            trimmed = trimmed.substring(2);
        }
        return trimmed;
    }

    private static String blankIfPlaceholder(String geneSymbol) {
        return geneSymbol == null || geneSymbol.equals(".") ? "" : geneSymbol;
    }

    private record Config(@JsonProperty("links") List<ReportLink> links) {
        private Config {
            links = links == null ? List.of() : links;
        }
    }

    /**
     * A link whose placeholders have been resolved for a specific variant. Exposes JavaBean getters
     * for use from the OGNL-based Thymeleaf report template.
     */
    public static final class ResolvedLink {

        private final String label;
        private final String url;

        ResolvedLink(String label, String url) {
            this.label = label;
            this.url = url;
        }

        public String getLabel() {
            return label;
        }

        public String getUrl() {
            return url;
        }
    }
}
