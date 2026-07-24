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

import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencyData;
import org.monarchinitiative.svart.CoordinateSystem;
import org.monarchinitiative.svart.Coordinates;
import org.monarchinitiative.svart.Strand;

import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
class ReportLinksTest {

    private VariantEvaluation variant(String rsId) {
        VariantEvaluation.Builder builder = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(10), Strand.POSITIVE, CoordinateSystem.oneBased(), 123256215, "T", "G")
                .geneSymbol("FGFR2");
        if (rsId != null) {
            builder.frequencyData(FrequencyData.of(rsId));
        }
        return builder.build();
    }

    private VariantEvaluation longDel() {
        // Sequence-resolved DEL: REF length must match coordinate span (svart validates this).
        String longRef = "A" + "T".repeat(100); // 101 bp
        int start = 52772578;
        int end = start + longRef.length() - 1;
        return VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(12), Strand.POSITIVE,
                        Coordinates.oneBased(start, end), longRef, "A", 1 - longRef.length())
                .geneSymbol("KRT2")
                .build();
    }

    @Test
    void interpolateReplacesPlaceholders() {
        String result = ReportLinks.interpolate("x/{assembly}/{chr}:{pos}:{ref}:{alt}",
                Map.of("assembly", "hg19", "chr", "10", "pos", "123", "ref", "T", "alt", "G"));
        assertEquals("x/hg19/10:123:T:G", result);
    }

    @Test
    void interpolateReturnsNullWhenPlaceholderMissing() {
        assertNull(ReportLinks.interpolate("x/Rs{rsId}", Map.of("rsId", "")));
        assertNull(ReportLinks.interpolate("x/{gene}", Map.of()));
    }

    @Test
    void rsIdNumberStripsRsPrefix() {
        assertEquals("141235720", ReportLinks.rsIdNumber("rs141235720"));
        assertEquals("141235720", ReportLinks.rsIdNumber("RS141235720"));
        assertEquals("141235720", ReportLinks.rsIdNumber("141235720"));
        assertEquals("", ReportLinks.rsIdNumber(""));
        assertEquals("", ReportLinks.rsIdNumber(null));
    }

    @Test
    void ensemblHostByAssembly() {
        assertEquals("grch37", ReportLinks.ensemblHost(GenomeAssembly.HG19));
        assertEquals("www", ReportLinks.ensemblHost(GenomeAssembly.HG38));
    }

    @Test
    void ucscChromAddsPrefix() {
        assertEquals("chr12", ReportLinks.ucscChrom("12"));
        assertEquals("chrX", ReportLinks.ucscChrom("X"));
        assertEquals("chrM", ReportLinks.ucscChrom("MT"));
        assertEquals("chr12", ReportLinks.ucscChrom("chr12"));
    }

    @Test
    void alleleForUrlRejectsLong() {
        assertEquals("T", ReportLinks.alleleForUrl("T"));
        assertEquals("", ReportLinks.alleleForUrl("A" + "C".repeat(50)));
        assertEquals("", ReportLinks.alleleForUrl(""));
        assertEquals("", ReportLinks.alleleForUrl(null));
    }

    @Test
    void resolveOmitsVarSomeForLongAllelesButEmitsRegionLink() {
        ReportLinks links = new ReportLinks(List.of(
                new ReportLink("VarSome", "https://varsome.com/variant/{assembly}/{chr}:{pos}:{ref}:{alt}"),
                new ReportLink("Ensembl", "https://{ensemblHost}.ensembl.org/Homo_sapiens/Location/View?r={chr}:{pos}-{end}")
        ));
        List<ReportLinks.ResolvedLink> resolved = links.resolve(longDel());
        VariantEvaluation del = longDel();
        assertEquals(1, resolved.size());
        assertEquals("Ensembl", resolved.get(0).getLabel());
        assertEquals(
                "https://grch37.ensembl.org/Homo_sapiens/Location/View?r=12:" + del.start() + "-" + del.end(),
                resolved.get(0).getUrl());
    }

    @Test
    void resolveDefaultLinksIncludeRegionForLongDel() {
        VariantEvaluation del = longDel();
        List<ReportLinks.ResolvedLink> resolved = ReportLinks.loadDefault().resolve(del);
        assertTrue(resolved.stream().noneMatch(l -> l.getLabel().equals("VarSome")));
        assertTrue(resolved.stream().anyMatch(l ->
                l.getLabel().equals("Ensembl")
                        && l.getUrl().equals("https://grch37.ensembl.org/Homo_sapiens/Location/View?r=12:"
                        + del.start() + "-" + del.end())));
    }

    @Test
    void resolveUcscRegionTemplate() {
        VariantEvaluation del = longDel();
        ReportLinks links = new ReportLinks(List.of(
                new ReportLink("UCSC", "https://genome.ucsc.edu/cgi-bin/hgTracks?db={ucscDb}&position={chrUcsc}:{pos}-{end}")
        ));
        List<ReportLinks.ResolvedLink> resolved = links.resolve(del);
        assertEquals(1, resolved.size());
        assertEquals(
                "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr12:" + del.start() + "-" + del.end(),
                resolved.get(0).getUrl());
    }

    @Test
    void loadMissingResourceGivesEmptyLinks() {
        assertTrue(ReportLinks.load("/no-such-report-links.yml").getLinks().isEmpty());
    }

    @Test
    void loadDefaultLoadsBundledLinks() {
        List<ReportLink> links = ReportLinks.loadDefault().getLinks();
        assertTrue(links.stream().anyMatch(link -> link.label().equals("VarSome")));
        assertTrue(links.stream().anyMatch(link -> link.label().equals("Ensembl")));
        assertTrue(links.stream().anyMatch(link -> link.label().equals("SNPedia")));
    }

    @Test
    void resolveFillsAllPlaceholdersForVariantWithRsId() {
        ReportLinks reportLinks = new ReportLinks(List.of(
                new ReportLink("VarSome", "https://varsome.com/variant/{assembly}/{chr}:{pos}:{ref}:{alt}"),
                new ReportLink("Ensembl", "https://{ensemblHost}.ensembl.org/Homo_sapiens/Location/View?r={chr}:{pos}-{end}"),
                new ReportLink("SNPedia", "https://www.snpedia.com/index.php/Rs{rsId}")
        ));

        List<ReportLinks.ResolvedLink> resolved = reportLinks.resolve(variant("rs141235720"));

        assertEquals(3, resolved.size());
        assertEquals("VarSome", resolved.get(0).getLabel());
        assertEquals("https://varsome.com/variant/hg19/10:123256215:T:G", resolved.get(0).getUrl());
        assertEquals("Ensembl", resolved.get(1).getLabel());
        assertEquals(
                "https://grch37.ensembl.org/Homo_sapiens/Location/View?r=10:123256215-123256215",
                resolved.get(1).getUrl());
        assertEquals("SNPedia", resolved.get(2).getLabel());
        assertEquals("https://www.snpedia.com/index.php/Rs141235720", resolved.get(2).getUrl());
    }

    @Test
    void resolveOmitsLinksWithUnresolvablePlaceholders() {
        ReportLinks reportLinks = new ReportLinks(List.of(
                new ReportLink("VarSome", "https://varsome.com/variant/{assembly}/{chr}:{pos}:{ref}:{alt}"),
                new ReportLink("SNPedia", "https://www.snpedia.com/index.php/Rs{rsId}")
        ));

        // no rsId -> only the VarSome link (which needs no rsId) is emitted
        List<ReportLinks.ResolvedLink> resolved = reportLinks.resolve(variant(null));

        assertEquals(1, resolved.size());
        assertEquals("VarSome", resolved.get(0).getLabel());
    }
}
