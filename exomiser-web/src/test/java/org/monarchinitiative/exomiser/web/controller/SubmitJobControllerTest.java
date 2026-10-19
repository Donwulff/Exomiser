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
package org.monarchinitiative.exomiser.web.controller;

import config.TestConfig;
import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.jsoup.Jsoup;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.monarchinitiative.exomiser.core.Exomiser;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.analysis.acmg.AcmgAssignment;
import org.monarchinitiative.exomiser.core.analysis.acmg.AcmgClassification;
import org.monarchinitiative.exomiser.core.analysis.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.analysis.sample.Sample;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.model.Gene;
import org.monarchinitiative.exomiser.core.model.GeneScore;
import org.monarchinitiative.exomiser.core.model.SampleGenotype;
import org.monarchinitiative.exomiser.core.model.SampleGenotypes;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.prioritisers.model.Disease;
import org.monarchinitiative.exomiser.test.ExomiserStubDataConfig;
import org.monarchinitiative.exomiser.web.ExomiserWebApp;
import org.monarchinitiative.svart.Coordinates;
import org.monarchinitiative.svart.Strand;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.AutoConfigureMockMvc;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.mock.web.MockMultipartFile;
import org.springframework.test.context.bean.override.mockito.MockitoSpyBean;
import org.springframework.test.context.junit.jupiter.SpringExtension;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.request.MockMvcRequestBuilders;

import java.nio.charset.StandardCharsets;
import java.util.EnumSet;
import java.util.List;

import static org.hamcrest.Matchers.containsString;
import static org.hamcrest.Matchers.not;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.doReturn;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.model;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.view;

/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
@ExtendWith(SpringExtension.class)
@AutoConfigureMockMvc
@SpringBootTest(classes = {ExomiserWebApp.class, ExomiserStubDataConfig.class, TestConfig.class})
public class SubmitJobControllerTest {

    @Autowired
    private MockMvc mockMvc;

    @MockitoSpyBean
    private Exomiser exomiser;

    @Test
    public void submitJobControllerGetReturnsOkStatusAndSubmitPage() throws Exception {
        mockMvc.perform(MockMvcRequestBuilders.get("/submit"))
                .andExpect(status().isOk())
                .andExpect(view().name("submit"));
    }

    @Test
    void submittedResultsRenderKeyLinksAndCompactVariants() throws Exception {
        VariantEvaluation snv = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(10), Strand.POSITIVE,
                        Coordinates.oneBased(123256215, 123256215), "T", "G", 0)
                .geneSymbol("FGFR2").sampleGenotypes(SampleGenotypes.of("sample", SampleGenotype.het())).build();
        snv.setContributesToGeneScoreUnderMode(ModeOfInheritance.AUTOSOMAL_DOMINANT);
        String longRef = "A" + "T".repeat(100);
        VariantEvaluation deletion = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(10), Strand.POSITIVE,
                        Coordinates.oneBased(123256220, 123256320), longRef, "A", -100)
                .geneSymbol("FGFR2").sampleGenotypes(SampleGenotypes.of("sample", SampleGenotype.het())).build();
        Gene gene = new Gene("FGFR2", 2263);
        gene.addVariant(snv);
        gene.addVariant(deletion);
        gene.setCompatibleInheritanceModes(EnumSet.of(ModeOfInheritance.AUTOSOMAL_DOMINANT));
        AcmgAssignment assignment = AcmgAssignment.of(snv, gene.geneIdentifier(), ModeOfInheritance.AUTOSOMAL_DOMINANT,
                Disease.builder().build(), AcmgEvidence.empty(), AcmgClassification.LIKELY_PATHOGENIC);
        gene.addGeneScore(new GeneScore(gene.geneIdentifier(), ModeOfInheritance.AUTOSOMAL_DOMINANT,
                1, 1, 1, 0, List.of(snv), List.of(), List.of(assignment)));
        var results = AnalysisResults.builder().sample(Sample.builder().build()).analysis(Analysis.builder().build())
                .sampleNames(List.of("sample")).genes(List.of(gene)).variantEvaluations(List.of(snv, deletion)).build();
        doReturn(results).when(exomiser).run(any(Sample.class), any(Analysis.class));

        String vcf = """
                ##fileformat=VCFv4.2
                #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
                10	123256215	.	T	G	.	PASS	.
                """;
        var upload = new MockMultipartFile("vcf", "synthetic.vcf", "text/plain", vcf.getBytes(StandardCharsets.UTF_8));
        var response = mockMvc.perform(MockMvcRequestBuilders.multipart("/submit").file(upload)
                        .param("phenotypes", "HP:0000001").param("frequency", "1")
                        .param("remove-dbsnp", "false").param("keep-non-pathogenic", "true")
                        .param("keep-off-target", "true").param("inheritance", "ANY")
                        .param("prioritiser", "OMIM_PRIORITY"))
                .andExpect(status().isOk())
                .andExpect(view().name("results"))
                .andExpect(model().attributeExists("reportLinks"))
                .andExpect(content().string(containsString(snv.toDisplayString())))
                .andExpect(content().string(containsString(deletion.toDisplayString())))
                .andExpect(content().string(not(containsString(longRef))))
                .andReturn().getResponse();

        var links = Jsoup.parse(response.getContentAsString()).select("a.external-link");
        assertEquals(List.of("VarSome", "gnomAD"), links.eachText());
        assertEquals("https://gnomad.broadinstitute.org/variant/10-123256215-T-G?dataset=gnomad_r2_1",
                links.get(1).attr("href"));
    }

}
