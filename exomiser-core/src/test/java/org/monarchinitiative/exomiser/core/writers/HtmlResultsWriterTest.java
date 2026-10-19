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

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.monarchinitiative.exomiser.core.writers;

import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.jsoup.Jsoup;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import org.monarchinitiative.exomiser.api.v1.OutputProto;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.analysis.acmg.AcmgAssignment;
import org.monarchinitiative.exomiser.core.analysis.acmg.AcmgClassification;
import org.monarchinitiative.exomiser.core.analysis.acmg.AcmgEvidence;
import org.monarchinitiative.exomiser.core.analysis.sample.Sample;
import org.monarchinitiative.exomiser.core.filters.*;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.genome.TestFactory;
import org.monarchinitiative.exomiser.core.genome.TestVariantFactory;
import org.monarchinitiative.exomiser.core.model.Gene;
import org.monarchinitiative.exomiser.core.model.GeneScore;
import org.monarchinitiative.exomiser.core.model.SampleGenotype;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.frequency.Frequency;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencyData;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.monarchinitiative.exomiser.core.model.pathogenicity.*;
import org.monarchinitiative.exomiser.core.prioritisers.OmimPriorityResult;
import org.monarchinitiative.exomiser.core.prioritisers.PhivePriority;
import org.monarchinitiative.exomiser.core.prioritisers.model.Disease;
import org.monarchinitiative.exomiser.core.prioritisers.service.TestPriorityServiceFactory;
import org.monarchinitiative.exomiser.core.proto.ProtoParser;
import org.monarchinitiative.svart.Coordinates;
import org.monarchinitiative.svart.Strand;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.EnumSet;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class HtmlResultsWriterTest {

    private final HtmlResultsWriter instance = new HtmlResultsWriter();

    private static Path testOutDir;

    private VariantEvaluation unAnnotatedVariantEvaluation1;
    private VariantEvaluation unAnnotatedVariantEvaluation2;

    private Gene fgfr2Gene;
    private Gene shhGene;

    @BeforeAll
    public static void makeTempDir() throws IOException {
        testOutDir = Files.createTempDirectory("exomiser_html_writer_test");
    }

    @AfterAll
    public static void tearDown() throws IOException {
        Files.delete(testOutDir);
    }

    @BeforeEach
    public void setUp() {
        TestVariantFactory varFactory = new TestVariantFactory();

        VariantEvaluation fgfr2MissenseVariantEvaluation = TestVariantFactory.buildVariant(10, 123256215, "T", "G", SampleGenotype.het(), 30, 2.2);
        fgfr2MissenseVariantEvaluation.setFrequencyData(FrequencyData.of("rs123456", Frequency.of(FrequencySource.THOUSAND_GENOMES, 0.01f)));
        fgfr2MissenseVariantEvaluation.setPathogenicityData(PathogenicityData.of(PolyPhenScore.of(1f), MutationTasterScore
                .of(1f), SiftScore.of(0f), CaddScore.of(1f)));
        fgfr2MissenseVariantEvaluation.addFilterResult(FilterResult.pass(FilterType.FREQUENCY_FILTER));
        fgfr2MissenseVariantEvaluation.addFilterResult(FilterResult.pass(FilterType.VARIANT_EFFECT_FILTER));

        //TODO: make a few canned variants from the smallTest.vcf and Pfeiffer.vcf these can then be used to run the full system without the vcf if required.
        fgfr2Gene = TestFactory.newGeneFGFR2();
        fgfr2Gene.addVariant(fgfr2MissenseVariantEvaluation);

        VariantEvaluation shhIndelVariantEvaluation = TestVariantFactory.buildVariant(7, 155604800, "C", "CTT", SampleGenotype.het(), 30, 1.0);
        shhGene = TestFactory.newGeneSHH();
        shhGene.addVariant(shhIndelVariantEvaluation);

        fgfr2Gene.addPriorityResult(new OmimPriorityResult(fgfr2Gene.entrezGeneId(), fgfr2Gene.geneSymbol(), 1f, Collections.emptyList(), Collections.emptyMap()));
        fgfr2Gene.addGeneScore(new GeneScore(fgfr2Gene.geneIdentifier(), ModeOfInheritance.AUTOSOMAL_DOMINANT, 1.0, 1.0, 1.0, 0.0, List.of(fgfr2MissenseVariantEvaluation), List.of(), List.of()));
        shhGene.addPriorityResult(new OmimPriorityResult(shhGene.entrezGeneId(), shhGene.geneSymbol(), 1f, Collections.emptyList(), Collections.emptyMap()));

        unAnnotatedVariantEvaluation1 = TestVariantFactory.buildVariant(5, 10, "C", "T", SampleGenotype.het(), 30, 1.0);
        unAnnotatedVariantEvaluation2 = TestVariantFactory.buildVariant(5, 10, "C", "T", SampleGenotype.het(), 30, 1.0);
    }

    private AnalysisResults buildAnalysisResults(Sample sample, Analysis analysis, List<Gene> genes, List<VariantEvaluation> variantEvaluations) {
        return AnalysisResults.builder()
                .sample(sample)
                .analysis(analysis)
                .sampleNames(List.of("Slartibartfast"))
                .genes(genes)
                .variantEvaluations(variantEvaluations)
                .build();
    }

    @Test
    public void testWriteTemplateWithEmptyData() {
        Sample sample = Sample.builder().build();
        Analysis analysis = Analysis.builder().build();
        AnalysisResults analysisResults = buildAnalysisResults(sample, analysis, Collections.emptyList(), Collections.emptyList());

        String outputPrefix = testOutDir.resolve("testWriteTemplateWithEmptyData").toString();
        OutputSettings settings = OutputSettings.builder()
                .outputPrefix(outputPrefix)
                .build();

        instance.writeFile(analysisResults, settings);
        Path testOutFile = Paths.get(outputPrefix + ".html");
        assertTrue(testOutFile.toFile().exists());
        assertTrue(testOutFile.toFile().delete());
    }

    @ParameterizedTest
    @EnumSource(value = OutputFormat.class, names = {"HTML", "HTML_LEGACY"})
    void bothTemplatesShowUnannotatedVariantsOutsideGeneResults(OutputFormat format) {
        AnalysisResults results = buildAnalysisResults(Sample.builder().build(), Analysis.builder().build(),
                List.of(), List.of(unAnnotatedVariantEvaluation1));
        assertEquals(1, results.unAnnotatedVariantEvaluations().size());

        String output = ResultsWriterFactory.getResultsWriter(format).writeString(results, OutputSettings.builder().build());

        assertTrue(output.contains("Unanalysed Variants"));
        assertTrue(output.contains(unAnnotatedVariantEvaluation1.toDisplayString()));
    }

    @Test
    public void testWriteTemplateWithUnAnnotatedVariantData() {
        Sample sample = Sample.builder().build();
        Analysis analysis = Analysis.builder().build();

        List<VariantEvaluation> variants = List.of(unAnnotatedVariantEvaluation1, unAnnotatedVariantEvaluation2);
        AnalysisResults analysisResults = buildAnalysisResults(sample, analysis, Collections.emptyList(), variants);

        String testOutFilePrefix = testOutDir.resolve("testWriteTemplateWithUnAnnotatedVariantData").toString();
        OutputSettings settings = OutputSettings.builder()
                .outputPrefix(testOutFilePrefix)
                .build();

        instance.writeFile(analysisResults, settings);

        Path testOutFile = Paths.get(testOutFilePrefix + ".html");
        assertTrue(testOutFile.toFile().exists());
        assertTrue(testOutFile.toFile().delete());
    }

    @Test
    public void testWriteTemplateWithUnAnnotatedVariantDataAndGenes() throws Exception {
        Sample sample = Sample.builder().build();
        Analysis analysis = Analysis.builder().build();

        List<VariantEvaluation> variants = List.of(unAnnotatedVariantEvaluation1, unAnnotatedVariantEvaluation2);
        List<Gene> genes = List.of(fgfr2Gene, shhGene);
        AnalysisResults analysisResults = buildAnalysisResults(sample, analysis, genes, variants);

        String testOutFilePrefix = testOutDir.resolve("testWriteTemplateWithUnAnnotatedVariantDataAndGenes").toString();
        OutputSettings settings = OutputSettings.builder().outputPrefix(testOutFilePrefix).build();

        instance.writeFile(analysisResults, settings);
        Path testOutFile = Paths.get(testOutFilePrefix + ".html");

        List<String> lines = Files.readAllLines(testOutFile);
        assertFalse(lines.isEmpty());

        assertTrue(testOutFile.toFile().exists());
        assertTrue(testOutFile.toFile().delete());
    }

    @Test
    public void testWriteTemplateWithEmptyDataAndFullAnalysis() throws Exception {

        Sample sample = Sample.builder()
                .hpoIds(List.of("HP:000001", "HP:000002"))
                .build();

        Analysis analysis = Analysis.builder()
                .addStep(new RegulatoryFeatureFilter())
                .addStep(new FrequencyFilter(0.1f))
                .addStep(new PathogenicityFilter(true))
                .addStep(new PhivePriority(TestPriorityServiceFactory.stubPriorityService()))
                .build();

        AnalysisResults analysisResults = buildAnalysisResults(sample, analysis, Collections.emptyList(), Collections.emptyList());
        OutputSettings settings = OutputSettings.builder().build();

        String output = instance.writeString(analysisResults, settings);
        assertFalse(output.isEmpty());
    }

    @Test
    public void testWriteTemplateWithFullDataAndFullAnalysis() throws Exception {

        unAnnotatedVariantEvaluation1.setContributesToGeneScoreUnderMode(ModeOfInheritance.ANY);

        List<VariantEvaluation> variants = List.of(unAnnotatedVariantEvaluation1, unAnnotatedVariantEvaluation2);
        List<Gene> genes = List.of(fgfr2Gene, shhGene);

        Sample sample = Sample.builder()
                .hpoIds(List.of("HP:000001", "HP:000002"))
                .build();

        Analysis analysis = Analysis.builder()
                .addStep(new RegulatoryFeatureFilter())
                .addStep(new FrequencyFilter(0.1f))
                .addStep(new PathogenicityFilter(true))
                .addStep(new PhivePriority(TestPriorityServiceFactory.stubPriorityService()))
                .build();
        AnalysisResults analysisResults = buildAnalysisResults(sample, analysis, genes, variants);

        OutputSettings settings = OutputSettings.builder().build();

        String output = instance.writeString(analysisResults, settings);
        assertTrue(output.contains("Gene Symbol"));
        assertTrue(output.contains("FGFR2"));
    }

    @Test
    void bothHtmlFormatsSurviveYamlParsingAndWriteSeparateReports(@TempDir Path outputDir) throws IOException {
        var proto = ProtoParser.parseFromJsonOrYaml(OutputProto.OutputOptions.newBuilder(), """
                outputFileName: report
                outputFormats: [HTML, HTML_LEGACY]
                """).build();
        var settings = new OutputSettingsProtoConverter().toDomain(proto.toBuilder()
                .setOutputDirectory(outputDir.toString()).build());
        AnalysisResults results = buildAnalysisResults(Sample.builder().build(), Analysis.builder().build(),
                List.of(fgfr2Gene, shhGene), List.of());

        AnalysisResultsWriter.writeToFile(results, settings);

        String current = Files.readString(outputDir.resolve("report.html"));
        String legacy = Files.readString(outputDir.resolve("report.legacy.html"));
        assertEquals(instance.writeString(results, settings), current);
        assertEquals(HtmlResultsWriter.legacy().writeString(results, settings), legacy);
        assertFalse(current.equals(legacy));
    }

    @ParameterizedTest
    @EnumSource(value = OutputFormat.class, names = {"HTML", "HTML_LEGACY"})
    void bothTemplatesRenderContributingVariantWithAcmgAssignment(OutputFormat format) {
        VariantEvaluation contributingVariant = TestVariantFactory.buildVariant(10, 123256215, "T", "G", SampleGenotype.het(), 30, 2.2);
        var context = contributingVariant.variantContext();
        var filteredContext = new VariantContextBuilder(context).filter("VQLOW")
                .genotypes(new GenotypeBuilder(context.getGenotype(0)).GQ(78).DP(42).make()).make();
        contributingVariant = contributingVariant.toBuilder().variantContext(filteredContext).quality(-10).build();
        contributingVariant.setFrequencyData(FrequencyData.of("rs123456", Frequency.of(FrequencySource.THOUSAND_GENOMES, 0.01f)));
        contributingVariant.setPathogenicityData(PathogenicityData.of(PolyPhenScore.of(1f)));
        contributingVariant.setContributesToGeneScoreUnderMode(ModeOfInheritance.AUTOSOMAL_DOMINANT);

        Gene gene = TestFactory.newGeneFGFR2();
        gene.addVariant(contributingVariant);
        gene.setCompatibleInheritanceModes(EnumSet.of(ModeOfInheritance.AUTOSOMAL_DOMINANT));
        AcmgAssignment acmgAssignment = AcmgAssignment.of(contributingVariant, gene.geneIdentifier(),
                ModeOfInheritance.AUTOSOMAL_DOMINANT, Disease.builder().build(), AcmgEvidence.empty(),
                AcmgClassification.LIKELY_PATHOGENIC);
        gene.addGeneScore(new GeneScore(gene.geneIdentifier(), ModeOfInheritance.AUTOSOMAL_DOMINANT, 1.0, 1.0, 1.0, 0.0,
                List.of(contributingVariant), List.of(), List.of(acmgAssignment)));

        AnalysisResults analysisResults = buildAnalysisResults(Sample.builder().build(), Analysis.builder().build(),
                List.of(gene), List.of(contributingVariant));

        String output = ResultsWriterFactory.getResultsWriter(format).writeString(analysisResults, OutputSettings.builder().build());

        // the contributing-variant section is driven by the ACMG assignments - check it rendered with the
        // display string, report links and caller call-quality line
        assertTrue(output.contains(format == OutputFormat.HTML ? "LIKELY PATHOGENIC" : "LIKELY_PATHOGENIC"));
        assertTrue(output.contains("10-123256215-T-G"));
        assertTrue(output.contains("varsome.com"));
        assertTrue(output.contains("VQLOW"));
        String renderedText = Jsoup.parse(output).text();
        assertTrue(renderedText.contains("GQ 78"));
        assertTrue(renderedText.contains("DP 42"));
        assertFalse(output.matches("(?s).*?>\\s*-10(?:\\.0+)?\\s*<.*"));
    }

    @ParameterizedTest
    @EnumSource(value = OutputFormat.class, names = {"HTML", "HTML_LEGACY"})
    void bothTemplatesRenderLongVariantsWithoutAlleleLinks(OutputFormat format) {
        String ref = "A" + "T".repeat(100);
        VariantEvaluation deletion = VariantEvaluation.builder()
                .variant(GenomeAssembly.HG19.getContigById(10), Strand.POSITIVE,
                        Coordinates.oneBased(123256215, 123256315), ref, "A", -100)
                .geneSymbol("FGFR2")
                .build();
        Gene gene = TestFactory.newGeneFGFR2();
        gene.addVariant(deletion);
        gene.addGeneScore(GeneScore.builder().geneIdentifier(gene.geneIdentifier()).combinedScore(1).build());
        AnalysisResults results = buildAnalysisResults(Sample.builder().build(), Analysis.builder().build(),
                List.of(gene), List.of(deletion));

        String output = ResultsWriterFactory.getResultsWriter(format).writeString(results, OutputSettings.builder().build());

        assertTrue(output.contains(deletion.toDisplayString()));
        assertTrue(output.contains("grch37.ensembl.org/Homo_sapiens/Location/View?r=10:123256215-123256315"));
        assertFalse(output.contains(ref));
        assertFalse(output.contains("varsome.com/variant"));
        assertFalse(output.contains("gnomad.broadinstitute.org/variant"));
    }

    @Test
    void testLegacyTemplateRendersVariantDataAndReportLinks() {
        Sample sample = Sample.builder()
                .hpoIds(List.of("HP:000001", "HP:000002"))
                .build();
        Analysis analysis = Analysis.builder()
                .addStep(new FrequencyFilter(0.1f))
                .addStep(new PathogenicityFilter(true))
                .build();
        List<Gene> genes = List.of(fgfr2Gene, shhGene);
        AnalysisResults analysisResults = buildAnalysisResults(sample, analysis, genes, List.of());

        String output = HtmlResultsWriter.legacy().writeString(analysisResults, OutputSettings.builder().build());

        assertTrue(output.contains("FGFR2"));
        // compact/gnomAD-style display identifier rather than the old hgvsGenome form
        assertTrue(output.contains("10-123256215-T-G"));
        // configurable external-information links (report-links.yml)
        assertTrue(output.contains("varsome.com"));
    }

    @Test
    void testWriteTemplateWithPhenotypeOnlyAnalysis() throws Exception {
        Gene rbma8Gene = TestFactory.newGeneRBM8A();
        rbma8Gene.addGeneScore(GeneScore.builder().geneIdentifier(rbma8Gene.geneIdentifier()).modeOfInheritance(ModeOfInheritance.AUTOSOMAL_DOMINANT).combinedScore(0.85).build());

        List<Gene> genes = List.of(fgfr2Gene, rbma8Gene, shhGene);

        Sample sample = Sample.builder()
                .hpoIds(List.of("HP:000001", "HP:000002"))
                .build();

        Analysis analysis = Analysis.builder()
                .addStep(new PhivePriority(TestPriorityServiceFactory.stubPriorityService()))
                .build();
        AnalysisResults analysisResults = buildAnalysisResults(sample, analysis, genes, List.of());

        OutputSettings settings = OutputSettings.builder().build();

        String output = instance.writeString(analysisResults, settings);
        assertTrue(output.contains("Gene Symbol"));
        assertTrue(output.contains("FGFR2"));
        assertTrue(output.contains("RBM8A"));
    }
}
