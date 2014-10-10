/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.charite.compbio.exomiser.core.model;

import de.charite.compbio.exomiser.core.filter.FilterResult;
import de.charite.compbio.exomiser.core.filter.FilterResultStatus;
import static de.charite.compbio.exomiser.core.filter.FilterResultStatus.FAIL;
import static de.charite.compbio.exomiser.core.filter.FilterResultStatus.PASS;
import de.charite.compbio.exomiser.core.filter.FilterType;
import de.charite.compbio.exomiser.core.filter.FrequencyFilterResult;
import de.charite.compbio.exomiser.core.filter.QualityFilterResult;
import de.charite.compbio.exomiser.core.frequency.FrequencyData;
import de.charite.compbio.exomiser.core.pathogenicity.PathogenicityData;
import jannovar.common.Genotype;
import jannovar.exome.Variant;
import jannovar.genotype.GenotypeCall;
import java.util.EnumSet;
import java.util.Set;
import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.notNullValue;
import static org.hamcrest.CoreMatchers.nullValue;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.Test;

/**
 * Tests for non-bean (i.e. logic-containing) methods in
 * {@code VariantEvaluation} class
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class VariantEvaluationTest {

    VariantEvaluation instance;
    private static final Integer QUALITY = 2;
    private static final Integer READ_DEPTH = 6;
    private static final Genotype HETEROZYGOUS = Genotype.HETEROZYGOUS;

    private static final FilterResult FAIL_FREQUENCY_RESULT = new FrequencyFilterResult(0.1f, FAIL);
    private static final FilterResult PASS_FREQUENCY_RESULT = new FrequencyFilterResult(0.1f, PASS);

    private static final FilterResult PASS_QUALITY_RESULT = new QualityFilterResult(0.45f, PASS);

    public VariantEvaluationTest() {

    }

    @Before
    public void setUp() {
        GenotypeCall genotypeCall = new GenotypeCall(HETEROZYGOUS, QUALITY, READ_DEPTH);
        byte chr = 1;
        Variant variant = new Variant(chr, 1, "A", "T", genotypeCall, 2.2f, "");
        instance = new VariantEvaluation(variant);
    }

    @Test
    public void testGetRef() {
        assertThat(instance.getRef(), equalTo("A"));
    }

    @Test
    public void testGetAlt() {
        assertThat(instance.getAlt(), equalTo("T"));

    }

    @Test
    public void testGetVariantReadDepth() {
        assertThat(instance.getVariantReadDepth(), equalTo(READ_DEPTH));
    }

    @Test
    public void testGetGeneSymbol() {
        assertThat(instance.getGeneSymbol(), equalTo("."));
    }

    @Test
    public void testGetVariantStartPosition() {
        assertThat(instance.getVariantStartPosition(), equalTo(1));
    }

    @Test
    public void testGetVariantEndPosition() {
        assertThat(instance.getVariantEndPosition(), equalTo(1));
    }

    @Test
    public void testThatTheConstructorDoesNotSetAFrequencyDataObject() {
        FrequencyData frequencyData = instance.getFrequencyData();
        assertThat(frequencyData, nullValue());
    }

    @Test
    public void testThatTheConstructorCreatesAnEmptyPathogenicityDataObject() {
        PathogenicityData pathogenicityData = instance.getPathogenicityData();
        assertThat(pathogenicityData, notNullValue());
        assertThat(pathogenicityData.getMutationTasterScore(), nullValue());
        assertThat(pathogenicityData.getPolyPhenScore(), nullValue());
        assertThat(pathogenicityData.getSiftScore(), nullValue());
        assertThat(pathogenicityData.getCaddScore(), nullValue());
        assertThat(pathogenicityData.hasPredictedScore(), is(false));
    }

    @Test
    public void testThatAddingAFilterResultUpdatesVariantScore() {

        assertThat(instance.getVariantScore(), equalTo(1.0f));

        //adding a FilterResult also updates the score of the VariantEvaluation 
        instance.addFilterResult(PASS_QUALITY_RESULT);

        assertThat(instance.getFilterResults().size(), equalTo(1));
        assertThat(instance.getVariantScore(), equalTo(PASS_QUALITY_RESULT.getScore()));
    }

    @Test
    public void testThatAddingTwoPassFilterResultsUpdatesVariantScore() {

        assertThat(instance.getVariantScore(), equalTo(1.0f));

        float expectedScore = instance.getVariantScore();

        expectedScore *= PASS_QUALITY_RESULT.getScore();
        expectedScore *= PASS_FREQUENCY_RESULT.getScore();
        //adding a FilterResult also updates the score of the VariantEvaluation 
        instance.addFilterResult(PASS_QUALITY_RESULT);
        instance.addFilterResult(PASS_FREQUENCY_RESULT);

        assertThat(instance.getFilterResults().size(), equalTo(2));
        assertThat(instance.getVariantScore(), equalTo(expectedScore));
    }
    
    @Test
    public void testThatAddingPassAndFailFilterResultsUpdatesVariantScore() {

        assertThat(instance.getVariantScore(), equalTo(1.0f));

        float expectedScore = instance.getVariantScore();

        expectedScore *= PASS_QUALITY_RESULT.getScore();
        expectedScore *= FAIL_FREQUENCY_RESULT.getScore();
        //adding a FilterResult also updates the score of the VariantEvaluation 
        instance.addFilterResult(PASS_QUALITY_RESULT);
        instance.addFilterResult(FAIL_FREQUENCY_RESULT);

//        assertThat(instance.getFilterResults().size(), equalTo(2));
        assertThat(instance.getVariantScore(), equalTo(expectedScore));
    }

    @Test
    public void testGetFailedFilterTypes() {
        Set<FilterType> expectedFilters = EnumSet.of(FAIL_FREQUENCY_RESULT.getFilterType());

        instance.addFilterResult(FAIL_FREQUENCY_RESULT);
        assertThat(instance.getFailedFilterTypes(), equalTo(expectedFilters));
    }
    
    @Test
    public void testGetFailedFilterTypesDontContainPassedFilterTypes() {
        Set<FilterType> expectedFilters = EnumSet.of(FAIL_FREQUENCY_RESULT.getFilterType());

        instance.addFilterResult(FAIL_FREQUENCY_RESULT);
        instance.addFilterResult(PASS_QUALITY_RESULT);

        assertThat(instance.getFailedFilterTypes(), equalTo(expectedFilters));
    }
    
    @Test
    public void testFailsFiltersWhenFailedFilterResultAdded() {
        instance.addFilterResult(FAIL_FREQUENCY_RESULT);
        
        assertThat(instance.passedFilters(), is(false));
    }

    @Test
    public void testVariantScoreIsUpdatedWhenFailedFilterResultAdded() {

        instance.addFilterResult(FAIL_FREQUENCY_RESULT);

        assertThat(instance.getVariantScore(), equalTo(FAIL_FREQUENCY_RESULT.getScore()));
    }

    @Test
    public void testVariantScoreIsUpdatedWhenPassedAndFailedFilterResultAdded() {
        float qualScore = PASS_QUALITY_RESULT.getScore();
        //adding a FilterResult also updates the score of the VariantEvaluation 
        instance.addFilterResult(PASS_QUALITY_RESULT);
        assertThat(instance.getVariantScore(), equalTo(qualScore));

        float freqScore = 0.1f;
        FilterResult frequencyScore = new FrequencyFilterResult(freqScore, FAIL);
        //adding a failed FilterResult also updates the score of the VariantEvaluation         
        instance.addFilterResult(frequencyScore);

        assertThat(instance.getVariantScore(), equalTo(qualScore * freqScore));
    }

    @Test
    public void testPassesFiltersWhenNoFiltersHaveBeenApplied() {
        assertThat(instance.getFailedFilterTypes().isEmpty(), is(true));
        assertThat(instance.getFilterResults().isEmpty(), is(true));
        assertThat(instance.passedFilters(), is(true));
    }

    @Test
    public void testPassesFilterIsTrueWhenPassedFilterResultAdded() {
        FilterType passedFilterType = PASS_QUALITY_RESULT.getFilterType();

        instance.addFilterResult(PASS_QUALITY_RESULT);
        instance.addFilterResult(FAIL_FREQUENCY_RESULT);

        assertThat(instance.passedFilter(passedFilterType), is(true));
    }

    @Test
    public void testPassesFilterIsFalseWhenFailedFilterResultAdded() {
        FilterType filterType = FAIL_FREQUENCY_RESULT.getFilterType();

        instance.addFilterResult(FAIL_FREQUENCY_RESULT);

        assertThat(instance.passedFilter(filterType), is(false));
    }

}