/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.charite.compbio.exomiser.cli.options;

import de.charite.compbio.exomiser.core.ExomiserSettings;
import static de.charite.compbio.exomiser.core.ExomiserSettings.DISEASE_ID_OPTION;
import org.apache.commons.cli.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class DiseaseIdOptionMarshaller extends AbstractOptionMarshaller {

    private static final Logger logger = LoggerFactory.getLogger(DiseaseIdOptionMarshaller.class);

    public DiseaseIdOptionMarshaller() {
        option = new Option("D", DISEASE_ID_OPTION, true, "OMIM ID for disease being sequenced. e.g. OMIM:101600");
    }

    @Override
    public void applyValuesToSettingsBuilder(String[] values, ExomiserSettings.SettingsBuilder settingsBuilder) {
        settingsBuilder.diseaseId(values[0]);
    }  
    
}
