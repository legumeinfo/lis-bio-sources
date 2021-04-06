package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2015-2016 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.Reader;
import java.io.IOException;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Loads data from LIS datastore files about organisms and strains.
 *
 * @author Sam Hokin
 */
public class AboutFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(AboutFileConverter.class);

    /**
     * Create a new DatastoreFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public AboutFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().contains("description_")) {
            processDescriptionFile(reader);
        } else if (getCurrentFile().getName().contains("strains_")) {
            processStrainsFile(reader);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
	store(organisms.values());
        store(strains.values());
	
    }

    /**
     * Process an organism description file, which is in YAML format:
     * 0           1         2
     * description_Phaseolus_vulgaris.yml
     *
     * %YAML 1.2
     * ##### Phaseolus vulgaris	
     * organism.taxid:	3885
     * organism.genus:	Phaseolus
     * organism.species:	vulgaris
     * organism.abbrev:	phavu
     * organism.commonName:	common bean
     * organism.description:	Common bean was likely domesticated independently both in Central America and in the Andes....
     */
    void processDescriptionFile(Reader reader) throws IOException {
        // get the organism
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String[] threeparts = dotparts[0].split("_");
        String genus = threeparts[1];
        String species = threeparts[2];
        Item organism = getOrganism(genus, species);
        // now load the attributes
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("%")) continue;
            if (line.startsWith("#")) continue;
            String[] parts = line.split("\t");
            if (parts.length>1) {
                String attributeName = parts[0].replace("organism.","").replace(":","");
                String attributeValue = parts[1].trim();
		if (attributeValue.length()>0) {
		    if (attributeName.equals("taxid")) {
			// do nothing, we get taxonId in getOrganism()
		    } else if (attributeName.equals("abbrev")) {
			organism.setAttribute("abbreviation", attributeValue);
		    } else {
			organism.setAttribute(attributeName, attributeValue);
		    }
                }
            }
        }
        br.close();
    }

    /**
     * Process a strains file, which is in YAML format:
     * 0       1         2
     * strains_Phaseolus_vulgaris.yml
     *
     * %YAML 1.2
     * #####
     * strain.identifier:	G19833
     * strain.accession:	
     * strain.name:	G19833
     * strain.origin:	Peru
     * strain.description:	Andean landrace G19833 was selected for genome sequencing partly due to its resistance to numerous diseases...
     * #####
     * strain.identifier:	BAT93
     * strain.accession:	PI 633451
     * strain.name:	BAT93
     * strain.origin:	CIAT
     * strain.description:	Accession BAT93 is a Mesomarican line that has been used in numerous breeding projects and trait-mapping studies.
     */
    void processStrainsFile(Reader reader) throws IOException {
        // get the organism
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String[] threeparts = dotparts[0].split("_");
	String genus = threeparts[1];
	String species = threeparts[2];
        Item organism = getOrganism(genus, species);
	Item strain = null; // the current strain in the strains file
        // spin through the strain sections
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("%") || line.startsWith("#")) continue;
	    if (line.startsWith("strain.identifier")) {
		// top of new strain section
		// strain.identifier:	W05
                String[] parts = line.split("\t");
		String strainId = parts[1].trim();
		strain = getStrain(strainId, organism);
            } else {
		// other strain attributes
		// strain.origin:	Shanxi Sheng, China
                String[] parts = line.split("\t");
                if (parts.length>1) {
                    String attributeName = parts[0].replace("strain.","").replace(":","");
                    String attributeValue = parts[1].trim();
		    strain.setAttribute(attributeName, attributeValue);
                }
            }
        }
        br.close();
    }
}
