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

import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Loads data from LIS datastore files, using the file names to detect what sort of files they are, and,
 * sometimes, what the organism is.
 *
 * @author Sam Hokin
 */
public class AboutFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(AboutFileConverter.class);

    Item dataSet;
    Item organism;
    Map<String,Item> strainMap = new HashMap<>();

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
	// DataSet attributes are given entirely in project.xml.
	dataSource = getDataSource();
	dataSet = getDataSet();
        // process the file
        if (getCurrentFile().getName().contains("description_")) {
            // description_Phaseolus_vulgaris.yml
            processDescriptionFile(reader);
        } else if (getCurrentFile().getName().contains("strains_")) {
            // strains_Phaseolus_vulgaris.yml
            processStrainsFile(reader);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
	store(dataSource);
	store(dataSet);
	store(organism);
        store(strainMap.values());
    }

    /**
     * Process an organism description file, which is in YAML format:
     *
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
        String[] dashparts = dotparts[0].split("_");
        String genus = null;   // for forming Organism.name
        String species = null; // for forming Organism.name
	String taxonId = null;
        Item organism = getOrganism(dashparts);
	organism.addToCollection("dataSources", dataSource);
	organism.addToCollection("dataSets", dataSet);
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
		// munge
                if (attributeName.equals("taxid")) {
		    attributeName = "taxonId";
		    taxonId = attributeValue;
		}
                if (attributeName.equals("abbrev")) attributeName = "abbreviation";
                if (attributeValue.length()>0) {
                    organism.setAttribute(attributeName, attributeValue);
                }
            }
        }
        br.close();
    }

    /**
     * Process a strains file, which is in YAML format:
     *
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
        String[] dashparts = dotparts[0].split("_");
        organism = getOrganism(dashparts);
        // spin through the strain sections
        BufferedReader br = new BufferedReader(reader);
        Map<String,String> attributes = new HashMap<>();
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#####")) {
                // new strain section, store previous strain
                if (attributes.size()>0) {
                    String strainId = attributes.get("identifier");
                    Item strain = getStrain(strainId, organism);
		    strain.addToCollection("dataSources", dataSource);
		    strain.addToCollection("dataSets", dataSet);
                    for (String name : attributes.keySet()) {
                        String value = attributes.get(name);
                        strain.setAttribute(name, value);
                    }
		    strainMap.put(strainId, strain);
                }
                // clear attributes map
                attributes = new HashMap<>();
            } else if (line.startsWith("#") || line.startsWith("%")) {
                // comment
                continue;
            } else {
                // put strain attribute into map
                String[] parts = line.split("\t");
                if (parts.length>1) {
                    String attributeName = parts[0].replace("strain.","").replace(":","");
                    String attributeValue = parts[1].trim();
                    attributes.put(attributeName, attributeValue);
                }
            }
        }
        br.close();
        // last one
        if (attributes.size()>0) {
            String strainId = attributes.get("identifier");
            Item strain = getStrain(strainId, organism);
	    strain.addToCollection("dataSources", dataSource);
	    strain.addToCollection("dataSets", dataSet);
            for (String name : attributes.keySet()) {
                String value = attributes.get(name);
                strain.setAttribute(name, value);
            }
	    strainMap.put(strainId, strain);
        }
    }
}
