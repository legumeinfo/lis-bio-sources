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
import java.io.IOException;
import java.io.Reader;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Store Phenotype/ontology term data from a tab-delimited file.
 *
 * #phenotype	ontologyId
 * Seed protein	SOY:0001676
 * Seed oil	SOY:0001668
 *
 * @author Sam Hokin
 */
public class PhenotypeFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(PhenotypeFileConverter.class);

    // things to store
    List<Item> dataSets = new ArrayList<>();
    List<Item> phenotypes = new ArrayList<>();
    List<Item> ontologyAnnotations = new ArrayList<>();
    Map<String,Item> ontologyTermMap = new HashMap<>();

    /**
     * Create a new PhenotypeFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public PhenotypeFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }
    
    /**
     * {@inheritDoc}
     * Process the phenotype-ontology term annotations by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (dataSource==null) {
	    dataSource = getDataSource();
        }

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

	Item dataSet = getDataSet();
	dataSets.add(dataSet);
	
        BufferedReader bufferedReader = new BufferedReader(reader);
	String line;
        while ((line=bufferedReader.readLine())!=null) {
            String[] parts = line.split("\t");
            if (line.startsWith("#") || line.trim().length()==0 || parts.length<2) {
                continue;
            }
            String phenotypeId = parts[0];
            String ontologyTermId = parts[1];
	    Item phenotype = createItem("Phenotype");
	    phenotype.setAttribute("primaryIdentifier", phenotypeId);
	    phenotype.setReference("dataSource", dataSource);
	    phenotype.setReference("dataSet", dataSet);
	    phenotypes.add(phenotype);
	    Item ontologyTerm = ontologyTermMap.get(ontologyTermId);
	    if (ontologyTerm==null) {
		ontologyTerm = createItem("OntologyTerm");
		ontologyTerm.setAttribute("identifier", ontologyTermId);
		ontologyTermMap.put(ontologyTermId, ontologyTerm);
	    }
	    Item ontologyAnnotation = createItem("OntologyAnnotation");
	    ontologyAnnotation.setReference("ontologyTerm", ontologyTerm);
	    ontologyAnnotation.setReference("subject", phenotype);
	    ontologyAnnotations.add(ontologyAnnotation);
        }
        bufferedReader.close();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
	store(dataSource);
	store(dataSets);
	store(phenotypes);
	store(ontologyTermMap.values());
	store(ontologyAnnotations);
    }
}
