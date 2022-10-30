package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;

import java.util.Arrays;
import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Loads Phytozome gene family identifiers and descriptions from a tab-delimited file.
 *
 * There is no README, so we use project.xml properties to define the DataSet.
 *
 * NOTE: not parsing out the IPR and GO terms into OntologyAnnotations yet!
 *
 * @author Sam Hokin
 */
public class PhytozomeFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(PhytozomeFileConverter.class);

    // local Items to store
    List<Item> ontologyAnnotations = new LinkedList<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> geneFamilies = new HashMap<>();

    /**
     * Create a new PhytozomeFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public PhytozomeFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        // DataSet
        if (dataSetName==null || dataSetDescription==null) {
            throw new RuntimeException("ERROR: dataSetName and dataSetDescription must be set in project.xml.");
        }
        dataSet = createItem("DataSet");
        dataSet.setAttribute("name", dataSetName);
        dataSet.setAttribute("description", dataSetDescription);
        if (dataSetUrl!=null) {
            dataSet.setAttribute("url", dataSetUrl);
        }
        if (dataSetLicence!=null) {
            dataSet.setAttribute("licence", dataSetLicence);
        } else {
            dataSet.setAttribute("licence", DatastoreFileConverter.DEFAULT_DATASET_LICENCE);
        }
        dataSet.setReference("dataSource", dataSource);
        // process files
        if (getCurrentFile().getName().startsWith("phytozome") &&
            getCurrentFile().getName().endsWith("descriptors.txt")) {
            processPhytozomeFile(reader);
	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws Exception {
        store(dataSource);
        store(dataSet);
	store(geneFamilies.values());
	// store(ontologyTerms.values());
	// store(ontologyAnnotations);
    }

    /**
     * Process a Phytozome descriptor file, which is simply ID and comment.
     * 0
     * phytozome_10_2.59086617
     * 1
     * L-lactate dehydrogenase A-like [Glycine max]
     * IPR001557 (L-lactate/malate dehydrogenase)
     * GO:0003824 (catalytic activity),
     * GO:0005975 (carbohydrate metabolic process),
     * GO:0044262 (cellular carbohydrate metabolic process),
     * GO:0055114 (oxidation-reduction process)
     * *-*- gi|356497151|ref|XP_003517426.1|
     *
     * 0
     * phytozome_10_2.59026833
     * 1
     * 30S ribosomal protein S13
     * IPR001892 (Ribosomal protein S13), 
     * IPR010979 (Ribosomal protein S13-like, H2TH),
     * IPR027437 (30s ribosomal protein S13, C-terminal)
     * GO:0003676 (nucleic acid binding),
     * GO:0003723 (RNA binding),
     * GO:0003735 (structural constituent of ribosome),
     * GO:0005622 (intracellular),
     * GO:0005840 (ribosome),
     * GO:0006412 (translation)
     * ***- Medtr3g073140.1
     */
    void processPhytozomeFile(Reader reader) throws IOException {
        // spin through the file
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\\t");
            String identifier = fields[0];
            String description = fields[1];
            Item geneFamily = getGeneFamily(identifier);
            geneFamily.setAttribute("description", description);
        }
    }

    /**
     * Get/add a GeneFamily Item.
     */
    public Item getGeneFamily(String primaryIdentifier) {
        if (geneFamilies.containsKey(primaryIdentifier)) {
            return geneFamilies.get(primaryIdentifier);
        } else {
            Item geneFamily = createItem("GeneFamily");
            geneFamily.setAttribute("primaryIdentifier", primaryIdentifier);
            geneFamilies.put(primaryIdentifier, geneFamily);
            return geneFamily;
        }
    }

    /**
     * Get/add an OntologyTerm Item, keyed by identifier
     */
    public Item getOntologyTerm(String identifier) {
        if (ontologyTerms.containsKey(identifier)) {
            return ontologyTerms.get(identifier);
        } else {
            Item ontologyTerm = createItem("OntologyTerm");
            ontologyTerm.setAttribute("identifier", identifier);
            ontologyTerms.put(identifier, ontologyTerm);
            return ontologyTerm;
        }
    }

}
