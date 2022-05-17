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
 * Load gene family data from LIS datastore files.
 *
 * There is no README, so we use project.xml properties to define the DataSet.
 *
 * @author Sam Hokin
 */
public class GeneFamilyFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GeneFamilyFileConverter.class);

    // local Items to store
    List<Item> ontologyAnnotations = new LinkedList<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> proteins = new HashMap<>();
    Map<String,Item> genes = new HashMap<>();
    Map<String,Item> proteinDomains = new HashMap<>();
    Map<String,Item> geneFamilies = new HashMap<>();

    /**
     * Create a new GeneFamilyFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GeneFamilyFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        // there is no README (yet)
        // DataSet
        if (dataSetName==null || dataSetUrl==null || dataSetDescription==null) {
            throw new RuntimeException("ERROR: dataSetName, dataSetUrl, and dataSetDescription must be set in project.xml.");
        }
        dataSet = createItem("DataSet");
        dataSet.setAttribute("name", dataSetName);
        dataSet.setAttribute("url", dataSetUrl);
        dataSet.setAttribute("description", dataSetDescription);
        if (dataSetLicence!=null) {
            dataSet.setAttribute("licence", dataSetLicence);
        } else {
            dataSet.setAttribute("licence", DatastoreFileConverter.DEFAULT_DATASET_LICENCE);
        }
        dataSet.setReference("dataSource", dataSource);
        // process files
        if (getCurrentFile().getName().endsWith(".info_annot_ahrd.tsv")) {
            processInfoAnnotAhrdFile(reader);
	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws Exception {
	store(geneFamilies.values());
	store(ontologyTerms.values());
	store(ontologyAnnotations);
	store(proteins.values());
	store(genes.values());
	store(proteinDomains.values());
    }

    /**
     * Get/add a GeneFamily Item.
     */
    public Item getGeneFamily(String identifier) {
        if (geneFamilies.containsKey(identifier)) {
            return geneFamilies.get(identifier);
        } else {
            Item geneFamily = createItem("GeneFamily");
            geneFamily.setAttribute("identifier", identifier);
            geneFamilies.put(identifier, geneFamily);
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

    /**
     * Get/add a Protein Item, keyed by primaryIdentifier
     */
    public Item getProtein(String primaryIdentifier) {
        if (proteins.containsKey(primaryIdentifier)) {
            return proteins.get(primaryIdentifier);
        } else {
            Item protein = createItem("Protein");
            protein.setAttribute("primaryIdentifier", primaryIdentifier);
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, true);
	    if (secondaryIdentifier!=null) protein.setAttribute("secondaryIdentifier", secondaryIdentifier);
            proteins.put(primaryIdentifier, protein);
            return protein;
        }
    }

    /**
     * Get/add a Gene Item, keyed by primaryIdentifier
     */
    public Item getGene(String primaryIdentifier) {
        if (genes.containsKey(primaryIdentifier)) {
            return genes.get(primaryIdentifier);
        } else {
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", primaryIdentifier);
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, true);
	    if (secondaryIdentifier!=null) gene.setAttribute("secondaryIdentifier", secondaryIdentifier);
            genes.put(primaryIdentifier, gene);
            return gene;
        }
    }

    /**
     * Get/add a ProteinDomain Item.
     */
    public Item getProteinDomain(String identifier) {
        if (proteinDomains.containsKey(identifier)) {
            return proteinDomains.get(identifier);
        } else {
            Item proteinDomain = createItem("ProteinDomain");
            proteinDomain.setAttribute("primaryIdentifier", identifier);
	    proteinDomains.put(identifier, proteinDomain);
            return proteinDomain;
        }
    }

    /**
     * Process an info_annot_ahrd.tsv file which contains gene families and semi-colon separated groups of ontology terms.
     * 0      1       2    3    4               5
     * lis.genefam.fam1.M65K.info_annot_ahrd.tsv
     * legfed_v1_0.L_LFXSXJ splicing factor 3B subunit 3-like isoform X2 [Glycine max];
     *                      IPR004871 (Cleavage/polyadenylation specificity factor, A subunit, C-terminal); 
     *                      GO:0003676 (nucleic acid binding), GO:0005634 (nucleus)
     *
     * legume.genefam.fam1.M65K.family_fasta/legfed_v1_0.L_LFXSXJ
     *                                    ^^^^^^^^
     * >lupan.Lup015831.1
     * ----CASFAKLT--TLSPHWIGNNSFSSRRGGSSPLTATRRVSLPIRASSYSDELVQTAK
     * TIASPGRGILAIDESNATCGKRLASIGLDNTEVNRQAYRQLLLTTPGLGEYISGAILFEE
     * ...
     * >phavu.Phvul.007G033800.1
     * -----------------------------------TFSPRRVSLPIRASSYQQELVQTAK
     * SIASPGRGILAIDESNATCGKRLASIGLDNTEVNRQAYRQLLLTTPGLGEYISGAILFEE
     * ...
     */
    void processInfoAnnotAhrdFile(Reader reader) throws IOException {
        // spin through the AHRD file lines
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment line
            // parse record and create items
            InfoAnnotAhrdRecord record = new InfoAnnotAhrdRecord(line);
            Item geneFamily = getGeneFamily(record.identifier);
            geneFamily.setAttribute("version", record.version);
            geneFamily.setAttribute("description", record.description);
            // name is needed for merge with Phylotree
            // legfed_v1_0.L_LFXSXJ name=L_LFXSXJ
            String[] parts = record.identifier.split("\\.");
            if (parts.length==2) {
                geneFamily.setAttribute("name", parts[1]);
            } else {
                throw new RuntimeException("Gene family identifier "+record.identifier+" is not in dot-separated format.");
            }
            // GO terms
            for (String identifier : record.go.keySet()) {
                String description = record.go.get(identifier);
                Item goTerm = getOntologyTerm(identifier);
                goTerm.setAttribute("description", description);
                Item goAnnotation = createItem("OntologyAnnotation");
                goAnnotation.setReference("subject", geneFamily);
                goAnnotation.setReference("ontologyTerm", goTerm);
		ontologyAnnotations.add(goAnnotation);
            }
            // interpro protein domains
            for (String identifier : record.interpro.keySet()) {
                Item proteinDomain = getProteinDomain(identifier);
                String description = record.interpro.get(identifier);
                proteinDomain.setAttribute("description", description);
                proteinDomain.addToCollection("geneFamilies", geneFamily);
            }
        }
        br.close();
    }
}
