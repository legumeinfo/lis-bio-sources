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

import org.ncgr.zip.GZIPBufferedReader;

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

    boolean ahrdFileFound = false;

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
        dataSet.setAttribute("synopsis", dataSetDescription);
        if (dataSetLicence!=null) {
            dataSet.setAttribute("licence", dataSetLicence);
        } else {
            dataSet.setAttribute("licence", DatastoreFileConverter.DEFAULT_DATASET_LICENCE);
        }
        dataSet.setReference("dataSource", dataSource);
        // process files
        if (getCurrentFile().getName().endsWith(".info_annot_ahrd.tsv.gz")) {
            processInfoAnnotAhrdFile();
            ahrdFileFound = true;
	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws Exception {
        if (!ahrdFileFound) {
            throw new RuntimeException("File ending in .info_annot_ahrd.tsv.gz not found. Aborting");
        }
        store(dataSource);
        store(dataSet);
	store(geneFamilies.values());
	store(ontologyTerms.values());
	store(ontologyAnnotations);
	store(proteins.values());
	store(genes.values());
	store(proteinDomains.values());
    }

    /**
     * Process an info_annot_ahrd.tsv.gz file which contains gene families and semi-colon separated groups of ontology terms.
     * 0   1       2    3    4               5   6
     * lis.genefam.fam1.M65K.info_annot_ahrd.tsv.gz
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
    void processInfoAnnotAhrdFile() throws IOException {
        // spin through the AHRD file lines
        BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment line
            // parse record and create items
            InfoAnnotAhrdRecord record = new InfoAnnotAhrdRecord(line);
            Item geneFamily = getGeneFamily(record.identifier);
            geneFamily.setAttribute("version", record.version);
            geneFamily.setAttribute("description", record.description);
            // GO terms - only store identifier
            for (String identifier : record.go.keySet()) {
                Item goTerm = getOntologyTerm(identifier);
                Item goAnnotation = createItem("OntologyAnnotation");
                goAnnotation.setReference("subject", geneFamily);
                goAnnotation.setReference("ontologyTerm", goTerm);
		ontologyAnnotations.add(goAnnotation);
            }
            // InterPro protein domains - only store identifier
            for (String identifier : record.interpro.keySet()) {
                Item proteinDomain = getProteinDomain(identifier);
                geneFamily.addToCollection("proteinDomains", proteinDomain);
            }
        }
        br.close();
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

    /**
     * Get/add a Protein Item, keyed by primaryIdentifier
     */
    public Item getProtein(String primaryIdentifier) {
        if (proteins.containsKey(primaryIdentifier)) {
            return proteins.get(primaryIdentifier);
        } else {
            Item protein = createItem("Protein");
            protein.setAttribute("primaryIdentifier", primaryIdentifier);
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
            genes.put(primaryIdentifier, gene);
            return gene;
        }
    }

    /**
     * Get/add a ProteinDomain Item.
     */
    public Item getProteinDomain(String primaryIdentifier) {
        if (proteinDomains.containsKey(primaryIdentifier)) {
            return proteinDomains.get(primaryIdentifier);
        } else {
            Item proteinDomain = createItem("ProteinDomain");
            proteinDomain.setAttribute("primaryIdentifier", primaryIdentifier);
	    proteinDomains.put(primaryIdentifier, proteinDomain);
            return proteinDomain;
        }
    }

}
