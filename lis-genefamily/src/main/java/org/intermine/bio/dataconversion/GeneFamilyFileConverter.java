package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
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
 * Load gene family data from LIS datastore files.
 *
 * @author Sam Hokin
 */
public class GeneFamilyFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GeneFamilyFileConverter.class);

    // Items to store
    List<Item> ontologyAnnotations = new ArrayList<>();
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
        if (getCurrentFile().getName().endsWith(".info_annot_ahrd.tsv")) {
	    // legume.genefam.fam1.M65K.info_annot_ahrd.tsv
            String fastaDirname = getCurrentFile().getParent()+"/"+getCurrentFile().getName().replace("info_annot_ahrd.tsv", "family_fasta");
            printInfoBlurb(fastaDirname);
            processInfoAnnotAhrdFile(reader);
	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws Exception {
	store(dataSource);
	store(dataSets.values());
	store(organisms.values());
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
        Item geneFamily;
        if (geneFamilies.containsKey(identifier)) {
            geneFamily = geneFamilies.get(identifier);
        } else {
            geneFamily = createItem("GeneFamily");
            geneFamily.setAttribute("identifier", identifier);
            geneFamilies.put(identifier, geneFamily);
        }
        return geneFamily;
    }

    /**
     * Get/add an OntologyTerm Item, keyed by identifier
     */
    public Item getOntologyTerm(String identifier) {
        Item ontologyTerm;
        if (ontologyTerms.containsKey(identifier)) {
            ontologyTerm = ontologyTerms.get(identifier);
        } else {
            ontologyTerm = createItem("OntologyTerm");
            ontologyTerm.setAttribute("identifier", identifier);
            ontologyTerms.put(identifier, ontologyTerm);
        }
        return ontologyTerm;
    }

    /**
     * Get/add a Protein Item, keyed by primaryIdentifier
     */
    public Item getProtein(String primaryIdentifier) {
        Item protein;
        if (proteins.containsKey(primaryIdentifier)) {
            protein = proteins.get(primaryIdentifier);
        } else {
            protein = createItem("Protein");
            protein.setAttribute("primaryIdentifier", primaryIdentifier);
	    String secondaryIdentifier = extractSecondaryIdentifier(primaryIdentifier, true);
	    if (secondaryIdentifier!=null) protein.setAttribute("secondaryIdentifier", secondaryIdentifier);
            proteins.put(primaryIdentifier, protein);
        }
        return protein;
    }

    /**
     * Get/add a Gene Item, keyed by primaryIdentifier
     */
    public Item getGene(String primaryIdentifier) {
        Item gene;
        if (genes.containsKey(primaryIdentifier)) {
            gene = genes.get(primaryIdentifier);
        } else {
            gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", primaryIdentifier);
	    String secondaryIdentifier = extractSecondaryIdentifier(primaryIdentifier, true);
	    if (secondaryIdentifier!=null) gene.setAttribute("secondaryIdentifier", secondaryIdentifier);
            genes.put(primaryIdentifier, gene);
        }
        return gene;
    }

    /**
     * Get/add a ProteinDomain Item.
     */
    public Item getProteinDomain(String identifier) {
        Item proteinDomain;
        if (proteinDomains.containsKey(identifier)) {
            proteinDomain = proteinDomains.get(identifier);
        } else {
            proteinDomain = createItem("ProteinDomain");
            proteinDomain.setAttribute("primaryIdentifier", identifier);
	    proteinDomains.put(identifier, proteinDomain);
        }
        return proteinDomain;
    }

    /**
     * Process an info_annot_ahrd.tsv file which contains gene families and semi-colon separated groups of ontology terms.
     * 0      1       2    3    4               5
     * lis.genefam.fam1.M65K.info_annot_ahrd.tsv
     * legfed_v1_0.L_LFXSXJ-consensus  splicing factor 3B subunit 3-like isoform X2 [Glycine max]; 
     *                                 IPR004871 (Cleavage/polyadenylation specificity factor, A subunit, C-terminal); 
     *                                 GO:0003676 (nucleic acid binding), GO:0005634 (nucleus)
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
	Item dataSet = getDataSet();
	// get the FASTA directory
        String fastaDirname = getCurrentFile().getParent()+"/"+getCurrentFile().getName().replace("info_annot_ahrd.tsv", "family_fasta");
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
            geneFamily.setReference("dataSet", dataSet);
            // GO terms
            for (String identifier : record.go.keySet()) {
                String description = record.go.get(identifier);
                Item goTerm = getOntologyTerm(identifier);
                goTerm.setAttribute("description", description);
                Item goAnnotation = createItem("OntologyAnnotation");
                goAnnotation.setReference("subject", geneFamily);
                goAnnotation.setReference("ontologyTerm", goTerm);
                goAnnotation.addToCollection("dataSets", dataSet);
		ontologyAnnotations.add(goAnnotation);
            }
            // interpro protein domains
            for (String identifier : record.interpro.keySet()) {
                Item proteinDomain = getProteinDomain(identifier);
                String description = record.interpro.get(identifier);
                proteinDomain.setAttribute("description", description);
                proteinDomain.addToCollection("geneFamilies", geneFamily);
            }
            // load the gene family FASTA if present to link proteins
            String fastaFilename = fastaDirname+"/"+record.identifier;
            File fastaFile = new File(fastaFilename);
            if (fastaFile.exists()) {
                BufferedReader fbr = new BufferedReader(new FileReader(fastaFile));
                String fline = null;
                while ((fline=fbr.readLine())!=null) {
                    if (fline.startsWith(">")) {
                        String name = fline.substring(1);
                        String[] parts = name.split("\\.");
                        String gensp = parts[0];
                        Item organism = getOrganism(gensp);
                        Item protein = getProtein(name);
			String geneIdentifier = extractGeneIdentifierFromProteinIdentifier(name);
			Item gene = getGene(geneIdentifier);
                        protein.setReference("geneFamily", geneFamily);
			protein.addToCollection("genes", gene);
                        protein.addToCollection("dataSets", dataSet);
			gene.setReference("geneFamily", geneFamily);
			gene.addToCollection("proteins", protein);
                        gene.addToCollection("dataSets", dataSet);
                    }
                }
                fbr.close();
            }
        }
        br.close();
    }
}
