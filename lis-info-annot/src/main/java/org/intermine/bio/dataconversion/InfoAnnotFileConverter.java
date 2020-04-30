package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.io.IOException;
import java.io.InputStream;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.Properties;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Loads data from an LIS datastore info_annot.txt file, e.g. phavu.G19833.gnm2.ann1.PB8d.info_annot.txt.
 *
 * @author Sam Hokin
 */
public class InfoAnnotFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(InfoAnnotFileConverter.class);

    // local things to store
    Map<String,Item> genes = new HashMap<>();
    Map<String,Item> proteins = new HashMap<>();
    Map<String,Item> mRNAs = new HashMap<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> ontologyAnnotations = new HashMap<>(); // keyed by identifier_version_subject

    // we store some ontologies here
    Item geneOntology;
    Item pfamOntology;
    Item pantherOntology;
    Item kogOntology;
    Item ecOntology;
    Item koOntology;

    /**
     * Create a new InfoAnnotFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public InfoAnnotFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
	dataSource = getDataSource();
	// GO
        geneOntology = createItem("Ontology");
        geneOntology.setAttribute("name", "GO");
        geneOntology.setAttribute("url", "http://www.geneontology.org");
        // Pfam
        pfamOntology = createItem("Ontology");
        pfamOntology.setAttribute("name", "Pfam");
        pfamOntology.setAttribute("url", "https://pfam.xfam.org/");
        // PANTHER
        pantherOntology = createItem("Ontology");
        pantherOntology.setAttribute("name", "PANTHER");
        pantherOntology.setAttribute("url", "http://www.pantherdb.org/");
        // KOG
        kogOntology = createItem("Ontology");
        kogOntology.setAttribute("name", "KOG");
        kogOntology.setAttribute("url", "https://genome.jgi.doe.gov/Tutorial/tutorial/kog.html");
        // ENZYME
        ecOntology = createItem("Ontology");
        ecOntology.setAttribute("name", "ENZYME");
        ecOntology.setAttribute("url", "https://enzyme.expasy.org/");
        // KEGG Ontology
        koOntology = createItem("Ontology");
        koOntology.setAttribute("name", "KEGG Orthology");
        koOntology.setAttribute("url", "https://www.genome.jp/kegg/ko.html");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().endsWith(".info_annot.txt")) {
            processInfoAnnotFile(reader);
	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
	store(dataSource);
	store(dataSets.values());
	store(organisms.values());
	store(strains.values());
        store(geneOntology);
        store(pfamOntology);
        store(pantherOntology);
        store(kogOntology);
        store(ecOntology);
        store(koOntology);
	store(ontologyTerms.values());
        store(ontologyAnnotations.values());
        store(genes.values());
        store(proteins.values());
        store(mRNAs.values());
    }

    /**
     * Get/add a Gene Item, keyed by primaryIdentifier
     */
    public Item getGene(String primaryIdentifier) {
        Item gene;
        if (genes.containsKey(primaryIdentifier)) {
            gene = genes.get(primaryIdentifier);
        } else {
            // phavu.Phvul.002G040500
            gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", primaryIdentifier);
	    String secondaryIdentifier = extractSecondaryIdentifier(primaryIdentifier, true);
	    if (secondaryIdentifier!=null) gene.setAttribute("secondaryIdentifier", secondaryIdentifier);
            genes.put(primaryIdentifier, gene);
        }
        return gene;
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
     * Get/add an MRNA Item, keyed by primaryIdentifier (!)
     */
    public Item getMRNA(String primaryIdentifier) {
        Item mRNA;
        if (mRNAs.containsKey(primaryIdentifier)) {
            mRNA = mRNAs.get(primaryIdentifier);
        } else {
            // Phvul.002G040500.1
            mRNA = createItem("MRNA");
            mRNA.setAttribute("primaryIdentifier", primaryIdentifier);
	    String secondaryIdentifier = extractSecondaryIdentifier(primaryIdentifier, true);
	    if (secondaryIdentifier!=null) mRNA.setAttribute("secondaryIdentifier", secondaryIdentifier);
            mRNAs.put(primaryIdentifier, mRNA);
        }
        return mRNA;
    }

    /**
     * Form a key for an OntologyAnnotation for dupe avoidance.
     */
    public static String formAnnotKey(String termIdentifier, String subjectId) { 
        return termIdentifier+"_"+subjectId;
    }

    /**
     * Process an info_annot.txt file which contains relationships between genes, transcripts, proteins and ontology terms.
     * This will also link genes to proteins. The file name starts with gensp.strain.assembly.annotation.
     * 0     1      2    3    4    5          6
     * phavu.G19833.gnm2.ann1.PB8d.info_annot.txt
     *
     * Note: gensp.strain.assembly.annotation must be prepended to names in this file, which only contain the core IDs.
     *
     * pacId locusName transcriptName peptideName Pfam Panther KOG ec KO GO Best-hit-arabi-name arabi-symbol arabi-defline
     * 37170591 Phvul.001G000400 Phvul.001G000400.1 Phvul.001G000400.1.p PF00504 PTHR21649,PTHR21649:SF24 1.10.3.9 K14172 GO:0016020,GO:0009765 AT1G76570.1 Chlorophyll family protein
     */
    void processInfoAnnotFile(Reader reader) throws IOException {
	Item dataSet = getDataSet();
        String assemblyVersion = extractAssemblyVersion(getCurrentFile().getName());
        String annotationVersion = extractAnnotationVersion(getCurrentFile().getName());
	String gensp = extractGensp(getCurrentFile().getName());
	String strainId = extractStrainIdentifier(getCurrentFile().getName());
	// organism and strain
        Item organism = getOrganism(gensp);
        Item strain = getStrain(strainId, organism);
        // spin through the file
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            // comment line
            if (line.startsWith("#")) continue;
            // String bestHitAtName;
            // String bestHitAtSymbol;
            // String bestHitAtDefline;
            InfoAnnotRecord record = new InfoAnnotRecord(line);
            if (record.pacId!=null) {
                // the gene
                String geneIdentifier = formPrimaryIdentifier(gensp, strainId, assemblyVersion, annotationVersion, record.locusName);
                Item gene = getGene(geneIdentifier);
                gene.setAttribute("assemblyVersion", assemblyVersion);
                gene.setAttribute("annotationVersion", annotationVersion);
                gene.setReference("organism", organism);
                gene.setReference("strain", strain);
                gene.addToCollection("dataSets", dataSet);
                // the protein
                String proteinIdentifier = formPrimaryIdentifier(gensp, strainId, assemblyVersion, annotationVersion, record.peptideName);
                Item protein = getProtein(proteinIdentifier);
                protein.setAttribute("assemblyVersion", assemblyVersion);
                protein.setAttribute("annotationVersion", annotationVersion);
                protein.setReference("organism", organism);
                protein.setReference("strain", strain);
                protein.addToCollection("genes", gene);
                protein.addToCollection("dataSets", dataSet);
                // the transcript = mRNA
                String mRNAIdentifier = formPrimaryIdentifier(gensp, strainId, assemblyVersion, annotationVersion, record.transcriptName);
                Item mRNA = getMRNA(mRNAIdentifier);
                mRNA.setAttribute("assemblyVersion", assemblyVersion);
                mRNA.setAttribute("annotationVersion", annotationVersion);
                mRNA.setReference("gene", gene); 
                mRNA.setReference("organism", organism);
                mRNA.setReference("strain", strain);
                mRNA.addToCollection("dataSets", dataSet);
                mRNA.setReference("protein", protein);
                // GO terms
                for (String identifier : record.GO) {
                    Item term = getOntologyTerm(identifier);
                    term.setReference("ontology", geneOntology);
                    String annotKey = formAnnotKey(identifier, geneIdentifier);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item annotation = createItem("OntologyAnnotation");
                        annotation.setReference("subject", gene);
                        annotation.setReference("ontologyTerm", term);
                        annotation.addToCollection("dataSets", dataSet);
                        ontologyAnnotations.put(annotKey, annotation);
                    }
                }
                // Pfam terms
                for (String identifier : record.Pfam) {
                    Item term = getOntologyTerm(identifier);
                    term.setReference("ontology", pfamOntology);
                    String annotKey = formAnnotKey(identifier, proteinIdentifier);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item annotation = createItem("OntologyAnnotation");
                        annotation.setReference("subject", protein);
                        annotation.setReference("ontologyTerm", term);
                        annotation.addToCollection("dataSets", dataSet);
                        ontologyAnnotations.put(annotKey, annotation);
                    }
                }
                // Panther terms
                for (String identifier : record.Panther) {
                    Item term = getOntologyTerm(identifier);
                    term.setReference("ontology", pantherOntology);
                    String annotKey = formAnnotKey(identifier, proteinIdentifier);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item annotation = createItem("OntologyAnnotation");
                        annotation.setReference("subject", protein);
                        annotation.setReference("ontologyTerm", term);
                        annotation.addToCollection("dataSets", dataSet);
                        ontologyAnnotations.put(annotKey, annotation);
                    }                }
                // KOG terms
                for (String identifier : record.KOG) {
                    Item term = getOntologyTerm(identifier);
                    term.setReference("ontology", kogOntology);
                    String annotKey = formAnnotKey(identifier, proteinIdentifier);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item annotation = createItem("OntologyAnnotation");
                        annotation.setReference("subject", protein);
                        annotation.setReference("ontologyTerm", term);
                        annotation.addToCollection("dataSets", dataSet);
                        ontologyAnnotations.put(annotKey, annotation);
                    }
                }
                // ec terms
                for (String identifier : record.ec) {
                    Item term = getOntologyTerm(identifier);
                    term.setReference("ontology", ecOntology);
                    String annotKey = formAnnotKey(identifier, proteinIdentifier);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item annotation = createItem("OntologyAnnotation");
                        annotation.setReference("subject", protein);
                        annotation.setReference("ontologyTerm", term); 
                        annotation.addToCollection("dataSets", dataSet);
                        ontologyAnnotations.put(annotKey, annotation);
                    }
                }
                // KO terms
                for (String identifier : record.KO) {
                    Item term = getOntologyTerm(identifier);
                    term.setReference("ontology", koOntology);
                    String annotKey = formAnnotKey(identifier, geneIdentifier);
                    if (!ontologyAnnotations.containsKey(annotKey)) {
                        Item annotation = createItem("OntologyAnnotation");
                        annotation.setReference("subject", gene);
                        annotation.setReference("ontologyTerm", term);
                        annotation.addToCollection("dataSets", dataSet);
                        ontologyAnnotations.put(annotKey, annotation);
                    }
                }
            }
        }
        br.close();
    }
}
