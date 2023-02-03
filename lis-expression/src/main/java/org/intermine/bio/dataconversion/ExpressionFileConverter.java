package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.ncgr.datastore.Readme;
import org.ncgr.datastore.validation.ExpressionCollectionValidator;
import org.ncgr.zip.GZIPBufferedReader;

/**
 * DataConverter to create ExpressionSource, ExpressionSample and ExpressionValue items from a single set of datastore expression files.
 *
 * expression/G19833.gnm1.ann1.expr.4ZDQ/
 * ├── phavu.G19833.gnm1.ann1.expr.4ZDQ.obo.tsv.gz
 * ├── phavu.G19833.gnm1.ann1.expr.4ZDQ.samples.tsv.gz
 * ├── phavu.G19833.gnm1.ann1.expr.4ZDQ.values.tsv.gz
 * └── README.G19833.gnm1.ann1.expr.4ZDQ.yml
 *
 * @author Sam Hokin
 */
public class ExpressionFileConverter extends DatastoreFileConverter {
    
    private static final Logger LOG = Logger.getLogger(ExpressionFileConverter.class);

    // README items
    Item expressionSource;

    // Lists
    List<Item> expressionValues = new ArrayList<>();
    List<Item> ontologyAnnotations = new ArrayList<>();

    // Maps
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> samples = new HashMap<>();
    Map<String,Item> genes = new HashMap<>();

    // validate the collection first by storing a flag
    boolean collectionValidated = false;

    /**
     * Constructor.
     *
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     * @throws ObjectStoreException os
     */
    public ExpressionFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * Called for each file found.
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws IOException {
        if (!collectionValidated) {
            ExpressionCollectionValidator validator = new ExpressionCollectionValidator(getCurrentFile().getParent());
            validator.validate();
            if (!validator.isValid()) {
                throw new RuntimeException("Collection "+getCurrentFile().getParent()+" does not pass validation.");
            }
            collectionValidated = true;
        }
	if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
            setStrain();
            // check required stuff for expression
            if (readme.expression_unit==null) {
                throw new RuntimeException("ERROR: a required field is missing from expression README "+getCurrentFile().getName()+": "+
                                           "Required fields are: expression_unit");
            }
            // ExpressionSource
            expressionSource = createItem("ExpressionSource");
            expressionSource.setAttribute("primaryIdentifier", readme.identifier);
            expressionSource.setAttribute("synopsis", readme.synopsis);
            expressionSource.setAttribute("description", readme.description);
            expressionSource.setAttribute("unit", readme.expression_unit);
            expressionSource.setReference("organism", organism);
            expressionSource.setReference("strain", strain);
            if (readme.geoseries!=null) expressionSource.setAttribute("geoSeries", readme.geoseries);
            if (readme.sraproject!=null) expressionSource.setAttribute("sra", readme.sraproject);
            if (readme.bioproject!=null) expressionSource.setAttribute("bioProject", readme.bioproject);
            if (publication!=null) expressionSource.addToCollection("publications", publication);
        } else if (getCurrentFile().getName().endsWith("samples.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
	    processSamplesFile();
        } else if (getCurrentFile().getName().endsWith("values.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processValuesFile();
        } else if (getCurrentFile().getName().endsWith("obo.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processOboFile();
        } else {
            System.out.println(" x skipping "+getCurrentFile().getName());
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        // DatastoreFileConverter items
        storeCollectionItems();
        // add references to samples
        for (Item sample : samples.values()) {
            sample.setReference("source", expressionSource);
        }
        // add references to genes
        for (Item gene : genes.values()) {
            gene.setReference("organism", organism);
            gene.setReference("strain", strain);
        }
        // local items
        store(expressionSource);
	store(genes.values());
        store(samples.values());
        store(ontologyTerms.values());
        store(ontologyAnnotations);
	store(expressionValues);
    }

    /**
     * Process the file which describes the samples.
     *
     * Only the first two columns are required.
     * 0            1     2            3          4       5                  6        7         8                9          10
     * #identifier  name  description  treatment  tissue  development_stage  species  genotype  replicate_group  biosample  sra_experiment
     */
    void processSamplesFile() throws IOException {
        // map the optional header attributes to the corresponding ExpressionSample attributes
        // description  treatment  tissue  development_stage  species  genotype  replicate_group  biosample  sra_experiment
        HashMap<String,String> optional = new HashMap<>();
        optional.put("description", "description");
        optional.put("treatment", "treatment");
        optional.put("tissue", "tissue");
        optional.put("development_stage", "developmentStage");
        optional.put("species", "species");
        optional.put("genotype", "genotype");
        optional.put("replicate_group", "replicateGroup");
        optional.put("biosample", "bioSample");
        optional.put("sra_experiment", "sraExperiment");
        String[] colnames = null;
        int num = 0;
	BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#identifier")) {
                // use the column header names to identify column content
                colnames = line.split("\t");
            } else if (line.startsWith("#") || line.trim().length()==0) {
                // comment or blank
                continue;
            } else {
                // sample data row, increment num as we go down the file
                num++;
                String[] parts = line.split("\t");
                String identifier = parts[0];
                String name = parts[1];
                Item sample = getSample(identifier);
                sample.setAttribute("name", name);
                sample.setAttribute("num", String.valueOf(num));
                // optional attributes are given by column names
                Map<String,String> attributes = new HashMap<>();
                for (int i=0; i<colnames.length; i++) {
                    attributes.put(colnames[i], parts[i]);
                }
                for (String key : optional.keySet()) {
                    if (attributes.get(key)!=null && attributes.get(key).trim().length()>0) {
                        sample.setAttribute(optional.get(key), attributes.get(key));
                    }
                }
            }
        }
	br.close();
    }

    /**
     * Process a gene expression file. Each gene-sample entry creates an ExpressionValue.
     *
     * Header line must start with gene_id with tab-separated sample identifiers.
     * (If not, NumberFormatException will occur to let you know!)
     *
     * gene_id	SRR5199304	SRR5199305	SRR5199306	SRR5199307	...
     * cajca.ICPL87119.gnm1.ann1.C.cajan_00002	0	0	0	0	...
     */
    void processValuesFile() throws IOException {
        List<Item> sampleList = new ArrayList<>();
	BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            if (parts[0].equals("gene_id")) {
                // header line gives samples in order so add to list
                for (int i=1; i<parts.length; i++) {
                    sampleList.add(getSample(parts[i]));
                }
            } else {
                // a gene expression values line
                Item gene = getGene(parts[0]);
                for (int i=1; i<parts.length; i++) {
                    double value = Double.parseDouble(parts[i]);
                    Item sample = sampleList.get(i-1);
                    Item expressionValue = createItem("ExpressionValue");
                    expressionValue.setAttribute("value", String.valueOf(value));
                    expressionValue.setReference("feature", gene);
                    expressionValue.setReference("sample", sample);
		    expressionValues.add(expressionValue);
                }
            }
        }
	br.close();
    }

    /**
     * Process a file relating samples to ontology terms.
     */
    void processOboFile() throws IOException {
	BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            String sampleId = parts[0];
            String termId = parts[1];
            Item sample = getSample(sampleId);
            Item ontologyTerm = getOntologyTerm(termId);
            Item ontologyAnnotation = createItem("OntologyAnnotation");
            ontologyAnnotation.setReference("subject", sample);
            ontologyAnnotation.setReference("ontologyTerm", ontologyTerm);
            ontologyAnnotations.add(ontologyAnnotation);
        }
    }

    /**
     * Get or create a Sample.
     */
    Item getSample(String id) {
        if (samples.containsKey(id)) {
            return samples.get(id);
        } else {
            Item sample = createItem("ExpressionSample");
            sample.setAttribute("primaryIdentifier", id);
            samples.put(id, sample);
            return sample;
        }
    }

    /**
     * Get or create an OntologyTerm.
     */
    Item getOntologyTerm(String id) {
        if (ontologyTerms.containsKey(id)) {
            return ontologyTerms.get(id);
        } else {
            Item ontologyTerm = createItem("OntologyTerm");
            ontologyTerm.setAttribute("identifier", id);
            ontologyTerms.put(id, ontologyTerm);
            return ontologyTerm;
        }
    }

    /**
     * Get or create a Gene.
     */
    Item getGene(String id) {
        if (genes.containsKey(id)) {
            return genes.get(id);
        } else {
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", id);
            genes.put(id, gene);
            return gene;
        }
    }



}
