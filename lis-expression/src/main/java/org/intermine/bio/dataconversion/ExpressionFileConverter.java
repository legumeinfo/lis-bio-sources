package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.io.IOException;
import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.ncgr.datastore.Readme;
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
    Item bioProject;
    String expressionUnit;

    // Lists
    List<Item> expressionValues = new LinkedList<>();
    List<Item> annotations = new LinkedList<>();

    // Maps
    Map<String,Item> terms = new HashMap<>();
    Map<String,Item> samples = new HashMap<>();
    Map<String,Item> genes = new HashMap<>();

    // keep track of files read in case they're gzipped!
    boolean samplesRead = false;
    boolean valuesRead = false;
    boolean oboRead = false;
        
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
            expressionSource.setReference("organism", organism);
            expressionSource.setReference("strain", strain);
            expressionUnit = readme.expression_unit; // applied to ExpressionValue.unit in close()
            if (readme.geoseries!=null) expressionSource.setAttribute("geoSeries", readme.geoseries);
            if (readme.sraproject!=null) expressionSource.setAttribute("sra", readme.sraproject);
            if (readme.bioproject!=null) {
                bioProject = createItem("BioProject");
                bioProject.setAttribute("accession", readme.bioproject);
                expressionSource.setReference("bioProject", bioProject);
            }
            if (publication!=null) expressionSource.addToCollection("publications", publication);
        } else if (getCurrentFile().getName().endsWith("samples.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
	    processSamples();
            samplesRead = true;
        } else if (getCurrentFile().getName().endsWith("values.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processExpression();
            valuesRead = true;
        } else if (getCurrentFile().getName().endsWith("obo.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processOboFile();
            oboRead = true;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        if (!samplesRead || !valuesRead || !oboRead) {
            throw new RuntimeException("One of samples.tsv.gz, values.tsv.gz, and/or obo.tsv.gz files not read. Aborting.");
        }
        // DatastoreFileConverter items
        storeCollectionItems();
        // add references to samples
        for (Item sample : samples.values()) {
            sample.setReference("organism", organism);
            sample.setReference("strain", strain);
            sample.setReference("source", expressionSource);
            if (publication!=null) sample.addToCollection("publications", publication);
        }
        // add references to genes
        for (Item gene : genes.values()) {
            gene.setReference("organism", organism);
            gene.setReference("strain", strain);
        }
        // local items
        store(expressionSource);
        if (bioProject!=null) store(bioProject);
	store(genes.values());
        store(samples.values());
        store(terms.values());
        store(annotations);
	store(expressionValues);
    }

    /**
     * Process the file which describes the samples.
     *
     * cajca.ICPL87119.gnm1.ann1.expr.KEY4.samples.tsv.gz
     * 
     * 0sample_name                       1key       sample_uniquename                  description                                    treatment             tissue
     * Mature seed at reprod (SRR5199304) SRR5199304 Mature seed at reprod (SRR5199304) Mature seed at Reproductive stage (SRR5199304) Mature seed at reprod Mature seed 
     *
     * dev_stage          age      organism      infraspecies cultivar        other sra_run    biosample_accession sra_accession bioproject_accession sra_study
     * Reproductive stage          Cajanus cajan ICPL87119    Asha(ICPL87119)       SRR5199304 SAMN06264156        SRS1937936    PRJNA354681          SRP097728
     */
    void processSamples() throws IOException {
        String[] colnames = null;
        int num = 0;
	BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            // header contains colnames
            if (line.startsWith("sample_name")) {
                colnames = parts;
            } else {
                // increment sample number
                num++;
                // 0sample_name = Sample.name
                // 1key = Sample.primaryIdentifier
                String name = parts[0];
                String id = parts[1];
                Item sample = samples.get(id);
                if (sample==null) {
                    sample = createItem("ExpressionSample");
                    sample.setAttribute("primaryIdentifier", id);
                    samples.put(id, sample);
                }
                sample.setAttribute("name", name);
                sample.setAttribute("num", String.valueOf(num));
		// sample-specific attributes
                for (int i=0; i<colnames.length; i++) {
                    switch(colnames[i]) {
                    case "description" :
                        sample.setAttribute("description", parts[i]);
                        break;
                    case "biosample_accession" :
                        // bioSample accession is all we need, everything else is on NCBI
                        sample.setAttribute("bioSample", parts[i]);
                        break;
                    }
                }
            }
        }
	br.close();
    }

    /**
     * Process a gene expression file. Each gene-sample entry creates an ExpressionValue.
     *
     * cajca.ICPL87119.gnm1.ann1.expr.KEY4.values.tsv
     * --------------------------------------------------
     * geneID	SRR5199304	SRR5199305	SRR5199306	SRR5199307	...
     * cajca.ICPL87119.gnm1.ann1.C.cajan_00002	0	0	0	0	...
     */
    void processExpression() throws IOException {
        List<Item> sampleList = new LinkedList<>();
	BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            if (parts[0].toLowerCase().equals("geneid")) {
                // header line gives samples in order so add to list
                for (int i=1; i<parts.length; i++) {
                    String sampleId = parts[i];
                    Item sample = samples.get(sampleId);
                    if (sample==null) {
                        sample = createItem("ExpressionSample");
                        sample.setAttribute("primaryIdentifier", sampleId);
                        samples.put(sampleId, sample);
                    }
                    sampleList.add(sample);
                }
            } else {
                // a gene line
                String geneId = parts[0];
                Item gene = createItem("Gene");
                gene.setAttribute("primaryIdentifier", geneId);
                genes.put(geneId, gene);
                for (int i=1; i<parts.length; i++) {
                    double value = Double.parseDouble(parts[i]);
                    Item sample = sampleList.get(i-1);
                    Item expressionValue = createItem("ExpressionValue");
                    expressionValue.setAttribute("unit", expressionUnit);
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
            Item sample = samples.get(sampleId);
            if (sample==null) {
                sample = createItem("ExpressionSample");
                sample.setAttribute("primaryIdentifier", sampleId);
                samples.put(sampleId, sample);
            }
            Item term = terms.get(termId);
            if (term==null) {
                term = createItem("OntologyTerm");
                term.setAttribute("identifier", termId);
                terms.put(termId, term);
            }
            Item annotation = createItem("OntologyAnnotation");
            annotation.setReference("subject", sample);
            annotation.setReference("ontologyTerm", term);
            annotations.add(annotation);
        }
    }
}
