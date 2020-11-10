package org.intermine.bio.dataconversion;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.Reader;
import java.io.IOException;
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
 * DataConverter to create ExpressionSource, ExpressionSample and ExpressionValue items from a single set of datastore expression files.
 *
 * strain.gnmN.annN.expr.KEY4
 * ├── gensp.strain.gnmN.annN.expr.KEY4.samples.tsv
 * ├── gensp.strain.gnmN.annN.expr.KEY4.source.tsv
 * └── gensp.strain.gnmN.annN.expr.KEY4.values.tsv
 *
 * @author Sam Hokin
 */
public class ExpressionFileConverter extends DatastoreFileConverter {
    
    private static final Logger LOG = Logger.getLogger(ExpressionFileConverter.class);

    Item organism;
    Item strain;
    Item expressionSource;
    Item bioProject;
    Item publication;
    String unit;
    List<Item> expressionValues = new ArrayList<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> samples = new HashMap<>();
    Map<String,Item> genes = new HashMap<>();

    /**
     * Constructor.
     *
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     * @throws ObjectStoreException os
     */
    public ExpressionFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
        expressionSource = createItem("ExpressionSource");
        bioProject = createItem("BioProject");
        publication = createItem("Publication");
        expressionSource.setReference("publication", publication);
        expressionSource.setReference("bioProject", bioProject);
    }

    /**
     * Called for each file found.
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws IOException {
	if (getCurrentFile().getName().endsWith("source.tsv")) {
            processSource(reader);
        } else if (getCurrentFile().getName().endsWith("samples.tsv")) {
	    processSamples(reader);
        } else if (getCurrentFile().getName().endsWith("values.tsv")) {
            processExpression(reader);
        } else if (getCurrentFile().getName().endsWith("ontology.tsv")) {
            processOntology(reader);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        // DatastoreFileConverter items
        store(dataSource);
        store(dataSets.values());
        store(organisms.values());
        store(strains.values());
        // local items
        store(expressionSource);
	store(publication);
        store(bioProject);
	store(genes.values());
        store(samples.values());
        store(ontologyTerms.values());
	store(expressionValues);
    }

    /**
     * Process the datasource meta data file and put the info into a DataSet.
     *
     * cajca.ICPL87119.gnm1.ann1.expr.KEY4.source.tsv
     *
     * NAME	Gene expression atlas of pigeonpea Asha(ICPL87119)
     * SHORTNAME	Pigeonpea gene expression atlas
     * ORIGIN	SRA
     * DESCRIPTION	To be able to link the genome sequence information to the phenotype, especially...
     * UNIT     TPM
     * BIOPROJ_ACC	PRJNA354681
     * BIOPROJ_TITLE	Gene expression atlas of pigeonpea
     * BIOPROJ_DESCRIPTION	Pigeonpea (Cajanus cajan) is an important grain legume of the semi-arid tropics, mainly used...
     * BIOPROJ_PUBLICATION	Pazhamala, L. T., Purohit, S., Saxena, R. K., Garg, V., Krishnamurthy, L., Verdier, J., & Varshney, R. K. (2017). Gene expression...
     * SRA_PROJ_ACC	SRP097728
     * GEO_SERIES	
     * PUB_PMID 123456
     * PUB_LINK	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5429002/
     * PUB_FULLLINK	https://academic.oup.com/jxb/article/68/8/2037/3051749/Gene-expression-atlas-of-pigeonpea-and-its
     */
    void processSource(Reader reader) throws IOException {
        organism = getOrganism();
        strain = getStrain(organism);
        Item dataSet = getDataSet();
        expressionSource.setReference("dataSet", dataSet);
        unit = "TPM"; // default
	BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            if (parts.length==2) {
                switch(parts[0]) {
                case "NAME" :
                    expressionSource.setAttribute("identifier", parts[1]);
                    break;
                case "SHORTNAME" :
                    expressionSource.setAttribute("shortName", parts[1]);
                    break;
                case "ORIGIN" :
                    expressionSource.setAttribute("origin", parts[1]);
                    break;
                case "DESCRIPTION" :
                    expressionSource.setAttribute("description", parts[1]);
                    break;
                case "UNIT" :
                    unit = parts[1];
                case "GEO_SERIES" :
                    expressionSource.setAttribute("geoSeries", parts[1]);
                    break;
                case "SRA_PROJ_ACC" :
                    expressionSource.setAttribute("sra", parts[1]);
                    break;
		case "PUB_PMID" :
                    publication.setAttribute("pubMedId", parts[1]);
		    break;
                case "PUB_LINK" :
                    publication.setAttribute("url", parts[1]);
                    break;
                case "PUB_FULLLINK" :
                    publication.setAttribute("url", parts[1]);
                    break;
                case "BIOPROJ_ACC" :
                    bioProject.setAttribute("accession", parts[1]);
                    break;
                case "BIOPROJ_TITLE" :
                    bioProject.setAttribute("title", parts[1]);
                    break;
                case "BIOPROJ_DESCRIPTION" :
                    bioProject.setAttribute("description", parts[1]);
                    break;
                case "BIOPROJ_PUBLICATION" :
                    // do nothing
                    break;
                default :
                    // log existence of unknown field
                    LOG.info(getCurrentFile().getName()+" contains unknown field:"+line);
                }
            }
        }
	br.close();
    }

    /**
     * Process the file which describes the samples.
     *
     * cajca.ICPL87119.gnm1.ann1.expr.KEY4.samples.tsv
     * 
     * sample_name                        key        sample_uniquename                  description                                    treatment             tissue
     * Mature seed at reprod (SRR5199304) SRR5199304 Mature seed at reprod (SRR5199304) Mature seed at Reproductive stage (SRR5199304) Mature seed at reprod Mature seed 
     *
     * dev_stage          age      organism      infraspecies cultivar        other sra_run    biosample_accession sra_accession bioproject_accession sra_study
     * Reproductive stage          Cajanus cajan ICPL87119    Asha(ICPL87119)       SRR5199304 SAMN06264156        SRS1937936    PRJNA354681          SRP097728
     */
    void processSamples(Reader reader) throws IOException {
        organism = getOrganism();
        strain = getStrain(organism);
        Item dataSet = getDataSet();
        String[] colnames = null;
        int num = 0;
	BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            // header contains colnames
            if (line.startsWith("sample_name")) {
                colnames = parts;
            } else if (colnames!=null) {
                if (parts.length!=colnames.length) {
                    System.err.println(line);
                    throw new RuntimeException("ERROR: colnames.length="+colnames.length+" but line has "+parts.length+" fields!");
                }
                // get/create the sample with primaryIdentifier in column 2
                String id = parts[1];
                Item sample = samples.get(id);
                if (sample==null) {
                    sample = createItem("ExpressionSample");
                    sample.setAttribute("primaryIdentifier", id);
                }
		// common references
		sample.setReference("organism", organism);
		sample.setReference("strain", strain);
		sample.setReference("source", expressionSource);
		sample.addToCollection("dataSets", dataSet);
                sample.addToCollection("publications", publication);
		// sample-specific attributes
                num++;
                sample.setAttribute("num", String.valueOf(num));
                for (int i=0; i<colnames.length; i++) {
                    switch(colnames[i]) {
                    case "sample_name" :
                        sample.setAttribute("name", parts[i]);
                        break;
                    case "sample_uniquename" :
                        break;
                    case "description" :
                        sample.setAttribute("description", parts[i]);
                        break;
                    case "treatment" :
                        break;
                    case "tissue" :
                        break;
                    case "dev_stage" :
                        break;
                    case "age" :
                        break;
                    case "organism" :
                        // we use the filename for organism
                        break;
                    case "infraspecies" :
                        // we use the filename for strain
                        break;
                    case "cultivar" :
                        // we use the filename for strain
                        break;
                    case "other" :
                        break;
                    case "sra_run" :
                        break;
                    case "biosample_accession" :
                        // bioSample accession is all we need, everything else is on NCBI
                        sample.setAttribute("bioSample", parts[i]);
                        break;
                    case "sra_accession" :
                        break;
                    case "bioproject_accession" :
                        break;
                    case "sra_study" :
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
    void processExpression(Reader reader) throws IOException {
        organism = getOrganism();
        strain = getStrain(organism);
        Item dataSet = getDataSet();
        List<Item> sampleList = new ArrayList<>();
	BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            if (parts[0].toLowerCase().equals("geneid")) {
                // header line gives samples in order so add to list
                for (int i=1; i<parts.length; i++) {
                    String sampleId = parts[i];
                    if (samples.containsKey(sampleId)) {
                        sampleList.add(samples.get(sampleId));
                    } else {
                        Item sample = createItem("ExpressionSample");
			// common references
			sample.setReference("source", expressionSource);
			sample.setReference("dataSet", dataSet);
			sample.setReference("organism", organism);
			sample.setReference("strain", strain);
                        sample.setAttribute("identifier", sampleId);
                        samples.put(sampleId, sample);
                        sampleList.add(sample);
                    }
                }
            } else {
                // a gene line
                String geneId = parts[0];
                Item gene = createItem("Gene");
                gene.setAttribute("primaryIdentifier", geneId);
		// common references
		gene.addToCollection("dataSets", dataSet);
		gene.setReference("organism", organism);
		gene.setReference("strain", strain);
                genes.put(geneId, gene);
                for (int i=1; i<parts.length; i++) {
                    double value = Double.parseDouble(parts[i]);
                    Item sample = sampleList.get(i-1);
                    Item expressionValue = createItem("ExpressionValue");
                    expressionValue.setAttribute("value", String.valueOf(value));
                    expressionValue.setAttribute("unit", unit);
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
    void processOntology(Reader reader) throws IOException {
	BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            String sampleId = parts[0];
            String ontologyTermId = parts[1];
            Item sample = samples.get(sampleId);
            if (sample==null) {
                sample = createItem("ExpressionSample");
                sample.setAttribute("primaryIdentifier", sampleId);
                samples.put(sampleId, sample);
            }
            Item ontologyTerm = ontologyTerms.get(ontologyTermId);
            if (ontologyTerm==null) {
                ontologyTerm = createItem("OntologyTerm");
                ontologyTerm.setAttribute("identifier", ontologyTermId);
                ontologyTerms.put(ontologyTermId, ontologyTerm);
            }
            Item ontologyAnnotation = createItem("OntologyAnnotation");
            ontologyAnnotation.setReference("subject", sample);
            ontologyAnnotation.setReference("ontologyTerm", ontologyTerm);
            try {
                store(ontologyAnnotation);
            } catch (ObjectStoreException e) {
                throw new RuntimeException(e);
            }
        }
    }
}
