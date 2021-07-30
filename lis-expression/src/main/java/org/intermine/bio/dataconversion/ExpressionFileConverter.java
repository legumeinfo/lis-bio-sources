package org.intermine.bio.dataconversion;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
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

/**
 * DataConverter to create ExpressionSource, ExpressionSample and ExpressionValue items from a single set of datastore expression files.
 *
 * expression/G19833.gnm1.ann1.expr.4ZDQ/
 * ├── phavu.G19833.gnm1.ann1.expr.4ZDQ.obo.tsv
 * ├── phavu.G19833.gnm1.ann1.expr.4ZDQ.samples.tsv
 * ├── phavu.G19833.gnm1.ann1.expr.4ZDQ.values.tsv
 * └── README.G19833.gnm1.ann1.expr.4ZDQ.yml
 *
 * @author Sam Hokin
 */
public class ExpressionFileConverter extends DatastoreFileConverter {
    
    private static final Logger LOG = Logger.getLogger(ExpressionFileConverter.class);

    // singletons
    Item expressionSource = createItem("ExpressionSource");
    Item bioProject;
    Item publication;
    String expressionUnit = "TPM";

    // Lists
    List<Item> expressionValues = new LinkedList<>();
    List<Item> annotations = new LinkedList<>();

    // Maps
    Map<String,Item> terms = new HashMap<>();
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
    }

    /**
     * Called for each file found.
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws IOException {
	if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
        } else if (getCurrentFile().getName().endsWith("samples.tsv")) {
	    processSamples(reader);
        } else if (getCurrentFile().getName().endsWith("values.tsv")) {
            processExpression(reader);
        } else if (getCurrentFile().getName().endsWith("obo.tsv")) {
            processOboFile(reader);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        // DatastoreFileConverter items
        store(organisms.values());
        store(strains.values());
        store(dataSource);
        store(dataSets.values());
        // local items
        store(expressionSource);
	store(publication);
        if (bioProject!=null) store(bioProject);
	store(genes.values());
        store(samples.values());
        store(terms.values());
        store(annotations);
	store(expressionValues);
    }

    /**
     * Process the README, which contains the expression source metadata.
     *
     * identifier: 4ZDQ
     * scientific_name: Phaseolus vulgaris
     * taxid: 3885
     * bioproject: PRJNA210619
     * scientific_name_abbrev: phavu
     * genotype:
     * - G19833
     * description: "Phaseolus vulgaris cv. Negro jamapa expression profiles in the P. vulgaris reference genome (G 19833) for 24 unique samples..."
     * publication_doi: 10.1186/1471-2164-15-866
     * publication_title:  "An RNA-Seq based gene expression atlas of the common bean"
     * contributors: "O'Rourke, Jamie A; Iniguez, Luis P; Fu, Fengli; ..."
     * expression_unit: TPM
     * sraproject: SRP046307
     * geoseries: GSEnnnnn
     */
    void processReadme(Reader reader) throws IOException {
        Readme readme = Readme.getReadme(reader);
        // check required stuff
        if (readme.identifier==null ||
            readme.taxid==null ||
            readme.synopsis==null ||
            readme.description==null ||
            readme.genotype==null ||
            readme.publication_doi==null ||
            readme.publication_title==null ||
            readme.expression_unit==null) {
            throw new RuntimeException("ERROR: a required field is missing from "+getCurrentFile().getName()+": "+
                                       "Required fields are: identifier, taxid, synopsis, description, genotype, publication_doi, publication_title, expression_unit");
        }
        // Organism from README taxid rather than filename
        Item organism = getOrganism(Integer.parseInt(readme.taxid));
        Item strain = getStrain(readme.genotype[0], organism);
        // ExpressionSource
        expressionSource.setAttribute("primaryIdentifier", readme.identifier);
        expressionSource.setAttribute("synopsis", readme.synopsis);
        expressionSource.setAttribute("description", readme.description);
        expressionSource.setReference("organism", organism);
        expressionSource.setReference("strain", strain);
        expressionUnit = readme.expression_unit; // applied to ExpressionValue.unit in close()
        if (readme.geoseries!=null) expressionSource.setAttribute("geoSeries", readme.geoseries);
        if (readme.sraproject!=null) expressionSource.setAttribute("sra", readme.sraproject);
        // BioProject
        if (readme.bioproject!=null) {
            bioProject = createItem("BioProject");
            bioProject.setAttribute("accession", readme.bioproject);
            expressionSource.setReference("bioProject", bioProject);
        }
        // Publication
        publication = createItem("Publication");
        publication.setAttribute("doi", readme.publication_doi);
        publication.setAttribute("title", readme.publication_title);
        expressionSource.addToCollection("publications", publication);
        // add to samples if we've already read the samples file
        if (samples.size()>0) {
            for (Item sample : samples.values()) {
                sample.addToCollection("publications", publication);
            }
        }
        
        // override DataSet.description from README
        Item dataSet = getDataSet();
        dataSet.setAttribute("description", readme.description);
        dataSet.setReference("publication", publication);
        expressionSource.setReference("dataSet", dataSet);
    }

    /**
     * Process the file which describes the samples.
     *
     * cajca.ICPL87119.gnm1.ann1.expr.KEY4.samples.tsv
     * 
     * 0sample_name                       1key       sample_uniquename                  description                                    treatment             tissue
     * Mature seed at reprod (SRR5199304) SRR5199304 Mature seed at reprod (SRR5199304) Mature seed at Reproductive stage (SRR5199304) Mature seed at reprod Mature seed 
     *
     * dev_stage          age      organism      infraspecies cultivar        other sra_run    biosample_accession sra_accession bioproject_accession sra_study
     * Reproductive stage          Cajanus cajan ICPL87119    Asha(ICPL87119)       SRR5199304 SAMN06264156        SRS1937936    PRJNA354681          SRP097728
     */
    void processSamples(Reader reader) throws IOException {
        Item organism = getOrganism();
        Item strain = getStrain(organism);
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
		// common references
		sample.setReference("organism", organism);
		sample.setReference("strain", strain);
		sample.setReference("source", expressionSource);
                // add publication if we've already read the README
                if (publication!=null) sample.addToCollection("publications", publication);
		// sample-specific attributes
                for (int i=0; i<colnames.length; i++) {
                    switch(colnames[i]) {
                    case "sample_uniquename" :
                        // ignored
                        break;
                    case "description" :
                        sample.setAttribute("description", parts[i]);
                        break;
                    case "treatment" :
                        // ignored
                        break;
                    case "tissue" :
                        // ignored
                        break;
                    case "dev_stage" :
                        // ignored
                        break;
                    case "age" :
                        // ignored
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
                        // ignored
                        break;
                    case "sra_run" :
                        // ignored
                        break;
                    case "biosample_accession" :
                        // bioSample accession is all we need, everything else is on NCBI
                        sample.setAttribute("bioSample", parts[i]);
                        break;
                    case "sra_accession" :
                        // ignored
                        break;
                    case "bioproject_accession" :
                        // ignored
                        break;
                    case "sra_study" :
                        // ignored
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
        Item organism = getOrganism();
        Item strain = getStrain(organism);
        Item dataSet = getDataSet();
        List<Item> sampleList = new LinkedList<>();
	BufferedReader br = new BufferedReader(reader);
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
		// common references
		gene.addToCollection("dataSets", dataSet);
		gene.setReference("organism", organism);
		gene.setReference("strain", strain);
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
    void processOboFile(Reader reader) throws IOException {
	BufferedReader br = new BufferedReader(reader);
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
