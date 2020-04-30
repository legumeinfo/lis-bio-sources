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
 * ICPL87119.gnm1.ann1.expr.KEY4
 * ├── cajca.ICPL87119.gnm1.ann1.KEY4.exprSamples.tsv
 * ├── cajca.ICPL87119.gnm1.ann1.KEY4.exprSource.tsv
 * └── cajca.ICPL87119.gnm1.ann1.KEY4.genesSamplesTpm.tsv
 *
 * cajca.ICPL87119.gnm1.ann1.KEY4.exprSource.tsv
 * ---------------------------------------------
 * #DATASOURCE	
 * ##	
 * ##Name, descriptive (max 250 char)	
 * NAME	Gene expression atlas of pigeonpea Asha(ICPL87119)
 * #Shortname: max 64 char	
 * SHORTNAME	Pigeonpea gene expression atlas
 * #Origin: ex. SRA, GEO, A labname, etc.	
 * ORIGIN	SRA
 * #Description: our description if needed	
 * DESCRIPTION	To be able to link the genome sequence information to the phenotype, especially...
 * #NCBI BioProj acc, PRJNA num	
 * BIOPROJ_ACC	PRJNA354681
 * #NCBI BioProj Title	
 * BIOPROJ_TITLE	Gene expression atlas of pigeonpea
 * #NCBI BioProj description/ abstract(?)	
 * BIOPROJ_DESCRIPTION	Pigeonpea (Cajanus cajan) is an important grain legume of the semi-arid tropics, mainly used...
 * #NCBI SRP number	
 * SRA_PROJ_ACC	SRP097728
 * # GEO series if exists	
 * GEO_SERIES	
 * #Associated publication at NCBI	
 * BIOPROJ_PUBLICATION	Pazhamala, L. T., Purohit, S., Saxena, R. K., Garg, V., Krishnamurthy, L., Verdier, J., & Varshney, R. K. (2017). Gene expression...
 * # PubMed ID
 * PUB_PMID     123456
 * #Link to Pub	
 * PUB_LINK	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5429002/
 * #Link to full publication	
 * PUB_FULLLINK	https://academic.oup.com/jxb/article/68/8/2037/3051749/Gene-expression-atlas-of-pigeonpea-and-its
 
 * cajca.ICPL87119.gnm1.ann1.KEY4.exprSamples.tsv
 * ----------------------------------------------
 * #SAMPLE																
 * #																
 * #Hints:																
 * #sample display name	Lookup key for data column (=SRR no)	unique name for LIS Chado	Descriptive	Exp design: factors-rep combinations	tissue/organ/plant part			Chado spelling	sub-spp, cv-type	strain, line, cultivar, genotype	Other_attributes, place-holder, add extra columns	SRR(not in recc)	SAMN number	SRS number	PRJN number	SRP number
 * sample_name	key	sample_uniquename	description	treatment	tissue	dev_stage	age	organism	infraspecies	cultivar	other	sra_run	biosample_accession	sra_accession	bioproject_accession	sra_study
 * Mature seed at reprod (SRR5199304)	SRR5199304	Mature seed at reprod (SRR5199304)	Mature seed at Reproductive stage (SRR5199304)	Mature seed at reprod	Mature seed	Reproductive stage	Cajanus cajan	ICPL87119	Asha(ICPL87119)		SRR5199304	SAMN06264156	SRS1937936	PRJNA354681	SRP097728
 *
 * cajca.ICPL87119.gnm1.ann1.KEY4.genesSamplesTpm.tsv
 * --------------------------------------------------
 * geneID	SRR5199304	SRR5199305	SRR5199306	SRR5199307	...
 * cajca.ICPL87119.gnm1.ann1.C.cajan_00002	0	0	0	0	...
 *
 * @author Sam Hokin
 */
public class ExpressionFileConverter extends DatastoreFileConverter {
    
    private static final Logger LOG = Logger.getLogger(ExpressionFileConverter.class);

    // local Items to store
    Item expressionSource;
    Item bioProject;
    List<Item> publications = new ArrayList<>();
    List<Item> expressionValues = new ArrayList<>();
    Map<String,Item> samples = new HashMap<>();
    Map<String,Item> genes = new HashMap<>();

    // we have to process three files in order, so get the Files from the Readers
    File sourceFile;
    File samplesFile;
    File expressionFile;

    /**
     * Constructor.
     *
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     * @throws ObjectStoreException os
     */
    public ExpressionFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
	dataSource = getDataSource();
    }

    /**
     * Called for each file found.
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().contains("README")) return;
	if (getCurrentFile().getName().endsWith("exprSource.tsv")) {
            // cajca.ICPL87119.gnm1.ann1.KEY4.exprSource.tsv
	    sourceFile = getCurrentFile();
        } else if (getCurrentFile().getName().endsWith("exprSamples.tsv")) {
            // cajca.ICPL87119.gnm1.ann1.KEY4.exprSamples.tsv
	    samplesFile = getCurrentFile();
        } else if (getCurrentFile().getName().endsWith("genesSamplesTpm.tsv")) {
            // cajca.ICPL87119.gnm1.ann1.KEY4.genesSamplesTpm.tsv
	    expressionFile = getCurrentFile();
        }
	// process once we have all three Files instantiated
	if (sourceFile!=null && samplesFile!=null && expressionFile!=null) {
            processSource(sourceFile);
	    processSamples(samplesFile);
            processExpression(expressionFile);
	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        // store everything that is global
        store(dataSource);
        store(dataSets.values());
        store(organisms.values());
        store(strains.values());
        store(bioProject);
        store(expressionSource);
	store(genes.values());
        store(samples.values());
	store(publications);
	store(expressionValues);
    }

    /**
     * Process the datasource meta data file and put the info into a DataSet.
     *
     * cajca.ICPL87119.gnm1.ann1.KEY4.exprSource.tsv
     *
     * NAME	Gene expression atlas of pigeonpea Asha(ICPL87119)
     * SHORTNAME	Pigeonpea gene expression atlas
     * ORIGIN	SRA
     * DESCRIPTION	To be able to link the genome sequence information to the phenotype, especially...
     * BIOPROJ_ACC	PRJNA354681
     * BIOPROJ_TITLE	Gene expression atlas of pigeonpea
     * BIOPROJ_DESCRIPTION	Pigeonpea (Cajanus cajan) is an important grain legume of the semi-arid tropics, mainly used...
     * SRA_PROJ_ACC	SRP097728
     * GEO_SERIES	
     * BIOPROJ_PUBLICATION	Pazhamala, L. T., Purohit, S., Saxena, R. K., Garg, V., Krishnamurthy, L., Verdier, J., & Varshney, R. K. (2017). Gene expression...
     * PUB_PMID 123456
     * PUB_LINK	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5429002/
     * PUB_FULLLINK	https://academic.oup.com/jxb/article/68/8/2037/3051749/Gene-expression-atlas-of-pigeonpea-and-its
     */
    void processSource(File file) throws IOException {
	BufferedReader br = new BufferedReader(new FileReader(file));
	Item dataSet = getDataSet(); // populates name (=filename), description and url (from project.xml)
	Item organism = getOrganism();
	Item strain = getStrain(organism);
        expressionSource = createItem("ExpressionSource");
        expressionSource.setAttribute("unit", "TPM"); // NOTE: assume TPM
        expressionSource.setReference("dataSet", dataSet);
        // we'll favor pubMedId
	String pubMedId = null;
        String pubLink = null;
        String pubFullLink = null;
        // create BioProject only if BIOPROJ lines exist
        bioProject = null;
        // spin through the file
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            if (parts.length==2) {
                switch(parts[0]) {
                case "NAME" :
                    // dataSet.setAttribute("name", parts[1]); // use the file name, not what's in the source file
                    expressionSource.setAttribute("primaryIdentifier", parts[1]);
                    break;
                case "SHORTNAME" :
                    dataSet.setAttribute("shortName", parts[1]);
                    break;
                case "ORIGIN" :
                    dataSet.setAttribute("origin", parts[1]);
                    break;
                case "DESCRIPTION" :
                    // dataSet.setAttribute("description", parts[1]); // use what's in project.xml
                    break;
                case "GEO_SERIES" :
                    dataSet.setAttribute("geoSeries", parts[1]);
                    break;
                case "SRA_PROJ_ACC" :
                    dataSet.setAttribute("sra", parts[1]);
                    break;
		case "PUB_PMID" :
		    pubMedId = parts[1];
		    break;
                case "PUB_LINK" :
                    pubLink = parts[1];
                    break;
                case "PUB_FULLLINK" :
                    pubFullLink = parts[1];
                    break;
                case "BIOPROJ_ACC" :
                    if (bioProject==null) bioProject = createItem("BioProject");
                    bioProject.setAttribute("accession", parts[1]);
                    break;
                case "BIOPROJ_TITLE" :
                    if (bioProject==null) bioProject = createItem("BioProject");
                    bioProject.setAttribute("title", parts[1]);
                    break;
                case "BIOPROJ_DESCRIPTION" :
                    if (bioProject==null) bioProject = createItem("BioProject");
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
        if (pubMedId!=null || pubLink!=null || pubFullLink!=null) {
            Item pub = createItem("Publication");
	    if (pubMedId!=null) {
		pub.setAttribute("pubMedId", pubMedId);
	    } else if (pubLink!=null) {
                pub.setAttribute("url", pubLink);
            } else if (pubFullLink!=null) {
                pub.setAttribute("url", pubFullLink);
            }
            dataSet.setReference("publication", pub);
	    publications.add(pub);
        }
        if (bioProject!=null) {
            dataSet.setReference("bioProject", bioProject);
        }
	br.close();
    }

    /**
     * Process the file which describes the samples.
     *
     * cajca.ICPL87119.gnm1.ann1.KEY4.exprSamples.tsv
     * 
     * sample_name                        key        sample_uniquename                  description                                    treatment             tissue
     * Mature seed at reprod (SRR5199304) SRR5199304 Mature seed at reprod (SRR5199304) Mature seed at Reproductive stage (SRR5199304) Mature seed at reprod Mature seed 
     *
     * dev_stage          age      organism      infraspecies cultivar        other sra_run    biosample_accession sra_accession bioproject_accession sra_study
     * Reproductive stage          Cajanus cajan ICPL87119    Asha(ICPL87119)       SRR5199304 SAMN06264156        SRS1937936    PRJNA354681          SRP097728
     */
    void processSamples(File file) throws IOException {
	BufferedReader br = new BufferedReader(new FileReader(file));
	Item dataSet = getDataSet();
	Item organism = getOrganism();
	Item strain = getStrain(organism);
        String[] colnames = null;
        int num = 0;
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
                String key = null;
                Item sample = createItem("ExpressionSample");
		// common references
		sample.setReference("source", expressionSource);
		sample.setReference("dataSet", dataSet);
		sample.setReference("organism", organism);
		sample.setReference("strain", strain);
		// sample-specific attributes
                num++;
                sample.setAttribute("num", String.valueOf(num));
                for (int i=0; i<colnames.length; i++) {
                    switch(colnames[i]) {
                    case "sample_name" :
                        sample.setAttribute("name", parts[i]);
                        break;
                    case "key" :
                        key = parts[i];
                        sample.setAttribute("primaryIdentifier", parts[i]);
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
                if (key!=null && sample!=null) {
                    samples.put(key, sample);
                }
            }
        }
	br.close();
    }

    /**
     * Process a gene expression file. Each gene-sample entry creates an ExpressionValue.
     *
     * cajca.ICPL87119.gnm1.ann1.KEY4.genesSamplesTpm.tsv
     * --------------------------------------------------
     * geneID	SRR5199304	SRR5199305	SRR5199306	SRR5199307	...
     * cajca.ICPL87119.gnm1.ann1.C.cajan_00002	0	0	0	0	...
     */
    void processExpression(File file) throws IOException {
	BufferedReader br = new BufferedReader(new FileReader(file));
	Item dataSet = getDataSet();
	Item organism = getOrganism();
	Item strain = getStrain(organism);
        List<Item> sampleList = new ArrayList<>();
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            if (parts[0].equals("geneID")) {
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
                        sample.setAttribute("primaryIdentifier", sampleId);
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
                    expressionValue.setReference("gene", gene);
                    expressionValue.setReference("sample", sample);
		    expressionValues.add(expressionValue);
                }
            }
        }
	br.close();
    }
}
