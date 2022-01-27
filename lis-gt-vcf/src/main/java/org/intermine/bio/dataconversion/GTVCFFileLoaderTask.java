package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2021 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 */

import java.io.File;
import java.io.IOException;

import java.util.Map;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.LinkedList;
import java.util.Set;
import java.util.HashSet;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.task.FileDirectDataLoaderTask;
import org.intermine.xml.full.Item;

import org.intermine.model.bio.Chromosome;
import org.intermine.model.bio.DataSet;
import org.intermine.model.bio.DataSource;
import org.intermine.model.bio.Location;
import org.intermine.model.bio.Organism;
import org.intermine.model.bio.Population;
import org.intermine.model.bio.Publication;
import org.intermine.model.bio.Strain;

import org.intermine.model.bio.GeneticMarker;
import org.intermine.model.bio.GenotypingStudy;
import org.intermine.model.bio.GenotypingSample;
import org.intermine.model.bio.GenotypingRecord;
import org.intermine.model.bio.Genotype;

import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import org.ncgr.datastore.Readme;

/**
 * Load genotyping study data from a VCF files.
 *
 * @author Sam Hokin
 */
public class GTVCFFileLoaderTask extends FileDirectDataLoaderTask {
    private static final Logger LOG = Logger.getLogger(GTVCFFileLoaderTask.class);

    // project.xml setters
    String dataSetUrl, dataSetLicence;

    // collection items
    Organism organism;
    Strain strain;
    DataSource dataSource;
    DataSet dataSet;
    Publication publication;

    // local items
    GenotypingStudy study;
    String assemblyVersion, annotationVersion; // attributes of markers
    String yuckyPrefix;                        // for prefixing marker identifiers

    /**
     * If a value is specified this url will used when a DataSet is created.
     * @param dataSetUrl the title of the DataSet of any new features
     */
    public void setDataSetUrl(String dataSetUrl) {
        this.dataSetUrl = dataSetUrl;
    }

    /**
     * dataSetLicence is set in project.xml
     */
    public void setDataSetLicence(String licence) {
	this.dataSetLicence = licence;
    }
    
    /**
     * Process and load all of the included files.
     */
    @Override
    public void process() {
        // process files
        super.process();
        try {
            // store collection items
            getDirectDataLoader().store(organism);
            getDirectDataLoader().store(strain);
            getDirectDataLoader().store(dataSource);
            getDirectDataLoader().store(dataSet);
            if (publication!=null) getDirectDataLoader().store(publication);
            // store local items
            getDirectDataLoader().store(study);
            // DO WE NEED THESE?
            // getIntegrationWriter().commitTransaction();
            // getIntegrationWriter().beginTransaction();
        } catch (ObjectStoreException ex) {
            throw new BuildException("Failed to store object", ex);
        }
    }

    /**
     * Be sure to close the data loader so the last batch gets stored. only needed for tests
     * since the data loading task usually does that for the live builds.
     * @throws ObjectStoreException if we can't store to db
     */
    public void close() throws ObjectStoreException {
        getDirectDataLoader().close();
    }

    /**
     * @throws BuildException if an ObjectStore method fails
     */
    @Override
    public void execute() {
        // don't configure dynamic attributes if this is a unit test!
        if (getProject()!=null) {
            configureDynamicAttributes(this);
        }
        // this will call processFile() for each file
        super.execute();
    }

    /**
     * Process a README and then a VCF file.
     *
     * @param file the File to process.
     * @throws BuildException if the is a problem
     */
    @Override
    public void processFile(File file) {
        if (file.getName().startsWith("README")) {
            try {
                processReadme(file);
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            } catch (ObjectStoreException ex) {
                throw new RuntimeException(ex);
            }
        } else if (file.getName().endsWith("vcf.gz")) {
            System.out.println("##########");
            System.out.println("# Processing "+file.getName());
            System.out.println("##########");
            try {
                processVCFFile(file);
            } catch (ObjectStoreException ex) {
                throw new RuntimeException(ex);
            }
        }
    }

    /**
     * Process the README file for all the collection data.
     * Note: we get strain from the collection identifier, not the free-form genotype entries.
     */
    void processReadme(File file) throws IOException, ObjectStoreException {
        Readme readme = Readme.parse(file);
        // project.xml check
        if (dataSetUrl==null) {
            throw new RuntimeException("ERROR: dataSetUrl must be set in project.xml.");
        }
        // check required stuff
        if (readme.identifier==null ||
            readme.taxid==0 ||
            readme.synopsis==null ||
            readme.description==null ||
            readme.scientific_name_abbrev==null ||
            readme.genotype==null
            ) {
            throw new RuntimeException("ERROR in README: a required field is missing. "+
                                       "Required fields are: identifier, taxid, synopsis, description, scientific_name_abbrev, genotype");
        }
        String collection = DatastoreUtils.extractCollectionFromReadme(file);
        // DataSet
        dataSet = getDirectDataLoader().createObject(DataSet.class);
        dataSet.setName(readme.identifier);
        dataSet.setSynopsis(readme.synopsis);
        dataSet.setDescription(readme.description);
        assemblyVersion = DatastoreUtils.extractAssemblyVersionFromCollection(readme.identifier);
        annotationVersion = DatastoreUtils.extractAnnotationVersionFromCollection(readme.identifier);
        if (assemblyVersion!=null && annotationVersion!=null) {
            dataSet.setVersion(assemblyVersion+"."+annotationVersion);
        } else if (assemblyVersion!=null) {
            dataSet.setVersion(assemblyVersion);
        }
        if (dataSetLicence!=null) {
            dataSet.setLicence(dataSetLicence);
        } else {
            dataSet.setLicence(DatastoreFileConverter.DEFAULT_DATASET_LICENCE);
        }
        dataSet.setUrl(dataSetUrl);
        // Publication
        if (readme.publication_doi!=null) {
            publication = getDirectDataLoader().createObject(Publication.class);
            publication.setDoi(readme.publication_doi);
            if (readme.publication_title!=null) publication.setTitle(readme.publication_title);
            dataSet.setPublication(publication);
        }
        // Organism
        organism = getDirectDataLoader().createObject(Organism.class);
        organism.setTaxonId(String.valueOf(readme.taxid));
        // Strain
        String strainIdentifier = DatastoreUtils.extractStrainIdentifierFromCollection(collection);
        if (strainIdentifier==null) {
            throw new RuntimeException("ERROR: could not extract strain identifier from "+collection+".");
        }
        strain = getDirectDataLoader().createObject(Strain.class);
        strain.setIdentifier(strainIdentifier);
        strain.setOrganism(organism);
        // GenotypingStudy
        study = getDirectDataLoader().createObject(GenotypingStudy.class);
        study.setPrimaryIdentifier(readme.identifier);
        study.setSynopsis(readme.synopsis);
        study.setOrganism(organism);
        study.setDataSet(dataSet);
        study.setDescription(readme.description);
        if (publication!=null) study.addPublications(publication);
        if (readme.genotyping_platform!=null) study.setGenotypingPlatform(readme.genotyping_platform);
        if (readme.genotyping_method!=null) study.setGenotypingMethod(readme.genotyping_method);
        if (readme.genbank_accession!=null) study.setGenbank(readme.genbank_accession);
        // Genotyping populations (from README.genotype)
        for (String genotype : readme.genotype) {
            Population population = getDirectDataLoader().createObject(Population.class);
            population.setIdentifier(genotype);
            getDirectDataLoader().store(population);
            study.addPopulations(population);
        }
        // form full-yuck prefix for markers
        String yuckyPrefix = readme.scientific_name_abbrev+"."+strainIdentifier+"."+assemblyVersion;
        if (annotationVersion!=null) yuckyPrefix += "."+annotationVersion;
    }

    /**
     * Process a VCF file.
     * NOTE: assumption is that variants are positioned on chromosomes, not supercontigs.
     */
    public void processVCFFile(File file) throws ObjectStoreException {
        // check that README has already been processed
        if (study==null) {
            throw new RuntimeException("ERROR: README must be before VCF file in project.xml.");
        }
        Map<String,Chromosome> chromosomes = new HashMap<>();
        List<GenotypingSample> samples = new LinkedList<>(); // has to support get(i)
        // spin through the tabix-indexed VCF file
        VCFFileReader vcfReader = new VCFFileReader(file);
        VCFHeader header = vcfReader.getFileHeader();
        // load samples first time
        if (samples.size()==0) {
            for (String sampleName : header.getGenotypeSamples()) {
                GenotypingSample sample = getDirectDataLoader().createObject(GenotypingSample.class);
                sample.setPrimaryIdentifier(sampleName);
                sample.setOrganism(organism);
                sample.setStudy(study);
                getDirectDataLoader().store(sample);
                samples.add(sample);
            }
            study.setSamples(Set.copyOf(samples));
            LOG.info("Loaded "+samples.size()+" samples from VCF header.");
        }
        // VCFRecords and VCFSampleRecords
        // #CHROM  POS ID         REF ALT QUAL FILTER INFO                       FORMAT  C01         C02        C08          C12          C14         ...
        // Chr01   5   marker123  T   G   97   .      DP=327;VDB=1.18178e-13;... GT:PL   0/0:0,12,85 0/0:0,9,61 0/0:0,25,126 0/0:0,60,170 0/0:0,3,124 ...
        Iterator<VariantContext> it = vcfReader.iterator();
        while (it.hasNext()) {
            VariantContext vc = it.next();
            String contig = vc.getContig();
            String identifier = vc.getID();
            Allele ref = vc.getReference();
            // ALTs
            List<Allele> alts = vc.getAlternateAlleles();
            String altString = alts.get(0).getBaseString();
            for (int i=1; i<alts.size(); i++) {
                altString += ","+alts.get(i).getBaseString();
            }
            // QUAL
            double qual = vc.getPhredScaledQual();
            // FILTER
            String filters = "";
            Set<String> filterSet = vc.getFilters();
            for (String filter : filterSet) {
                filters += filter;
            }
            // INFO data
            String info = "";
            Map<String,Object> attributes = vc.getCommonInfo().getAttributes();
            for (String key : attributes.keySet()) {
                info += key+"="+attributes.get(key).toString()+";";
            }
            // convert to full-yuck identifier if present and without dots
            boolean hasMarker = !identifier.equals(".");
            String secondaryIdentifier = null;
            if (hasMarker && !identifier.contains("\\.")) {
                // we have a non-full-yuck marker
                secondaryIdentifier = identifier;
                identifier = yuckyPrefix+"."+identifier;
            } else if (!hasMarker) {
                // create a marker with a custom identifier from the locationand ref/alts
                identifier = contig+":"+vc.getStart()+":"+ref.getBaseString()+"/"+altString;
            }
            // Chromosome
            Chromosome chromosome = chromosomes.get(vc.getContig());
            if (chromosome==null) {
                chromosome = getDirectDataLoader().createObject(Chromosome.class);
                chromosome.setPrimaryIdentifier(vc.getContig());
                chromosome.setOrganism(organism);
                chromosome.setStrain(strain);
                chromosome.addDataSets(dataSet);
                getDirectDataLoader().store(chromosome);
                chromosomes.put(vc.getContig(), chromosome);
            }
            // GeneticMarker
            GeneticMarker marker = getDirectDataLoader().createObject(GeneticMarker.class);
            marker.setPrimaryIdentifier(identifier);
            if (secondaryIdentifier!=null) marker.setSecondaryIdentifier(secondaryIdentifier);
            if (assemblyVersion!=null) marker.setAssemblyVersion(assemblyVersion);
            if (annotationVersion!=null) marker.setAnnotationVersion(annotationVersion);
            marker.setAlleles(ref.getBaseString()+"/"+altString);
            marker.setOrganism(organism);
            marker.setStrain(strain);
            marker.setChromosome(chromosome);
            marker.addDataSets(dataSet);
            // Location
            Location location = getDirectDataLoader().createObject(Location.class);
            location.setStart(vc.getStart());
            location.setEnd(vc.getEnd());
            location.setLocatedOn(chromosome);
            location.addDataSets(dataSet);
            // have to store location and marker together since mutually referencing
            location.setFeature(marker);
            marker.setChromosomeLocation(location);
            getDirectDataLoader().store(location);
            getDirectDataLoader().store(marker);
            // GenotypingRecord -- same identifier as marker
            GenotypingRecord record = getDirectDataLoader().createObject(GenotypingRecord.class);
            record.setIdentifier(identifier);
            record.setMarker(marker);
            record.setRef(ref.getBaseString());
            record.setAlt(altString);
            record.setQual(qual);
            record.setInfo(info);
            if (filters.length()>0) record.setFilter(filters);
            record.setStudy(study);
            record.setDataSet(dataSet);
            getDirectDataLoader().store(record);
            // Genotypes
            for (GenotypingSample sample : samples) {
                String sampleName = sample.getPrimaryIdentifier();
                String value = vc.getGenotype(sampleName).getGenotypeString();
                Genotype genotype = getDirectDataLoader().createSimpleObject(Genotype.class);
                genotype.setSample(sample);
                genotype.setValue(value);
                if (vc.getGenotype(sampleName).hasLikelihoods()) {
                    genotype.setLikelihoods(vc.getGenotype(sampleName).getLikelihoods().getAsString());
                }
                genotype.setRecord(record);
                getDirectDataLoader().store(genotype);
            }
        }
    }
}
