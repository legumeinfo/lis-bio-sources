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

/**
 * Load genotyping study data from a VCF files.
 *
 * @author Sam Hokin
 */
public class GTVCFFileLoaderTask extends FileDirectDataLoaderTask {
    private static final Logger LOG = Logger.getLogger(GTVCFFileLoaderTask.class);

    // singletons
    Organism organism;
    Strain strain;
    DataSource dataSource;
    DataSet dataSet;
    Publication publication;
    GenotypingStudy study;

    String taxonId;
    String dataSetUrl, dataSetVersion;

    /**
     * Process and load all of the included files.
     */
    @Override
    public void process() {
        DatastoreUtils datastoreUtils = new DatastoreUtils();
        try {
            // Organism
            organism = getDirectDataLoader().createObject(org.intermine.model.bio.Organism.class);
            organism.setTaxonId(taxonId);
            // Strain
            strain = getDirectDataLoader().createObject(org.intermine.model.bio.Strain.class);
            // DataSource
            dataSource = getDirectDataLoader().createObject(org.intermine.model.bio.DataSource.class);
            dataSource.setName(DatastoreFileConverter.DEFAULT_DATASOURCE_NAME);
            dataSource.setUrl(DatastoreFileConverter.DEFAULT_DATASOURCE_URL);
            dataSource.setDescription(DatastoreFileConverter.DEFAULT_DATASOURCE_DESCRIPTION);
            // DataSet - most set from README
            dataSet = getDirectDataLoader().createObject(org.intermine.model.bio.DataSet.class);
            dataSet.setDataSource(dataSource);
            dataSet.setLicence(DatastoreFileConverter.DEFAULT_DATASET_LICENCE);
            if (dataSetUrl!=null) dataSet.setUrl(dataSetUrl);
            if (dataSetVersion!=null) dataSet.setVersion(dataSetVersion);
            // Publication - set from README
            publication = getDirectDataLoader().createObject(org.intermine.model.bio.Publication.class);
            // Genotyping Study
            study = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingStudy.class);
        } catch (ObjectStoreException ex) {
            throw new BuildException(ex);
        }
        // process the files
        super.process();
        // wrap up
        try {
            // store global objects
            getDirectDataLoader().store(study);
            getDirectDataLoader().store(organism);
            getDirectDataLoader().store(dataSource);
            getDirectDataLoader().store(dataSet);
            getDirectDataLoader().store(publication);
            getIntegrationWriter().commitTransaction();
            getIntegrationWriter().beginTransaction();
        } catch (ObjectStoreException e) {
            throw new BuildException("Failed to store object", e);
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
     * Process a VCF file or a README.
     *
     * @param file the File to process.
     * @throws BuildException if the is a problem
     */
    @Override
    public void processFile(File file) {
        if (file.getName().endsWith("vcf.gz")) {
            try {
                System.err.println("Processing "+file.getName());
                processVCFFile(file);
            } catch (Exception ex) {
                throw new BuildException(ex);
            }
        } else if (file.getName().startsWith("README")) {
            try {
                System.err.println("Processing "+file.getName());
                processREADME(file);
            } catch (Exception ex) {
                throw new BuildException(ex);
            }
        } else {
            System.err.println("Skipping "+file.getName());
        }
    }

    /**
     * If a value is specified this url will used when a DataSet is created.
     * @param dataSetUrl the title of the DataSet of any new features
     */
    public void setDataSetUrl(String dataSetUrl) {
        this.dataSetUrl = dataSetUrl;
    }

    /**
     * If a value is specified this description will used when a DataSet is created.
     * @param dataSetVersion the version of the DataSet
     */
    public void setDataSetVersion(String dataSetVersion) {
        this.dataSetVersion = dataSetVersion;
    }

    /**
     * Set the Taxon ID to be used for the organism.
     * @param taxonId the NCBI taxonomy ID
     */
    public void setTaxonId(String value) {
        this.taxonId = value;
    }

    /**
     * Process the README, which contains the GenotypingStudy attributes.
     */
    void processREADME(File file) throws IOException, ObjectStoreException {
        Readme readme = Readme.getReadme(file);
        // Organism -- override taxonId that may have been set in project.xml
        organism.setTaxonId(readme.taxid);
        // DataSet
        dataSet.setName(readme.identifier);
        dataSet.setDescription(readme.description);
        // Publication
        publication.setDoi(readme.publication_doi);
        publication.setTitle(readme.publication_title);
        // Populations (from README.genotype)
        for (String genotype : readme.genotype) {
            Population population = getDirectDataLoader().createObject(org.intermine.model.bio.Population.class);
            population.setIdentifier(genotype);
            getDirectDataLoader().store(population);
            study.addPopulations(population);
        }
        // GenotypingStudy
        study.setPrimaryIdentifier(readme.identifier);
        study.setSynopsis(readme.synopsis);
        study.setDescription(readme.description);
        if (readme.genotyping_platform!=null) study.setGenotypingPlatform(readme.genotyping_platform);
        if (readme.genotyping_method!=null) study.setGenotypingMethod(readme.genotyping_method);
        if (readme.genbank_accession!=null) study.setGenbank(readme.genbank_accession);
        study.setOrganism(organism);
        study.setDataSet(dataSet);
        study.addPublications(publication);
    }

    /**
     * Process a VCF file.
     */
    public void processVCFFile(File file) throws ObjectStoreException {
        Map<String,Chromosome> chromosomes = new HashMap<>();
        List<GenotypingSample> samples = new LinkedList<>(); // has to support get(i)
        // Strain is for chromosomes and markers only
        strain.setIdentifier(DatastoreUtils.extractStrainIdentifier(file.getName()));
        strain.setOrganism(organism);
        getDirectDataLoader().store(strain);
        // form full-yuck prefix for marker 
        String gensp = DatastoreUtils.extractGensp(file.getName());
        String strainName = DatastoreUtils.extractStrainIdentifier(file.getName());
        String assemblyVersion = DatastoreUtils.extractAssemblyVersion(file.getName());
        String annotationVersion = DatastoreUtils.extractAnnotationVersion(file.getName());
        String yuckyPrefix = gensp+"."+strainName+"."+assemblyVersion;
        if (annotationVersion!=null) yuckyPrefix += "."+annotationVersion;
        // spin through the tabix-indexed VCF file
        VCFFileReader vcfReader = new VCFFileReader(file);
        VCFHeader header = vcfReader.getFileHeader();
        // load samples first time
        if (samples.size()==0) {
            for (String sampleName : header.getGenotypeSamples()) {
                GenotypingSample sample = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingSample.class);
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
                chromosome = getDirectDataLoader().createObject(org.intermine.model.bio.Chromosome.class);
                chromosome.setPrimaryIdentifier(vc.getContig());
                chromosome.setOrganism(organism);
                chromosome.setStrain(strain);
                chromosome.addDataSets(dataSet);
                getDirectDataLoader().store(chromosome);
                chromosomes.put(vc.getContig(), chromosome);
            }
            // GeneticMarker
            GeneticMarker marker = getDirectDataLoader().createObject(org.intermine.model.bio.GeneticMarker.class);
            marker.setPrimaryIdentifier(identifier);
            if (secondaryIdentifier!=null) marker.setSecondaryIdentifier(secondaryIdentifier);
            marker.setAssemblyVersion(assemblyVersion);
            marker.setAnnotationVersion(annotationVersion);
            marker.setAlleles(ref.getBaseString()+"/"+altString);
            marker.setOrganism(organism);
            marker.setStrain(strain);
            marker.setChromosome(chromosome);
            marker.addDataSets(dataSet);
            // Location
            Location location = getDirectDataLoader().createObject(org.intermine.model.bio.Location.class);
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
            GenotypingRecord record = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingRecord.class);
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
                Genotype genotype = getDirectDataLoader().createSimpleObject(org.intermine.model.bio.Genotype.class);
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
