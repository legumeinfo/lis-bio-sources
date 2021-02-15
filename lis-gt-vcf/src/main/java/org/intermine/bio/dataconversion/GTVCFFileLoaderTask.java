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
import java.util.List;
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
import org.intermine.model.bio.Publication;
import org.intermine.model.bio.Strain;

import org.intermine.model.bio.GenotypingStudy;
import org.intermine.model.bio.GenotypingSample;
import org.intermine.model.bio.VCFRecord;
import org.intermine.model.bio.VCFSampleRecord;

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

    // used to get species, etc. from filename
    DatastoreUtils datastoreUtils;

    // up here since it's touched by both processREADME and processVCF
    GenotypingStudy study;

    // only want to store each chromosome and sample once
    Map<String,Chromosome> chromosomes = new HashMap<>();
    Map<String,GenotypingSample> samples = new HashMap<>();

    // set by setters
    String dataSourceName, dataSourceUrl, dataSourceDescription;
    String dataSetName, dataSetUrl, dataSetDescription, dataSetVersion, dataSetLicence;

    // assume singles of these per loader invocation
    Organism organism;
    Strain strain;
    DataSource dataSource;
    DataSet dataSet;

    /**
     * Process and load all of the fasta files.
     */
    @Override
    public void process() {
        try {
	    datastoreUtils = new DatastoreUtils();
            super.process();
            // store-once objects touched in multiple methods
            getDirectDataLoader().store(study);
            // wrap up
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
        System.err.println("********** processing "+file.getName()+" **********");
        if (file.getName().endsWith("vcf.gz")) {
            try {
                processVCFFile(file);
            } catch (Exception ex) {
                throw new BuildException(ex);
            }
        } else if (file.getName().startsWith("README")) {
            try {
                processREADME(file);
            } catch (Exception ex) {
                throw new BuildException(ex);
            }
        }
    }

    /**
     * DataSource.name for any bioentities created
     * @param dataSourceName name of datasource for items created
     */
    public void setDataSourceName(String dataSourceName) {
        this.dataSourceName = dataSourceName;
    }

    /**
     * DataSource.url for any bioentities created
     * @param dataSourceUrl url of datasource for items created
     */
    public void setDataSourceUrl(String dataSourceUrl) {
        this.dataSourceUrl = dataSourceUrl;
    }

    /**
     * DataSource.description for any bioentities created
     * @param dataSourceDescription description of datasource for items created
     */
    public void setDataSourceDescription(String dataSourceDescription) {
        this.dataSourceDescription = dataSourceDescription;
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
     * @param dataSetDescription the description of the DataSet of any new features
     */
    public void setDataSetDescription(String dataSetDescription) {
        this.dataSetDescription = dataSetDescription;
    }

    /**
     * If a value is specified this description will used when a DataSet is created.
     * @param dataSetVersion the version of the DataSet
     */
    public void setDataSetVersion(String dataSetVersion) {
        this.dataSetVersion = dataSetVersion;
    }

    /**
     * Process the README, which contains the GenotypingStudy attributes.
     *
     * README.0SZD.yml
     * ---------------
     * identifier: 0SZD
     * subject: VCF file containing genotype information for 31 wild and cultivated soybeans received from Meng Ni, Tin Hang, and Hon-Ming Lam.
     * description: VCF file from resequencing 31 wild and cultivated chinese soybean accessions ....
     * genbank_accession: SRA020131
     * publication_doi: 10.1038/ng.715
     */
    void processREADME(File file) throws IOException, ObjectStoreException {
        // GenotypingStudy:
        // <attribute name="identifier" type="java.lang.String"/>
        // <attribute name="subject" type="java.lang.String"/>
        // <attribute name="description" type="java.lang.String"/>
        // <attribute name="genbank" type=="java.lang.String"/>
        // <reference name="publication" referenced-type="Publication"/>
        Readme readme = Readme.getReadme(file);
        if (study==null) {
            study = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingStudy.class);
        }
        study.setPrimaryIdentifier(readme.identifier);
        study.setSubject(readme.subject);
        study.setDescription(readme.description);
        study.setGenbank(readme.genbank_accession);
        Publication publication = getDirectDataLoader().createObject(org.intermine.model.bio.Publication.class);
        publication.setDoi(readme.publication_doi);
        getDirectDataLoader().store(publication);
        study.setPublication(publication);
    }

    /**
     * Process a VCF file. We assume here that only one VCF/README is encountered per directory.
     *
     * 0     1    2    3   KEY4=identifier
     * glyma.Wm82.gnm2.div.0SZD.SNPData.vcf.gz
     */
    public void processVCFFile(File file) throws ObjectStoreException {
        // create and store the singleton objects
        Organism organism = getOrganismObject(file.getName());
        Strain strain = getStrainObject(file.getName(), organism);
        DataSource dataSource = getDataSourceObject();
        DataSet dataSet = getDataSetObject(dataSource);
        String assemblyVersion = DatastoreUtils.extractAssemblyVersion(file.getName());
        String annotationVersion = DatastoreUtils.extractAnnotationVersion(file.getName());
        // GenotypingStudy
        if (study==null) {
            study = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingStudy.class);
        }
        study.setOrganism(organism);
        study.setStrain(strain);
        study.setDataSet(dataSet);
        // spin through the VCF file
        // NOTE: index must be present!
        VCFFileReader vcfReader = new VCFFileReader(file);
        VCFHeader header = vcfReader.getFileHeader();
        // load samples first time
        if (samples.size()==0) {
            for (String sampleName : header.getGenotypeSamples()) {
                GenotypingSample sample = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingSample.class);
                sample.setPrimaryIdentifier(sampleName);
                sample.setOrganism(organism);
                sample.setStrain(strain);
                sample.setStudy(study);
                getDirectDataLoader().store(sample);
                samples.put(sampleName, sample);
            }
            study.setSamples(new HashSet(samples.values()));
            LOG.info("Loaded "+samples.size()+" samples from VCF header.");
        }
        // VCFRecords and VCFSampleRecords
        // #CHROM  POS ID REF ALT QUAL FILTER INFO                       FORMAT  C01         C02        C08          C12          C14         ...
        // Chr01   5   .  T   G   97   .      DP=327;VDB=1.18178e-13;... GT:PL   0/0:0,12,85 0/0:0,9,61 0/0:0,25,126 0/0:0,60,170 0/0:0,3,124 ...
        for (VariantContext vc : vcfReader.iterator().toList()) {
            String identifier = vc.getID();
            String contig = vc.getContig();
            Allele ref = vc.getReference();
            // ALT - if there is more than one, store the first
            List<Allele> alts = vc.getAlternateAlleles();
            Allele alt = alts.get(0);
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
            // cook an identifier like Gm03:544:A/C
            if (identifier==null || identifier.equals(".")) {
                identifier = contig+":"+vc.getStart()+":"+ref.getBaseString()+"/"+alt.getBaseString();
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
            // Location
            Location location = getDirectDataLoader().createObject(org.intermine.model.bio.Location.class);
            location.setStart(vc.getStart());
            location.setEnd(vc.getEnd());
            location.setLocatedOn(chromosome);
            location.addDataSets(dataSet);
            // VCFRecord
            VCFRecord vcfRecord = getDirectDataLoader().createObject(org.intermine.model.bio.VCFRecord.class);
            vcfRecord.setPrimaryIdentifier(identifier);
            vcfRecord.setOrganism(organism);
            vcfRecord.setStrain(strain);
            vcfRecord.setLength(vc.getEnd()-vc.getStart()+1);
            vcfRecord.setAssemblyVersion(assemblyVersion);
            vcfRecord.setAnnotationVersion(annotationVersion);
            vcfRecord.setRef(ref.getBaseString());
            vcfRecord.setAlt(alts.get(0).getBaseString());
            vcfRecord.setQual(qual);
            if (filters.length()>0) {
                vcfRecord.setFilter(filters);
            }
            vcfRecord.setInfo(info);
            vcfRecord.setGenotypingStudy(study);
            vcfRecord.setChromosome(chromosome);
            vcfRecord.setChromosomeLocation(location);
            vcfRecord.addDataSets(dataSet);
            // location.feature
            location.setFeature(vcfRecord);
            // le storage
            getDirectDataLoader().store(location);
            getDirectDataLoader().store(vcfRecord);
            // VCFSampleRecord
            for (String sampleName : samples.keySet()) {
                GenotypingSample sample = samples.get(sampleName);
                String genotype = vc.getGenotype(sampleName).getGenotypeString();
                String likelihoods = vc.getGenotype(sampleName).getLikelihoods().getAsString();
                VCFSampleRecord vcfSampleRecord = getDirectDataLoader().createSimpleObject(org.intermine.model.bio.VCFSampleRecord.class);
                vcfSampleRecord.setSample(sample);
                vcfSampleRecord.setGenotype(genotype);
                vcfSampleRecord.setLikelihoods(likelihoods);
                vcfSampleRecord.setVcfRecord(vcfRecord);
                getDirectDataLoader().store(vcfSampleRecord);
            }
        }
    }

    /**
     * Get/Create an Organism from a filename
     */
    Organism getOrganismObject(String filename) throws ObjectStoreException {
        if (organism!=null) {
            return organism;
        } else {
            organism = getDirectDataLoader().createObject(org.intermine.model.bio.Organism.class);
            String gensp = DatastoreUtils.extractGensp(filename);
            organism.setAbbreviation(gensp);
            organism.setTaxonId(datastoreUtils.getTaxonId(gensp));
            organism.setGenus(datastoreUtils.getGenus(gensp));
            organism.setSpecies(datastoreUtils.getSpecies(gensp));
            getDirectDataLoader().store(organism);
            return organism;
        }
    }

    /**
     * Get/Create a Strain from a filename and Organism
     */
    Strain getStrainObject(String filename, Organism organism) throws ObjectStoreException {
        if (strain!=null) {
            return strain;
        } else {
            strain = getDirectDataLoader().createObject(org.intermine.model.bio.Strain.class);
            strain.setIdentifier(DatastoreUtils.extractStrainIdentifier(filename));
            strain.setOrganism(organism);
            getDirectDataLoader().store(strain);
            return strain;
        }
    }

    /**
     * Get/Create a DataSource
     */
    DataSource getDataSourceObject() throws ObjectStoreException {
        if (dataSource!=null) {
            return dataSource;
        } else {
            dataSource = getDirectDataLoader().createObject(org.intermine.model.bio.DataSource.class);
            if (dataSourceName!=null) dataSource.setName(dataSourceName);
            if (dataSourceUrl!=null) dataSource.setUrl(dataSourceUrl);
            if (dataSourceDescription!=null) dataSource.setDescription(dataSourceDescription);
            getDirectDataLoader().store(dataSource);
            return dataSource;
        }
    }

    /**
     * Get/Create a DataSet
     */
    DataSet getDataSetObject(DataSource dataSource) throws ObjectStoreException {
        if (dataSet!=null) {
            return dataSet;
        } else {
            dataSet = getDirectDataLoader().createObject(org.intermine.model.bio.DataSet.class);
            if (dataSource!=null) dataSet.setDataSource(dataSource);
            if (dataSetName!=null) dataSet.setName(dataSetName);
            if (dataSetUrl!=null) dataSet.setUrl(dataSetUrl);
            if (dataSetDescription!=null) dataSet.setDescription(dataSetDescription);
            if (dataSetVersion!=null) dataSet.setVersion(dataSetVersion);
            if (dataSetLicence!=null) dataSet.setLicence(dataSetLicence);
            getDirectDataLoader().store(dataSet);
            return dataSet;
        }
    }
}
