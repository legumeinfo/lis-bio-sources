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

    // used to get species, etc. from filename
    DatastoreUtils datastoreUtils;

    // up here since it's touched by both processREADME and processVCF
    GenotypingStudy study;

    // only want to store these things once
    Map<String,Chromosome> chromosomes = new HashMap<>();
    Map<String,GenotypingSample> samples = new HashMap<>();

    // objects common to a load
    Organism organism;
    Strain strain;
    DataSource dataSource;
    DataSet dataSet;

    // attributes hopefully set by setters
    String dataSourceName, dataSourceUrl, dataSourceDescription;
    String dataSetName, dataSetUrl, dataSetVersion, dataSetLicence;

    /**
     * Process and load all of the included files.
     */
    @Override
    public void process() {
        try {
	    datastoreUtils = new DatastoreUtils();
            super.process();
            // store-once objects touched in multiple methods
            getDirectDataLoader().store(study);
            getDirectDataLoader().store(organism);
            getDirectDataLoader().store(strain);
            getDirectDataLoader().store(dataSource);
            getDirectDataLoader().store(dataSet);
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
    // public void close() throws ObjectStoreException {
    //     getDirectDataLoader().close();
    // }

    /**
     * @throws BuildException if an ObjectStore method fails
     */
    // @Override
    // public void execute() {
    //     // don't configure dynamic attributes if this is a unit test!
    //     if (getProject()!=null) {
    //         configureDynamicAttributes(this);
    //     }
    //     // this will call processFile() for each file
    //     super.execute();
    // }

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
     * synopsis: VCF file containing genotype information...
     * genbank_accession: SRA020131
     * publication_doi: 10.1038/ng.715
     */
    void processREADME(File file) throws IOException, ObjectStoreException {
        // GenotypingStudy:
        // <attribute name="identifier" type="java.lang.String"/>
        // <attribute name="synopsis" type="java.lang.String"/>
        // <attribute name="description" type="java.lang.String"/>
        // <attribute name="genbank" type=="java.lang.String"/>
        // <reference name="publication" referenced-type="Publication"/>
        Readme readme = Readme.getReadme(file);
        // check required stuff
        if (readme.identifier==null ||
            readme.synopsis==null ||
            readme.description==null ||
            readme.contributors==null ||
            readme.publication_doi==null ||
            readme.publication_title==null) {
            throw new BuildException("ERROR: a required field is missing from "+file.getName()+". "+
                                     "Required fields are: identifier, synopsis, description, contributors, publication_doi, publication_title");
        }
        // load required stuff
        if (study==null) {
            study = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingStudy.class);
        }
        study.setPrimaryIdentifier(readme.identifier);
        study.setSynopsis(readme.synopsis);
        study.setDescription(readme.description);
        study.setContributors(readme.contributors);
        Publication publication = getDirectDataLoader().createObject(org.intermine.model.bio.Publication.class);
        publication.setDoi(readme.publication_doi);
        publication.setTitle(readme.publication_title);
        getDirectDataLoader().store(publication);
        study.setPublication(publication);
        // optional fields
        if (readme.genbank_accession!=null) study.setGenbank(readme.genbank_accession);
        // dataSet.description
        dataSource = getDataSourceObject();
        dataSet = getDataSetObject(dataSource);
        dataSet.setDescription(readme.description);
    }

    /**
     * Process a VCF file. We assume here that all VCFs in a directory are for the same organism/strain.
     *
     * 0     1    2    3   KEY4=identifier
     * glyma.Wm82.gnm2.div.0SZD.SNPData.vcf.gz
     */
    public void processVCFFile(File file) throws ObjectStoreException {
        // singleton objects
        organism = getOrganismObject(file.getName());
        strain = getStrainObject(file.getName(), organism);
        dataSource = getDataSourceObject();
        dataSet = getDataSetObject(dataSource);
        // strings for full-yuckification of marker/ID
        String gensp = DatastoreUtils.extractGensp(file.getName());
        String strainName = DatastoreUtils.extractStrainIdentifier(file.getName());
        String assemblyVersion = DatastoreUtils.extractAssemblyVersion(file.getName());
        String annotationVersion = DatastoreUtils.extractAnnotationVersion(file.getName());
        String yuckyPrefix = gensp+"."+strainName+"."+assemblyVersion;
        if (annotationVersion!=null) yuckyPrefix += "."+annotationVersion;
        // GenotypingStudy
        if (study==null) {
            study = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingStudy.class);
        }
        study.setOrganism(organism);
        study.setStrain(strain);
        study.setDataSet(dataSet);
        // spin through the tabix-indexed VCF file
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
                samples.put(sampleName, sample);
                getDirectDataLoader().store(sample);
            }
            study.setSamples(new HashSet(samples.values()));
            LOG.info("Loaded "+samples.size()+" samples from VCF header.");
        }
        // VCFRecords and VCFSampleRecords
        // #CHROM  POS ID         REF ALT QUAL FILTER INFO                       FORMAT  C01         C02        C08          C12          C14         ...
        // Chr01   5   marker123  T   G   97   .      DP=327;VDB=1.18178e-13;... GT:PL   0/0:0,12,85 0/0:0,9,61 0/0:0,25,126 0/0:0,60,170 0/0:0,3,124 ...
        for (VariantContext vc : vcfReader.iterator().toList()) {
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
                // create a marker with a custom identifier from the location and ref/alts
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
                chromosomes.put(vc.getContig(), chromosome);
                getDirectDataLoader().store(chromosome);
            }
            // Location
            Location location = getDirectDataLoader().createObject(org.intermine.model.bio.Location.class);
            location.setStart(vc.getStart());
            location.setEnd(vc.getEnd());
            location.setLocatedOn(chromosome);
            location.addDataSets(dataSet);
            // GeneticMarker
            GeneticMarker marker = getDirectDataLoader().createObject(org.intermine.model.bio.GeneticMarker.class);
            location.setFeature(marker);
            marker.setPrimaryIdentifier(identifier);
            if (secondaryIdentifier!=null) marker.setSecondaryIdentifier(secondaryIdentifier);
            marker.setAssemblyVersion(assemblyVersion);
            marker.setAnnotationVersion(annotationVersion);
            marker.setType("VCF");
            marker.setAlleles(ref.getBaseString()+"/"+altString);
            marker.setOrganism(organism);
            marker.setStrain(strain);
            marker.setChromosome(chromosome);
            marker.setChromosomeLocation(location);
            marker.addDataSets(dataSet);
            // GenotypingRecord -- same identifier as marker
            GenotypingRecord genotypingRecord = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingRecord.class);
            genotypingRecord.setIdentifier(identifier);
            genotypingRecord.setMarker(marker);
            genotypingRecord.setRef(ref.getBaseString());
            genotypingRecord.setAlt(altString);
            genotypingRecord.setQual(qual);
            genotypingRecord.setInfo(info);
            if (filters.length()>0) genotypingRecord.setFilter(filters);
            genotypingRecord.setStudy(study);
            genotypingRecord.setDataSet(dataSet);
            // le storage
            getDirectDataLoader().store(marker);
            getDirectDataLoader().store(location);
            getDirectDataLoader().store(genotypingRecord);
            // Genotypes
            for (String sampleName : samples.keySet()) {
                GenotypingSample sample = samples.get(sampleName);
                String value = vc.getGenotype(sampleName).getGenotypeString();
                Genotype genotype = getDirectDataLoader().createSimpleObject(org.intermine.model.bio.Genotype.class);
                genotype.setSample(sample);
                genotype.setValue(value);
                if (vc.getGenotype(sampleName).hasLikelihoods()) {
                    genotype.setLikelihoods(vc.getGenotype(sampleName).getLikelihoods().getAsString());
                }
                genotype.setRecord(genotypingRecord);
                // le storage
                getDirectDataLoader().store(genotype);
            }
        }
    }

    /**
     * Get/Create an Organism from a filename
     */
    Organism getOrganismObject(String filename) throws ObjectStoreException {
        if (organism==null) {
            organism = getDirectDataLoader().createObject(org.intermine.model.bio.Organism.class);
            String gensp = DatastoreUtils.extractGensp(filename);
            organism.setAbbreviation(gensp);
            organism.setTaxonId(datastoreUtils.getTaxonId(gensp));
            organism.setGenus(datastoreUtils.getGenus(gensp));
            organism.setSpecies(datastoreUtils.getSpecies(gensp));
        }
        return organism;
    }

    /**
     * Get/Create a Strain from a filename and Organism
     */
    Strain getStrainObject(String filename, Organism organism) throws ObjectStoreException {
        if (strain==null) {
            strain = getDirectDataLoader().createObject(org.intermine.model.bio.Strain.class);
            String strainName = DatastoreUtils.extractStrainIdentifier(filename);
            strain.setIdentifier(strainName);
            strain.setOrganism(organism);
        }
        return strain;
    }

    /**
     * Get/Create a DataSource
     */
    DataSource getDataSourceObject() throws ObjectStoreException {
        if (dataSource==null) {
            dataSource = getDirectDataLoader().createObject(org.intermine.model.bio.DataSource.class);
            if (dataSourceName!=null) dataSource.setName(dataSourceName);
            if (dataSourceUrl!=null) dataSource.setUrl(dataSourceUrl);
            if (dataSourceDescription!=null) dataSource.setDescription(dataSourceDescription);
        }
        return dataSource;
    }

    /**
     * Get/Create a DataSet
     */
    DataSet getDataSetObject(DataSource dataSource) throws ObjectStoreException {
        if (dataSet==null) {
            dataSet = getDirectDataLoader().createObject(org.intermine.model.bio.DataSet.class);
            if (dataSource!=null) dataSet.setDataSource(dataSource);
            if (dataSetName!=null) dataSet.setName(dataSetName);
            if (dataSetUrl!=null) dataSet.setUrl(dataSetUrl);
            if (dataSetVersion!=null) dataSet.setVersion(dataSetVersion);
            if (dataSetLicence!=null) dataSet.setLicence(dataSetLicence);
        }
        return dataSet;
    }
}
