package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2021 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import java.util.LinkedList;
import java.util.HashSet;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.task.FileDirectDataLoaderTask;
import org.intermine.xml.full.Item;

import org.intermine.model.bio.DataSet;
import org.intermine.model.bio.DataSource;
import org.intermine.model.bio.Organism;
import org.intermine.model.bio.Publication;

import org.intermine.model.bio.GeneticMarker;
import org.intermine.model.bio.Genotype;
import org.intermine.model.bio.GenotypingStudy;
import org.intermine.model.bio.GenotypingSample;
import org.intermine.model.bio.GenotypingRecord;

import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;

/**
 * Load genotyping study data from an MSTmap file (UC-Riverside).
 *
 * @author Sam Hokin
 */
public class MSTMapFileLoaderTask extends FileDirectDataLoaderTask {
    private static final Logger LOG = Logger.getLogger(MSTMapFileLoaderTask.class);

    // used to get species, etc. from filename
    DatastoreUtils datastoreUtils;

    // ordered list of samples
    LinkedList<GenotypingSample> samples = new LinkedList<>();

    // up here since it's touched by both processREADME and processVCF
    GenotypingStudy study;

    // other objects common to any loader
    Organism organism;
    DataSource dataSource;
    DataSet dataSet;

    // attributes hopefully set by setters
    String dataSourceName, dataSourceUrl, dataSourceDescription;
    String dataSetName, dataSetUrl, dataSetVersion, dataSetLicence;

    /**
     * Process and load all of the fasta files.
     */
    @Override
    public void process() {
        datastoreUtils = new DatastoreUtils();
        super.process();
        try {
            // store global objects
            getDirectDataLoader().store(study);
            getDirectDataLoader().store(organism);
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
        if (file.getName().endsWith("MSTmap.txt")) {
            try {
                System.err.println("Processing "+file.getName());
                processMSTMapFile(file);
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
     * If a value is specified this version will used when a DataSet is created.
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
     * subject: MSTMap file containing genotype information for 31 wild and cultivated soybeans received from Meng Ni, Tin Hang, and Hon-Ming Lam.
     * description: MSTMap file from resequencing 31 wild and cultivated chinese soybean accessions ....
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
        study.setContributors(readme.contributors);
        Publication publication = getDirectDataLoader().createObject(org.intermine.model.bio.Publication.class);
        publication.setDoi(readme.publication_doi);
        publication.setTitle(readme.publication_title);
        getDirectDataLoader().store(publication);
        study.setPublication(publication);
        // set DataSet.description from README
        dataSource = getDataSourceObject();
        dataSet = getDataSetObject(dataSource);
        dataSet.setDescription(readme.description);
    }

    /**
     * Process a MSTMap file. We assume here that all MSTMaps in a directory are for the same organism
     *
     * 0     1    2    3   KEY4=identifier
     * glyma.Wm82.gnm2.div.0SZD.SNPData.vcf.gz
     */
    public void processMSTMapFile(File file) throws ObjectStoreException, IOException, FileNotFoundException {
        // create and store the singleton objects
        organism = getOrganismObject(file.getName());
        dataSource = getDataSourceObject();
        dataSet = getDataSetObject(dataSource);
        // GenotypingStudy
        if (study==null) {
            study = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingStudy.class);
        }
        study.setOrganism(organism);
        study.setDataSet(dataSet);
        study.setParents(extractParents(file.getName()));
        // spin through the MSTMap file
        boolean haveReadSamples = false;
        String line = null;
        BufferedReader reader = new BufferedReader(new FileReader(file));
        while ((line=reader.readLine())!=null) {
            if (line.length()==0 || line.startsWith("#")) continue;
            String[] fields = line.split("\t");
            if (!haveReadSamples) {
                // first field is something like "locus_name"
                for (int i=1; i<fields.length; i++) {
                    String sampleName = fields[i];
                    GenotypingSample sample = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingSample.class);
                    sample.setPrimaryIdentifier(sampleName);
                    sample.setOrganism(organism);
                    sample.setStudy(study);
                    samples.add(sample);
                    getDirectDataLoader().store(sample);
                }
                study.setSamples(new HashSet(samples));
                haveReadSamples = true;
                LOG.info("Loaded "+samples.size()+" samples from MSTMap genotyping file.");
            } else {
                // GeneticMarker
                String markerSecondaryIdentifier = fields[0];
                GeneticMarker marker = getDirectDataLoader().createObject(org.intermine.model.bio.GeneticMarker.class);
                marker.setSecondaryIdentifier(markerSecondaryIdentifier);
                marker.setOrganism(organism);
                marker.addDataSets(dataSet);
                getDirectDataLoader().store(marker);
                // GenotypingRecord
                GenotypingRecord genotypingRecord = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingRecord.class);
                genotypingRecord.setIdentifier(markerSecondaryIdentifier);
                genotypingRecord.setStudy(study);
                genotypingRecord.setMarker(marker);
                getDirectDataLoader().store(genotypingRecord);
                // Genotype
                for (int i=1; i<fields.length; i++) {
                    GenotypingSample sample = samples.get(i-1);
                    Genotype genotype = getDirectDataLoader().createSimpleObject(org.intermine.model.bio.Genotype.class);
                    genotype.setValue(fields[i]);
                    genotype.setSample(sample);
                    genotype.setRecord(genotypingRecord);
                    try {
                        getDirectDataLoader().store(genotype);
                    } catch (Exception ex) {
                        // DEBUG
                        System.err.println("Error storing "+genotype.toString());
                        System.err.println(ex);
                    }
                }
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

    /**
     * Extract the parents from the given filename
     * ex. vigun.CB27_x_IT82E-18.gt.KEY4.MSTmap.txt
     */
    static String extractParents(String filename) {
        String[] fields = filename.split("\\.");
        return fields[1];
    }
}
