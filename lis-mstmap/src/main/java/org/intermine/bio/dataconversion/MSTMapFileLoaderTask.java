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

import java.util.List;
import java.util.LinkedList;
import java.util.Set;
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

import org.intermine.model.bio.Genotype;
import org.intermine.model.bio.GenotypingStudy;
import org.intermine.model.bio.GenotypingSample;
import org.intermine.model.bio.GenotypingRecord;

import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;

import org.ncgr.datastore.Readme;
import org.ncgr.zip.GZIPBufferedReader;

/**
 * Load genotyping study data from an MSTmap file (UC-Riverside).
 *
 * @author Sam Hokin
 */
public class MSTMapFileLoaderTask extends FileDirectDataLoaderTask {
    private static final Logger LOG = Logger.getLogger(MSTMapFileLoaderTask.class);

    // globals
    Organism organism;
    DataSource dataSource;
    DataSet dataSet;
    Publication publication;
    GenotypingStudy study;

    String dataSetUrl, dataSetVersion;

    /**
     * Process the files. super.process() calls processFile for each one.
     */
    @Override
    public void process() throws BuildException {
        DatastoreUtils datastoreUtils = new DatastoreUtils();
        try {
            // Organism - note no Strains loaded here
            organism = getDirectDataLoader().createObject(org.intermine.model.bio.Organism.class);
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
     * Process an MSTmap file or a README.
     *
     * @param file the File to process.
     * @throws BuildException if the is a problem
     */
    @Override
    public void processFile(File file) {
        if (file.getName().endsWith("mstmap.tsv.gz")) {
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
     * synopsis: MSTMap file containing genotype information for 31 wild and cultivated soybeans received from Meng Ni, Tin Hang, and Hon-Ming Lam.
     * description: MSTMap file from resequencing 31 wild and cultivated chinese soybean accessions ....
     * genbank_accession: SRA020131
     * publication_doi: 10.1038/ng.715
     */
    void processREADME(File file) throws IOException, ObjectStoreException {
        Readme readme = Readme.parse(file);
        // Organism
        organism.setTaxonId(String.valueOf(readme.taxid));
        // DataSet
        dataSet.setName(readme.identifier);
        dataSet.setDescription(readme.description);
        // Publication
        publication.setDoi(readme.publication_doi);
        publication.setTitle(readme.publication_title);
        // genotypes
        String genotypes = "";
        for (String genotype : readme.genotype) {
            if (genotypes.length()>0) genotypes += "|";
            genotypes += genotype;
        }
        study.setGenotypes(genotypes);
        // GenotypingStudy
        study.setPrimaryIdentifier(readme.identifier);
        study.setSynopsis(readme.synopsis);
        study.setDescription(readme.description);
        if (readme.genotyping_platform!=null) study.setGenotypingPlatform(readme.genotyping_platform);
        if (readme.genotyping_method!=null) study.setGenotypingMethod(readme.genotyping_method);
        if (readme.genbank_accession!=null) study.setGenbank(readme.genbank_accession);
        // references
        study.setOrganism(organism);
        study.setDataSet(dataSet);
        study.addPublications(publication);
    }

    /**
     * Process a gzipped MSTMap file. We assume here that all MSTMaps in a directory are for the same organism
     */
    public void processMSTMapFile(File file) throws ObjectStoreException, IOException, FileNotFoundException {
        List<GenotypingSample> samples = new LinkedList<>(); // has to support get(i)
        // spin through the MSTMap file
        String line = null;
        BufferedReader reader = GZIPBufferedReader.getReader(file);
        while ((line=reader.readLine())!=null) {
            if (line.length()==0 || line.startsWith("#")) continue;
            String[] fields = line.split("\t");
            if (samples.size()==0) {
                // first field is something like "locus_name"
                for (int i=1; i<fields.length; i++) {
                    String sampleName = fields[i];
                    GenotypingSample sample = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingSample.class);
                    sample.setPrimaryIdentifier(sampleName);
                    sample.setOrganism(organism);
                    sample.setStudy(study);
                    getDirectDataLoader().store(sample);
                    samples.add(sample);
                }
                study.setSamples(Set.copyOf(samples));
                LOG.info("Loaded "+samples.size()+" samples from MSTMap genotyping file.");
            } else {
                // marker name
                String markerName = fields[0];
                // GenotypingRecord
                GenotypingRecord record = getDirectDataLoader().createObject(org.intermine.model.bio.GenotypingRecord.class);
                record.setMarkerName(markerName);
                record.setStudy(study);
                record.setDataSet(dataSet);
                getDirectDataLoader().store(record);
                // Genotype
                for (int i=1; i<fields.length; i++) {
                    GenotypingSample sample = samples.get(i-1);
                    Genotype genotype = getDirectDataLoader().createSimpleObject(org.intermine.model.bio.Genotype.class);
                    genotype.setValue(fields[i]);
                    genotype.setSample(sample);
                    genotype.setRecord(record);
                    getDirectDataLoader().store(genotype);
                }
            }
        }
    }
}
