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
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStream;
import java.io.IOException;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.NoSuchElementException;

import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;

import org.biojava.nbio.core.exceptions.ParserException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import org.intermine.metadata.Util;
import org.intermine.model.InterMineObject;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.objectstore.query.PendingClob;
import org.intermine.task.FileDirectDataLoaderTask;

import org.intermine.model.bio.BioEntity;
import org.intermine.model.bio.Chromosome;
import org.intermine.model.bio.DataSet;
import org.intermine.model.bio.DataSource;
import org.intermine.model.bio.Organism;
import org.intermine.model.bio.Publication;
import org.intermine.model.bio.Strain;
import org.intermine.model.bio.Sequence;
import org.intermine.model.bio.SequenceFeature;
import org.intermine.model.bio.Supercontig;

import org.ncgr.datastore.Readme;
import org.ncgr.datastore.validation.GenomeCollectionValidator;
import org.ncgr.zip.GZIPFastaReader;

/**
 * A task that can read a set of FASTA files and create the corresponding Sequence objects in an ObjectStore.
 *
 * @author Kim Rutherford
 * @author Peter Mclaren
 * @author Sam Hokin
 */

public class GenomeFastaLoaderTask extends FileDirectDataLoaderTask {
    static final Logger LOG = Logger.getLogger(GenomeFastaLoaderTask.class);
    static final String SEQUENCE_TYPE = "dna";

    boolean loadHeaderDescriptions = false;

    // project.xml setters
    String idAttribute, descriptionAttribute;
    String dataSourceName, dataSourceUrl, dataSourceDescription;
    String dataSetUrl, dataSetLicence;

    Readme readme;
    String gensp, strainIdentifier, assemblyVersion, annotationVersion;

    List<String> chromosomePrefixes = new ArrayList<>();
    List<String> supercontigPrefixes = new ArrayList<>();
    
    // collection items
    Organism organism;
    Strain strain;
    DataSource dataSource;
    DataSet dataSet;
    Publication publication;
    
    boolean fastaProcessed = false; // flag to indicate that we processed the FASTA
    boolean collectionValidated = false; // validate the collection first by storing a flag

    /**
     * dataSetUrl can be set in project.xml
     */
    public void setDataSetUrl(String url) {
        this.dataSetUrl = url;
    }
    
    /**
     * dataSetLicence can be set in project.xml
     */
    public void setDataSetLicence(String licence) {
	this.dataSetLicence = licence;
    }

    /**
     * dataSourceName can be set in project.xml
     */
    public void setDataSourceName(String name) {
	this.dataSourceName = name;
    }
    /**
     * dataSourceUrl can be set in project.xml
     */
    public void setDataSourceUrl(String url) {
	this.dataSourceUrl = url;
    }
    /**
     * dataSourceDescription can be set in project.xml
     */
    public void setDataSourceDescription(String description) {
	this.dataSourceDescription = description;
    }

    /**
     * Process and load all of the fasta files.
     */
    @Override
    public void process() {
        // process files, which stores the features and sequences directly
        super.process();
        // bail if we haven't processed the FASTA
        if (!fastaProcessed) {
            throw new BuildException("FASTA file not processed - aborting.");
        }
        // store the extra objects with direct data loader
        try {
            // store collection items
            getDirectDataLoader().store(organism);
            getDirectDataLoader().store(strain);
            getDirectDataLoader().store(dataSource);
            getDirectDataLoader().store(dataSet);
            if (publication!=null) getDirectDataLoader().store(publication);
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
     * Handles each file.
     *
     * @param file the File to process.
     * @throws BuildException if there is a problem
     */
    @Override
    public void processFile(File file) throws BuildException {
        if (!collectionValidated) {
            GenomeCollectionValidator validator = new GenomeCollectionValidator(file.getParent());
            validator.validate();
            if (!validator.isValid()) {
                throw new BuildException("Collection "+file.getParent()+" does not pass validation.");
            }
            collectionValidated = true;
        }
        if (file.getName().startsWith("README")) {
            try {
                processReadme(file);
            } catch (IOException ex) {
                throw new BuildException(ex);
            } catch (ObjectStoreException ex) {
                throw new BuildException(ex);
            }
        } else if (file.getName().endsWith(".fna.gz")) {
            // README must precede FASTA since we have direct storage
            if (organism==null) {
                throw new BuildException("README missing or not read before FASTA file. Add to includes or switch order in project.xml.");
            }
            System.out.println("## Reading "+SEQUENCE_TYPE+" sequences from: "+file.getName());
            processFasta(file);
        } else {
            System.out.println("## - Skipping "+file.getName());
	}
    }

    /**
     * Process the README file for all the collection data.
     * Note: we get strain from the collection identifier, not the free-form genotype entries.
     */
    void processReadme(File file) throws IOException, ObjectStoreException {
        readme = Readme.parse(file);
        // project.xml check
        if (dataSetUrl==null) {
            throw new BuildException("dataSetUrl must be set in project.xml.");
        }
        // chromosome, supercontig prefixes
        if (readme.chromosome_prefix!=null) {
            for (String prefix : readme.chromosome_prefix.split(",")) {
                chromosomePrefixes.add(prefix);
            }
        }
        if (readme.supercontig_prefix!=null) {
            for (String prefix : readme.supercontig_prefix.split(",")) {
                supercontigPrefixes.add(prefix);
            }
        }
        // local vars
        gensp = readme.scientific_name_abbrev;
        strainIdentifier = DatastoreUtils.extractStrainIdentifierFromCollection(readme.identifier);
        assemblyVersion = DatastoreUtils.extractAssemblyVersionFromCollection(readme.identifier);
        annotationVersion = DatastoreUtils.extractAnnotationVersionFromCollection(readme.identifier);
        // DataSource
        dataSource = getDirectDataLoader().createObject(DataSource.class);
        if (dataSourceName==null) {
            dataSource.setName(DatastoreFileConverter.DEFAULT_DATASOURCE_NAME);
        } else {
            dataSource.setName(dataSourceName);
        }
        if (dataSourceUrl==null) {
            dataSource.setUrl(DatastoreFileConverter.DEFAULT_DATASOURCE_URL);
        } else {
            dataSource.setUrl(dataSourceUrl);
        }
        if (dataSourceDescription==null) {
            dataSource.setDescription(DatastoreFileConverter.DEFAULT_DATASOURCE_DESCRIPTION);
        } else {
            dataSource.setDescription(dataSourceDescription);
        }
        // DataSet
        dataSet = getDirectDataLoader().createObject(DataSet.class);
        dataSet.setDataSource(dataSource);
        dataSet.setName(readme.identifier);
        dataSet.setSynopsis(readme.synopsis);
        dataSet.setDescription(readme.description);
        dataSet.setUrl(dataSetUrl); // required in project.xml
        if (assemblyVersion!=null && annotationVersion!=null) {
            dataSet.setVersion(assemblyVersion+"."+annotationVersion);
        } else if (assemblyVersion!=null) {
            dataSet.setVersion(assemblyVersion);
        }
        if (dataSetLicence==null) {
            dataSet.setLicence(DatastoreFileConverter.DEFAULT_DATASET_LICENCE);
        } else {
            dataSet.setLicence(dataSetLicence);
        }
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
        organism.addDataSets(dataSet);
        // Strain
        if (strainIdentifier==null) {
            throw new BuildException("ERROR: could not extract strain identifier from "+readme.identifier+".");
        }
        strain = getDirectDataLoader().createObject(Strain.class);
        strain.setIdentifier(strainIdentifier);
        strain.setOrganism(organism);
        strain.addDataSets(dataSet);
    }
    
    /**
     * Process the FASTA file.
     */
    void processFasta(File file) throws BuildException {
        try {
            LinkedHashMap<String,DNASequence> sequenceMap = GZIPFastaReader.readFastaDNASequence(file);
            for (DNASequence sequence : sequenceMap.values()) {
                processSequence(sequence);
            }
        } catch (ParserException e) {
            throw new BuildException("Sequence not in FASTA format or wrong alphabet for: "+file, e);
        } catch (NoSuchElementException e) {
            throw new BuildException("No FASTA sequences in: "+file, e);
        } catch (ObjectStoreException e) {
            throw new BuildException("ObjectStore problem while processing: "+file, e);
        } catch (IOException e) {
            throw new BuildException("Problem reading FASTA: "+file, e);
        }
        // signal that FASTA has been processed since we use direct inserts without any lists here.
        fastaProcessed = true;
    }

    /**
     * Create a Sequence and an object for the given BioJava Sequence.
     *
     * @param bioJavaSequence the AbstractSequence object, either DNASequence or ProteinSequence
     * @throws ObjectStoreException if store() fails
     */
    void processSequence(AbstractSequence bioJavaSequence) throws ObjectStoreException {
        // the sequence string and MD5
        String sequence = bioJavaSequence.getSequenceAsString();
        String md5checksum = Util.getMd5checksum(sequence);
        // an InterMine Sequence
        Sequence bioSequence = getDirectDataLoader().createObject(org.intermine.model.bio.Sequence.class);
        bioSequence.setResidues(new PendingClob(sequence));
        bioSequence.setLength(bioJavaSequence.getLength());
        bioSequence.setMd5checksum(md5checksum);
        // the feature identifier
        String identifier = getIdentifier(bioJavaSequence);
        // HACK: don't allow spaces or tabs in primary identifiers; set symbol=extra part
        String symbol = null;
        String[] spaceChunks = identifier.split(" ");
        if (spaceChunks.length>1) {
            identifier = spaceChunks[0];
            symbol = spaceChunks[1];
        }
        String[] tabChunks = identifier.split("\t");
        if (tabChunks.length>1) {
            identifier = tabChunks[0];
            symbol = tabChunks[1];
        }
        // Use prefix match to identifier to set the class to Chromosome or Supercontig.
        if (isChromosome(identifier)) {
            // store Chromosome
            Class<? extends InterMineObject> imClass;
            Class<?> c;
            try {
                c = Class.forName("org.intermine.model.bio.Chromosome");
                imClass = (Class<? extends InterMineObject>) c;
                Chromosome feature = (Chromosome) getDirectDataLoader().createObject(imClass);
                feature.setPrimaryIdentifier(identifier);
                setSecondaryIdentifier(feature, identifier, false);
                setName(feature, bioJavaSequence, idAttribute, identifier, false);
                setAssemblyVersion(feature);
                storeSequenceFeature(feature, bioSequence);
            } catch (ClassNotFoundException ex) {
                throw new BuildException(ex);
            }
        } else if (isSupercontig(identifier)) {
            // store Supercontig
            Class<? extends InterMineObject> imClass;
            Class<?> c;
            try {
                c = Class.forName("org.intermine.model.bio.Supercontig");
                imClass = (Class<? extends InterMineObject>) c;
                Supercontig feature = (Supercontig) getDirectDataLoader().createObject(imClass);
                feature.setPrimaryIdentifier(identifier);
                setSecondaryIdentifier(feature, identifier, false);
                setName(feature, bioJavaSequence, idAttribute, identifier, false);
                setAssemblyVersion(feature);
                storeSequenceFeature(feature, bioSequence);
            } catch (ClassNotFoundException ex) {
                throw new BuildException(ex);
            }
        } else {
            System.out.println("## skipping non-matching sequence: "+identifier);
            return;
        }
    }
    
    /**
     * Store a SequenceFeature and its sequence.
     */
    void storeSequenceFeature(SequenceFeature feature, Sequence bioSequence) throws ObjectStoreException {
        feature.addDataSets(dataSet);
        feature.setOrganism(organism);
        feature.setStrain(strain);
        feature.setSequence(bioSequence);
        feature.setLength(bioSequence.getLength());
        if (publication!=null) feature.addPublications(publication);
        getDirectDataLoader().store(bioSequence);
        getDirectDataLoader().store(feature);
    }

    /**
     * For the given BioJava Sequence object, return an identifier to be used when creating the corresponding BioEntity.
     * @param bioJavaSequence the Sequence
     * @return an identifier
     */
    protected String getIdentifier(AbstractSequence bioJavaSequence) {
        String identifier = null;
        String header = bioJavaSequence.getAccession().getID();
        String[] bits = header.split(" ");
        if (bits[0].contains("|")) {
            String[] subbits = bits[0].split("\\|");
            identifier = subbits[1];
        } else {
            identifier = bits[0];
        }
        return identifier;
    }

    /**
     * For the given BioJava Sequence object, return the description, defined as everything that follows the first space in the header.
     * @param bioJavaSequence the Sequence
     * @return a description
     */
    protected String getDescription(AbstractSequence bioJavaSequence) {
        String description = null;
        String header = bioJavaSequence.getAccession().getID();
        String[] bits = header.split(" ");
        if (bits.length>1) {
            description = bits[1];
            for (int i=2; i<bits.length; i++) description += " "+bits[i];
            description = DatastoreUtils.unescape(description);
        }
        return description;
    }
    
    /**
     * Return the attribute with the given field name; else null.
     * 0                                                                      1a               1b                 1c                     1d
     * vigun.IT97K-499-35.gnm1.ann1.VigunL081000.1 pacid=39013057 polypeptide=VigunL081000.1.p locus=VigunL081000 ID=VigunL081000.1.v1.1 annot-version=v1.1
     * 0                                                            1a                       1b                        1c     1d                   1e   1f          1g
     * lupal.Amiga.gnm1.ann0.mRNA:Lalb_Chr00c01g0403611.1 locus_tag=Lalb_Chr00c01g0403611 gn=Lalb_Chr00c01g0403611 len=96 chr=Lalb_Chr00c01 strand=1 sp=Unknown def=Putative RNA
     *
     * @param bioJavaSequence the AbstractSequence
     * @param field the name of the desired field, e.g. protein_id
     * @return an identifier
     */
    static String getAttribute(AbstractSequence bioJavaSequence, String field) {
        String header = bioJavaSequence.getAccession().getID();
        if (!header.contains(field)) return null;
        String[] split = header.split(" "+field+"=");
        if (split.length==1) return null;
        String secondHalf = split[1];
        String[] parts = secondHalf.split(" ");
        return DatastoreUtils.unescape(parts[0]);
    }

    /**
     * Attribute with an identifier that is stored as the name, typically
     */
    public void setIdAttribute(String s) {
        idAttribute = s;
    }

    /**
     * Attribute with a description
     */
    public void setDescriptionAttribute(String s) {
        descriptionAttribute = s;
    }

    /**
     * Boolean attribute instructs us to load descriptions from FASTA headers in the following format:
     * >ID description here with spaces
     */
    public void setLoadHeaderDescriptions(String s) {
        loadHeaderDescriptions = s.equals("true");
    }

    /**
     * Set secondaryIdentifier from identifier or throw error.
     */
    static void setSecondaryIdentifier(BioEntity feature, String identifier, boolean isAnnotationFeature) {
        String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(identifier, isAnnotationFeature);
        if (secondaryIdentifier==null) {
            throw new BuildException("Could not get secondaryIdentifier from "+identifier+".");
        }
        feature.setSecondaryIdentifier(secondaryIdentifier);
    }

    /**
     * Set name to idAttribute if exists, else to secondaryIdentifier from identifier.
     */
    static void setName(BioEntity feature, AbstractSequence bioJavaSequence, String idAttribute, String identifier, boolean isAnnotationFeature) {
        if (idAttribute!=null) {
            String name = getAttribute(bioJavaSequence, idAttribute);
            if (name!=null) {
                feature.setName(name);
                return;
            }
        }
        String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(identifier, isAnnotationFeature);
        if (secondaryIdentifier==null) {
            throw new BuildException("Could not get secondaryIdentifier from "+identifier+".");
        }
        feature.setName(secondaryIdentifier);
    }

    /**
     * Set the assembly version or throw an exception.
     */
    void setAssemblyVersion(BioEntity feature) {
        if (assemblyVersion!=null) {
            feature.setAssemblyVersion(assemblyVersion);
        } else {
            throw new BuildException("Assembly version is not set from README.");
        }
    }

    /**
     * Set the annotation version or throw an exception.
     */
    void setAnnotationVersion(BioEntity feature) {
        if (annotationVersion!=null) {
            feature.setAnnotationVersion(annotationVersion);
        } else {
            throw new BuildException("Annotation version is not set from README.");
        }
    }

    /**
     * Return true if the given primaryIdentifier is for a Chromosome based on chromosome_prefix in the README.
     */
    public boolean isChromosome(String primaryIdentifier) {
        if (chromosomePrefixes.size()>0) {
            for (String prefix : chromosomePrefixes) {
                if (primaryIdentifier.startsWith(gensp+"."+strainIdentifier+"."+assemblyVersion+"."+prefix)) return true;
            }
        }
        return false;
    }

    /**
     * Return true if the given primaryIdentifier is for a Supercontig based on supercontig_prefix in the README.
     */
    public boolean isSupercontig(String primaryIdentifier) {
        if (supercontigPrefixes.size()>0) {
            for (String prefix : supercontigPrefixes) {
                if (primaryIdentifier.startsWith(gensp+"."+strainIdentifier+"."+assemblyVersion+"."+prefix)) return true;
            }
        }
        return false;
    }
}
