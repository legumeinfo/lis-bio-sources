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
import java.io.IOException;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.NoSuchElementException;

import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;

import org.biojava.nbio.core.exceptions.ParserException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import org.intermine.metadata.Util;
import org.intermine.model.InterMineObject;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.objectstore.query.PendingClob;
import org.intermine.task.FileDirectDataLoaderTask;

import org.intermine.model.bio.BioEntity;
import org.intermine.model.bio.CDS;
import org.intermine.model.bio.Chromosome;
import org.intermine.model.bio.DataSet;
import org.intermine.model.bio.DataSource;
import org.intermine.model.bio.MRNA;
import org.intermine.model.bio.Organism;
import org.intermine.model.bio.Publication;
import org.intermine.model.bio.Protein;
import org.intermine.model.bio.Strain;
import org.intermine.model.bio.Sequence;
import org.intermine.model.bio.SequenceFeature;
import org.intermine.model.bio.Supercontig;

import org.ncgr.datastore.Readme;

/**
 * A task that can read a set of FASTA files and create the corresponding Sequence objects in an ObjectStore.
 *
 * @author Kim Rutherford
 * @author Peter Mclaren
 * @author Sam Hokin
 */

public class LISFastaLoaderTask extends FileDirectDataLoaderTask {
    static final Logger LOG = Logger.getLogger(LISFastaLoaderTask.class);

    String sequenceType; // "dna" or "protein" based on file extension
    String className;
    boolean loadHeaderDescriptions = false;

    // project.xml setters
    String idAttribute, descriptionAttribute;
    String dataSourceName, dataSourceUrl, dataSourceDescription;
    String dataSetUrl, dataSetLicence;

    String gensp, assemblyVersion, annotationVersion;

    // collection items
    Organism organism;
    Strain strain;
    DataSource dataSource;
    DataSet dataSet;
    Publication publication;
    
    DatastoreUtils dsu; // for determining supercontigs vs chromosomes

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
        dsu = new DatastoreUtils();
        // process files, which stores the features and sequences
        super.process();
        try {
            // store collection items
            getDirectDataLoader().store(organism);
            getDirectDataLoader().store(strain);
            getDirectDataLoader().store(dataSource);
            getDirectDataLoader().store(dataSet);
            if (publication!=null) getDirectDataLoader().store(publication);
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
        // required project.xml parameters
        if (className==null || className.trim().length()==0) {
            throw new RuntimeException("className must be set in project.xml.");
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
        if (file.getName().startsWith("README")) {
            try {
                processReadme(file);
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            } catch (ObjectStoreException ex) {
                throw new RuntimeException(ex);
            }
        } else if (file.getName().endsWith(".fna") || file.getName().endsWith(".faa")) {
            // README must precede FASTA
            if (organism==null) {
                throw new RuntimeException("ERROR: README not read before FASTA file. Switch order in project.xml.");
            }
            if (file.getName().endsWith(".fna")) {
                sequenceType = "dna";
            } else {
                sequenceType = "protein";
            }
            System.out.println("# Reading "+sequenceType+" sequences from: "+file.getName());
            processFasta(file);
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
            readme.scientific_name_abbrev==null
            ) {
            throw new RuntimeException("ERROR in README: a required field is missing. "+
                                       "Required fields are: identifier, taxid, synopsis, description, scientific_name_abbrev");
        }
        String collection = DatastoreUtils.extractCollectionFromReadme(file);
        // needed for DatastoreUtils.isSupercontig()
        gensp = readme.scientific_name_abbrev;
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
        assemblyVersion = DatastoreUtils.extractAssemblyVersionFromCollection(readme.identifier);
        annotationVersion = DatastoreUtils.extractAnnotationVersionFromCollection(readme.identifier);
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
        // Strain
        String strainIdentifier = DatastoreUtils.extractStrainIdentifierFromCollection(collection);
        if (strainIdentifier==null) {
            throw new RuntimeException("ERROR: could not extract strain identifier from "+collection+".");
        }
        strain = getDirectDataLoader().createObject(Strain.class);
        strain.setIdentifier(strainIdentifier);
        strain.setOrganism(organism);
    }
    
    /**
     * Process the FASTA file.
     */
    void processFasta(File file) throws BuildException {
        try {
            if (sequenceType.equals("dna")) {
                LinkedHashMap<String,DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(file);
                for (DNASequence sequence : sequenceMap.values()) {
                    processSequence(sequence);
                }
            } else if (sequenceType.equals("protein")) {
                LinkedHashMap<String,ProteinSequence> sequenceMap = FastaReaderHelper.readFastaProteinSequence(file);
                for (ProteinSequence sequence : sequenceMap.values()) {
                    processSequence(sequence);
                }
            } else {
                throw new RuntimeException("Sequence type set in project.xml is neither dna nor protein.");
            }
        } catch (ParserException e) {
            throw new BuildException("Sequence not in FASTA format or wrong alphabet for: "+file, e);
        } catch (NoSuchElementException e) {
            throw new BuildException("No FASTA sequences in: "+file, e);
        } catch (FileNotFoundException e) {
            throw new BuildException("Problem reading file - file not found: "+file, e);
        } catch (ObjectStoreException e) {
            throw new BuildException("ObjectStore problem while processing: "+file, e);
        } catch (IOException e) {
            throw new BuildException("Error while closing FileReader for: "+file, e);
        }
    }

    /**
     * Create a Sequence and an object of type className for the given BioJava Sequence.
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
        // HACK: set the className to "Chromosome" or "Supercontig" based on identifier and identifying supercontig matching strings.
        if (className.equals("org.intermine.model.bio.Chromosome") || className.equals("org.intermine.model.bio.Supercontig")) {
            if (gensp==null) {
                throw new RuntimeException("ERROR: gensp is not set in README, or README not read before processSequence().");
            }
            if (dsu.isSupercontig(identifier)) {
                className = "org.intermine.model.bio.Supercontig";
            } else {
                className = "org.intermine.model.bio.Chromosome";
            }
        }
        // create the feature class
        Class<? extends InterMineObject> imClass;
        Class<?> c;
        try {
            c = Class.forName(className);
            if (InterMineObject.class.isAssignableFrom(c)) {
                imClass = (Class<? extends InterMineObject>) c;
            } else {
                throw new RuntimeException("Feature className must be a valid class in the model that inherits from InterMineObject, but was "+className);
            }
        } catch (ClassNotFoundException e1) {
            throw new RuntimeException("Unknown class: "+className+" while creating new Sequence object");
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // do the work separately for each class since some objects have special attributes to assign from header data //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if (className.equals("org.intermine.model.bio.Chromosome")) {
            // vigun.IT97K-499-35.gnm1.Vu01 (old4)
            Chromosome feature = (Chromosome) getDirectDataLoader().createObject(imClass);
	    // primaryIdentifier
            feature.setPrimaryIdentifier(identifier);
	    // secondaryIdentifier
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(identifier, false);
	    if (secondaryIdentifier!=null) feature.setSecondaryIdentifier(secondaryIdentifier);
            // name from ID attribute
            if (idAttribute!=null) {
                String name = getAttribute(bioJavaSequence, idAttribute);
                if (name!=null) feature.setName(name);
            }
            if (assemblyVersion!=null) feature.setAssemblyVersion(assemblyVersion);
            storeSequenceFeature(feature, bioSequence);
        } else if (className.equals("org.intermine.model.bio.Supercontig")) {
            // vigun.IT97K-499-35.gnm1.contig_700
            Supercontig feature = (Supercontig) getDirectDataLoader().createObject(imClass);
	    // primaryIdentifier
            feature.setPrimaryIdentifier(identifier);
            // secondaryIdentifier
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(identifier, false);
	    if (secondaryIdentifier!=null) feature.setSecondaryIdentifier(secondaryIdentifier);
	    // name from idAttribute
            if (idAttribute!=null) {
                String name = getAttribute(bioJavaSequence, idAttribute);
                if (name!=null) feature.setName(name);
            }
            if (assemblyVersion!=null) feature.setAssemblyVersion(assemblyVersion);
            storeSequenceFeature(feature, bioSequence);
        } else if (className.equals("org.intermine.model.bio.MRNA")) {
            MRNA feature = (MRNA) getDirectDataLoader().createObject(imClass);
	    // primaryIdentifier
	    feature.setPrimaryIdentifier(identifier);
            // secondaryIdentifier
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(identifier, true);
	    if (secondaryIdentifier!=null) feature.setSecondaryIdentifier(secondaryIdentifier);
	    // name from idAttribute
            if (idAttribute!=null) {
                String name = getAttribute(bioJavaSequence, idAttribute);
                if (name!=null) feature.setName(name);
            }
            if (assemblyVersion!=null) feature.setAssemblyVersion(assemblyVersion);
            if (annotationVersion!=null) feature.setAnnotationVersion(annotationVersion);
            storeSequenceFeature(feature, bioSequence);
        } else if (className.equals("org.intermine.model.bio.CDS")) {
            CDS feature = (CDS) getDirectDataLoader().createObject(imClass);
	    // primaryIdentifier
            feature.setPrimaryIdentifier(identifier);
            // secondaryIdentifier
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(identifier, true);
	    if (secondaryIdentifier!=null) feature.setSecondaryIdentifier(secondaryIdentifier);
	    // name from idAttribute
            if (idAttribute!=null) {
                String name = getAttribute(bioJavaSequence, idAttribute);
                if (name!=null) feature.setName(name);
            }
            if (assemblyVersion!=null) feature.setAssemblyVersion(assemblyVersion);
            if (annotationVersion!=null) feature.setAnnotationVersion(annotationVersion);
            // store with the sequence
            storeCDS(feature, bioSequence);
        } else if (className.equals("org.intermine.model.bio.Protein")) {
            // lupal.Amiga.gnm1.ann0.mRNA:Lalb_Chr00c01g0403611.1 locus_tag=Lalb_Chr00c01g0403611 gn=Lalb_Chr00c01g0403611 len=96 chr=Lalb_Chr00c01 strand=1 sp=Unknown
            // def=Putative RNA-directed DNA polymerase
            Protein feature = (Protein) getDirectDataLoader().createObject(imClass);
	    // primaryIdentifier
            feature.setPrimaryIdentifier(identifier);
            // secondaryIdentifier
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(identifier, true);
	    if (secondaryIdentifier!=null) feature.setSecondaryIdentifier(secondaryIdentifier);
	    // name from idAttribute
            if (idAttribute!=null) {
                String name = getAttribute(bioJavaSequence, idAttribute);
                if (name!=null) feature.setName(name);
            }
            // description
            if (loadHeaderDescriptions) {
                String description = getDescription(bioJavaSequence);
                if (description!=null) feature.setDescription(description);
            } else if (descriptionAttribute!=null) {
                String description = getAttribute(bioJavaSequence, descriptionAttribute);
                if (description!=null) feature.setDescription(description);
            }
            if (assemblyVersion!=null) feature.setAssemblyVersion(assemblyVersion);
            if (annotationVersion!=null) feature.setAnnotationVersion(annotationVersion);
            // store with the sequence
            storeProtein(feature, bioSequence);
        } else {
            throw new RuntimeException("Loading of "+className+" from FASTA isn't currently supported.");
        }
    }

    /**
     * Store a CDS and its sequence.
     */
    void storeCDS(CDS feature, Sequence bioSequence) throws ObjectStoreException {
        feature.addDataSets(dataSet);
        feature.setOrganism(organism);
        feature.setStrain(strain);
        feature.setSequence(bioSequence);
        feature.setLength(bioSequence.getLength());
        getDirectDataLoader().store(bioSequence);
        getDirectDataLoader().store(feature);
    }

    /**
     * Store a Protein and its sequence.
     */
    void storeProtein(Protein feature, Sequence bioSequence) throws ObjectStoreException {
        feature.addDataSets(dataSet);
        feature.setOrganism(organism);
        feature.setStrain(strain);
        feature.setSequence(bioSequence);
        feature.setLength(bioSequence.getLength());
        getDirectDataLoader().store(bioSequence);
        getDirectDataLoader().store(feature);
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
    protected String getAttribute(AbstractSequence bioJavaSequence, String field) {
        String header = bioJavaSequence.getAccession().getID();
        if (!header.contains(field)) return null;
        String[] split = header.split(" "+field+"=");
        if (split.length==1) return null;
        String secondHalf = split[1];
        String[] parts = secondHalf.split(" ");
        return parts[0];
    }

    /**
     * The default class name to use for objects created during load.  Generally this is
     * "org.intermine.model.bio.Chromosome" or "org.intermine.model.bio.Protein"
     * @param className the class name
     */
    public void setClassName(String className) {
        this.className = className;
    }

    /**
     * Return the class name set with setClassName().
     * @return the class name
     */
    public String getClassName() {
        return className;
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
}
