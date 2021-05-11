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
import org.intermine.model.bio.Gene;
import org.intermine.model.bio.MRNA;
import org.intermine.model.bio.Organism;
import org.intermine.model.bio.Publication;
import org.intermine.model.bio.Protein;
import org.intermine.model.bio.Strain;
import org.intermine.model.bio.Sequence;
import org.intermine.model.bio.SequenceFeature;
import org.intermine.model.bio.Supercontig;

/**
 * A task that can read a set of FASTA files and create the corresponding Sequence objects in an ObjectStore.
 *
 * @author Kim Rutherford
 * @author Peter Mclaren
 * @author Sam Hokin
 */

public class LISFastaLoaderTask extends FileDirectDataLoaderTask {
    static final Logger LOG = Logger.getLogger(LISFastaLoaderTask.class);

    String sequenceType = "dna"; // default, or "protein"

    String gensp;             // extracted from filename
    String strainIdentifier;  // extracted from filename
    String assemblyVersion;   // extracted from filename
    String annotationVersion; // extracted from filename
    String className;         // set in project.xml

    String dataSetUrl, dataSetVersion;

    // global objects to be stored at end
    Organism organism;
    Strain strain;
    DataSource dataSource;
    DataSet dataSet;
    Publication publication;
    
    String idAttribute, geneAttribute, proteinAttribute, transcriptAttribute, descriptionAttribute;

    Map<String,Gene> genes = new HashMap<>();
    Map<String,Protein> proteins = new HashMap<>();
    Map<String,MRNA> mRNAs = new HashMap<>();

    DatastoreUtils datastoreUtils; // for determining supercontigs vs chromosomes

    /**
     * Process and load all of the fasta files.
     */
    @Override
    public void process() {
        datastoreUtils = new DatastoreUtils();
        try {
            // Organism
            organism = getDirectDataLoader().createObject(org.intermine.model.bio.Organism.class);
            // Strain
            strain = getDirectDataLoader().createObject(org.intermine.model.bio.Strain.class);
            strain.setIdentifier(strainIdentifier);
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
        } catch (ObjectStoreException ex) {
            throw new BuildException(ex);
        }
        // process files
        super.process();
        // wrap up globals
        try {
            getDirectDataLoader().store(organism);
            getDirectDataLoader().store(strain);
            getDirectDataLoader().store(dataSource);
            getDirectDataLoader().store(dataSet);
            if (publication!=null) getDirectDataLoader().store(publication);
            // DO WE NEED THESE?
            // getIntegrationWriter().commitTransaction();
            // getIntegrationWriter().beginTransaction();
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
        // required project.xml parameters
        if (className==null || className.trim().length()==0) {
            throw new RuntimeException("className must be set in project.xml.");
        }
        // optional but recommended project.xml parameters
        if (strainIdentifier==null || strainIdentifier.trim().length()==0) {
            System.out.println("NOTE: strainIdentifier is missing in project.xml.");
        }
        // this will call processFile() for each file
        super.execute();
    }

    /**
     * Handles each file.
     *
     * @param file the File to process.
     * @throws BuildException if the is a problem
     */
    @Override
    public void processFile(File file) throws BuildException {
        if (file.getName().startsWith("README")) {
            System.out.println("Reading metadata from "+file.getName());
            LOG.info("LISFastaLoaderTask loading file "+file.getName());
            processREADME(file);
        } else if (file.getName().endsWith(".fna") || file.getName().endsWith(".faa")) {
            System.out.println("Reading "+sequenceType+" sequences from: "+file.getName());
            LOG.info("LISFastaLoaderTask loading file "+file.getName());
            processFasta(file);
        }
    }

    /**
     * Process the README, which contains metadata.
     */
    void processREADME(File file) throws BuildException {
        try {
            Readme readme = Readme.getReadme(file);
            // Organism
            organism.setTaxonId(readme.taxid);
            // Strain
            strain.setOrganism(organism);
            // DataSet
            dataSet.setName(readme.identifier);
            dataSet.setDescription(readme.description);
            // Publication
            if (readme.publication_doi!=null) {
                publication = getDirectDataLoader().createObject(Publication.class);
                publication.setDoi(readme.publication_doi);
                if (readme.publication_title!=null) publication.setTitle(readme.publication_title);
            }
        } catch (IOException ex) {
            throw new BuildException(ex);
        } catch (ObjectStoreException ex) {
            throw new BuildException(ex);
        }
    }

    /**
     * Process the FASTA file.
     */
    void processFasta(File file) throws BuildException {
        gensp = DatastoreUtils.extractGensp(file.getName());
        assemblyVersion = DatastoreUtils.extractAssemblyVersion(file.getName());
        annotationVersion = DatastoreUtils.extractAnnotationVersion(file.getName());
        organism.setAbbreviation(gensp);
        try {
            if (sequenceType.equalsIgnoreCase("dna")) {
                LinkedHashMap<String,DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(file);
                for (DNASequence sequence : sequenceMap.values()) {
                    processSequence(sequence);
                }
            } else if (sequenceType.equalsIgnoreCase("protein")) {
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
            if (datastoreUtils.isSupercontig(gensp,strainIdentifier,identifier)) {
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
            setCommonAttributes(feature);
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
            // assemblyVersion
            if (assemblyVersion!=null) feature.setAssemblyVersion(assemblyVersion);
            storeSequenceFeature(feature, bioSequence);
        } else if (className.equals("org.intermine.model.bio.Supercontig")) {
            // vigun.IT97K-499-35.gnm1.contig_700
            Supercontig feature = (Supercontig) getDirectDataLoader().createObject(imClass);
            setCommonAttributes(feature);
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
            // assemblyVersion
            if (assemblyVersion!=null) feature.setAssemblyVersion(assemblyVersion);
            storeSequenceFeature(feature, bioSequence);
        } else if (className.equals("org.intermine.model.bio.MRNA")) {
            MRNA feature = (MRNA) getDirectDataLoader().createObject(imClass);
            setCommonAttributes(feature);
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
            // assemblyVersion
            if (assemblyVersion!=null) feature.setAssemblyVersion(assemblyVersion);
            // annotationVersion
            if (annotationVersion!=null) feature.setAnnotationVersion(annotationVersion);
            storeSequenceFeature(feature, bioSequence);
        } else if (className.equals("org.intermine.model.bio.CDS")) {
            CDS feature = (CDS) getDirectDataLoader().createObject(imClass);
            setCommonAttributes(feature);
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
            // assemblyVersion
            if (assemblyVersion!=null) feature.setAssemblyVersion(assemblyVersion);
            // annotationVersion
            if (annotationVersion!=null) feature.setAnnotationVersion(annotationVersion);
            storeSequenceFeature(feature, bioSequence);
        } else if (className.equals("org.intermine.model.bio.Protein")) {
            // lupal.Amiga.gnm1.ann0.mRNA:Lalb_Chr00c01g0403611.1 locus_tag=Lalb_Chr00c01g0403611 gn=Lalb_Chr00c01g0403611 len=96 chr=Lalb_Chr00c01 strand=1 sp=Unknown
            // def=Putative RNA-directed DNA polymerase
            Protein feature = (Protein) getDirectDataLoader().createObject(imClass);
            setCommonAttributes(feature);
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
            if (descriptionAttribute!=null) {
                String description = getAttribute(bioJavaSequence, descriptionAttribute);
                if (description!=null) feature.setDescription(description);
            }
            // assemblyVersion
            if (assemblyVersion!=null) feature.setAssemblyVersion(assemblyVersion);
            // annotationVersion
            if (annotationVersion!=null) feature.setAnnotationVersion(annotationVersion);
            storeProtein(feature, bioSequence);
        } else {
            throw new RuntimeException("Loading of "+className+" from FASTA isn't currently supported.");
        }
    }

    /**
     * Store a Protein and its sequence.
     */
    void storeProtein(Protein feature, Sequence bioSequence) throws ObjectStoreException {
        feature.setSequence(bioSequence);
        feature.setLength(new Integer(bioSequence.getLength()));
        getDirectDataLoader().store(bioSequence);
        getDirectDataLoader().store(feature);
    }
    
    /**
     * Store a SequenceFeature and its sequence.
     */
    void storeSequenceFeature(SequenceFeature feature, Sequence bioSequence) throws ObjectStoreException {
        feature.setSequence(bioSequence);
        feature.setLength(new Integer(bioSequence.getLength()));
        getDirectDataLoader().store(bioSequence);
        getDirectDataLoader().store(feature);
    }

    /**
     * Set the attributes and references common to all LIS BioEntity objects.
     */
    void setCommonAttributes(BioEntity feature) throws ObjectStoreException {
        feature.addDataSets(dataSet);
        feature.setOrganism(organism);
        feature.setStrain(strain);
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
     * Set the sequence type to be passed to the FASTA parser.  The default is "dna".
     * @param sequenceType the sequence type
     */
    public void setSequenceType(String sequenceType) {
        if ("${fasta.sequenceType}".equals(sequenceType)) {
            this.sequenceType = "dna";
        } else {
            this.sequenceType = sequenceType;
        }
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
     * If a value is specified this url will used when a DataSet is created.
     * @param dataSetUrl the title of the DataSet of any new features
     */
    public void setDataSetUrl(String value) throws ObjectStoreException {
        this.dataSetUrl = value;
    }

    /**
     * If a value is specified this description will used when a DataSet is created.
     * @param dataSetVersion the version of the DataSet
     */
    public void setDataSetVersion(String value) throws ObjectStoreException {
        this.dataSetVersion = value;
    }

    /**
     * Set the strain identifier.
     */
    public void setStrainIdentifier(String value) throws ObjectStoreException {
        this.strainIdentifier = value;
    }

    /**
     * Attribute with an identifier that is stored as the name, typically
     */
    public void setIdAttribute(String s) {
        idAttribute = s;
    }

    /**
     * Attribute with a gene identifier
     */
    public void setGeneAttribute(String s) {
        geneAttribute = s;
    }
    /**
     * Attribute with a protein identifier
     */
    public void setProteinAttribute(String s) {
        proteinAttribute = s;
    }
    /**
     * Attribute with a transcript identifier
     */
    public void setTranscriptAttribute(String s) {
        transcriptAttribute = s;
    }
    /**
     * Attribute with a description
     */
    public void setDescriptionAttribute(String s) {
        descriptionAttribute = s;
    }
}
