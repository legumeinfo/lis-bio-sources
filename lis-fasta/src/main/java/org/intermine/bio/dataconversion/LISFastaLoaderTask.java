package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2019 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
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
import org.intermine.model.bio.Protein;
import org.intermine.model.bio.Strain;
import org.intermine.model.bio.Sequence;
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

    String taxonId;
    String gensp;
    String strainIdentifier;
    String assemblyVersion;   // extracted from filename
    String annotationVersion; // extracted from filename
    String className;

    File fastaFile;           // the FASTA file being processed

    DataSource dataSource;
    String dataSourceName, dataSourceUrl, dataSourceDescription;

    DataSet dataSet;
    String dataSetName, dataSetUrl, dataSetDescription, dataSetVersion, dataSetLicence;

    Organism organism;
    Strain strain;

    String idAttribute, geneAttribute, proteinAttribute, transcriptAttribute, descriptionAttribute;

    Map<String,Gene> genes = new HashMap<>();
    Map<String,Protein> proteins = new HashMap<>();
    Map<String,MRNA> mRNAs = new HashMap<>();

    DatastoreUtils datastoreUtils; // for determining supercontigs

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
     * Set the organism taxon ID and gensp.
     */
    public void setTaxonId(String taxonId) {
        this.taxonId = taxonId;
    }

    /**
     * Set the strain identifier.
     */
    public void setStrainIdentifier(String strainIdentifier) {
        this.strainIdentifier = strainIdentifier;
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

    /**
     * Return the full-yuck prefix from the dataSetName, like vigun.IT97K-499-35.gnm1.ann1.zb5D or phalu.G27455.gnm1.7NXX.genome_main
     */
    public String getYuckyPrefix() {
        return dataSetName.substring(0, dataSetName.length()-5);
    }

    /**
     * Process and load all of the fasta files.
     */
    @Override
    public void process() {
        try {
	    datastoreUtils = new DatastoreUtils();
            gensp = datastoreUtils.getGensp(taxonId);
            super.process();
            getIntegrationWriter().commitTransaction();
            getIntegrationWriter().beginTransaction();
            getDirectDataLoader().close();
        } catch (ObjectStoreException e) {
            throw new BuildException("Failed to store object", e);
        }
    }

    /**
     * Be sure to close the data loader so the last batch gets stored. only needed for tests
     * since the data loading task usually does that for hte live builds.
     * @throws ObjectStoreException if we can't store to db
     */
    public void close() throws ObjectStoreException {
        // store any data left over
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
     * Handles each fasta file. Factored out so we can supply files for testing.
     *
     * @param file the File to process.
     * @throws BuildException if the is a problem
     */
    @Override
    public void processFile(File file) {
	String blurb = "Reading "+sequenceType+" sequences from: "+file.getName();
	String hashes = "";
	for (int i=0; i<blurb.length(); i++) hashes += "#";
        System.out.println(hashes);
        System.out.println(blurb);
        System.out.println(hashes);
        LOG.info("LISFastaLoaderTask loading file "+file.getName());
	this.fastaFile = file;
	assemblyVersion = DatastoreUtils.extractAssemblyVersion(file.getName());
	annotationVersion = DatastoreUtils.extractAnnotationVersion(file.getName());
        try {
            // process FASTA file
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
     * Get and store() the Organism object to reference when creating new objects.
     * @param bioJavaSequence the biojava sequence to be parsed
     * @throws ObjectStoreException if there is a problem
     * @return the new Organism
     */
    protected Organism getOrganism() throws ObjectStoreException {
        if (organism==null) {
            organism = getDirectDataLoader().createObject(Organism.class);
            organism.setTaxonId(taxonId);
            organism.setAbbreviation(gensp);
            getDirectDataLoader().store(organism);
        }
        return organism;
    }

    /**
     * Get and store() the Strain object to reference when creating new objects.
     * @param bioJavaSequence the biojava sequence to be parsed (not used)
     * @throws ObjectStoreException if there is a problem
     * @return the new Strain
     */
    protected Strain getStrain() throws ObjectStoreException {
        if (strain==null && strainIdentifier!=null) {
            strain = getDirectDataLoader().createObject(Strain.class);
            strain.setIdentifier(strainIdentifier);
            strain.setOrganism(getOrganism());
            getDirectDataLoader().store(strain);
        }
        return strain;
    }

    /**
     * Get/add the DataSet Item. dataSetUrl and dataSetDescription MUST be set in project.xml.
     * @return the DataSet
     * @throws ObjectStoreException if there is an ObjectStore problem
     */
    DataSet getDataSet() throws ObjectStoreException {
	if (dataSet==null) {
	    if (dataSetUrl==null || dataSetDescription==null || dataSetVersion==null) {
		throw new RuntimeException("You must set lis-fasta.dataSetUrl, lis-fasta.dataSetDescription and lis-fasta.dataSetVersion in project.xml.");
	    }
	    dataSet = getDirectDataLoader().createObject(DataSet.class);
            dataSet.setName(fastaFile.getName());
            dataSet.setDataSource(getDataSource());
            dataSet.setUrl(dataSetUrl);
            dataSet.setDescription(dataSetDescription);
	    dataSet.setVersion(dataSetVersion);
            getDirectDataLoader().store(dataSet);
	}
	return dataSet;
    }

    /**
     * Store and/or return the DataSource set in project.xml.
     * @return the DataSource
     * @throws ObjectStoreException if there is a problem storing it.
     */
    DataSource getDataSource() throws ObjectStoreException {
        if (dataSource==null) {
            dataSource = getDirectDataLoader().createObject(DataSource.class);
	    if (dataSourceName==null) {
		dataSourceName = DatastoreFileConverter.DEFAULT_DATASOURCE_NAME;
	    }
	    if (dataSourceUrl==null) {
		dataSourceUrl = DatastoreFileConverter.DEFAULT_DATASOURCE_URL;
	    }
	    if (dataSourceDescription==null) {
		dataSourceDescription = DatastoreFileConverter.DEFAULT_DATASOURCE_DESCRIPTION;
	    }
            dataSource.setName(dataSourceName);
            dataSource.setUrl(dataSourceUrl);
            dataSource.setDescription(dataSourceDescription);
            getDirectDataLoader().store(dataSource);
        }
        return dataSource;
    }

    /**
     * Store and/or return the Gene given by the identifier.
     * @return the Gene
     * @throws ObjectStoreException
     */
    Gene getGene(String identifier) throws ObjectStoreException {
        // HACK: Prepend the full-yuck prefix to the identifier.
        identifier = getYuckyPrefix() + "." + identifier;
        if (genes.containsKey(identifier)) {
            return genes.get(identifier);
        } else {
            Gene gene = getDirectDataLoader().createObject(Gene.class);
            gene.setPrimaryIdentifier(identifier);
            getDirectDataLoader().store(gene);
            genes.put(identifier, gene);
            return gene;
        }
    }

    /**
     * Store and/or return the MRNA given by the identifier, plus set the protein reference.
     * @return the MRNA
     * @throws ObjectStoreException
     */
    MRNA getMRNA(String identifier, Protein protein) throws ObjectStoreException {
        // HACK: Prepend the full-yuck prefix to the identifier.
        identifier = getYuckyPrefix() + "." + identifier;
        if (mRNAs.containsKey(identifier)) {
            return mRNAs.get(identifier);
        } else {
            MRNA mRNA = getDirectDataLoader().createObject(MRNA.class);
            mRNA.setPrimaryIdentifier(identifier);
            mRNA.setProtein(protein);
            getDirectDataLoader().store(mRNA);
            mRNAs.put(identifier, mRNA);
            return mRNA;
        }
    }

    /**
     * Store and/or return the Protein given by the identifier.
     * @return the Protein
     * @throws ObjectStoreException
     */
    Protein getProtein(String identifier) throws ObjectStoreException {
        // HACK 1: Drop stupid JGI .p suffixes.
        if (identifier.endsWith(".p")) identifier = identifier.substring(0, identifier.length()-2);
        // HACK 2: Prepend the full-yuck prefix.
        identifier = getYuckyPrefix() + "." + identifier;
        if (proteins.containsKey(identifier)) {
            return proteins.get(identifier);
        } else {
            Protein protein = getDirectDataLoader().createObject(Protein.class);
            protein.setPrimaryIdentifier(identifier);
            getDirectDataLoader().store(protein);
            proteins.put(identifier, protein);
            return protein;
        }
    }

    /**
     * Store and/or return the Protein given by the identifier, plus add the given Gene to the genesproteins collection.
     * @return the Protein
     * @throws ObjectStoreException
     */
    Protein getProtein(String identifier, Gene gene) throws ObjectStoreException {
        // HACK 1: Drop stupid JGI .p suffixes.
        if (identifier.endsWith(".p")) identifier = identifier.substring(0, identifier.length()-2);
        // HACK 2: Prepend the full-yuck prefix.
        identifier = getYuckyPrefix() + "." + identifier;
        if (proteins.containsKey(identifier)) {
            return proteins.get(identifier);
        } else {
            Protein protein = getDirectDataLoader().createObject(Protein.class);
            protein.setPrimaryIdentifier(identifier);
            protein.addGenes(gene);
            getDirectDataLoader().store(protein);
            proteins.put(identifier, protein);
            return protein;
        }
    }
    
    /**
     * Create a Sequence and an object of type className for the given BioJava Sequence.
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

	// feature name -- not being used
        // String name = getName(bioJavaSequence);

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
            store(feature, bioSequence);

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
            store(feature, bioSequence);

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
            store(feature, bioSequence);

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
            store(feature, bioSequence);

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
            store(feature, bioSequence);
        } else {
            throw new RuntimeException("Loading of "+className+" from FASTA isn't currently supported.");
        }
    }

    /**
     * Store a feature and its sequence.
     */
    void store(BioEntity feature, Sequence bioSequence) throws ObjectStoreException {
        feature.setFieldValue("sequence", bioSequence);
        feature.setFieldValue("length", new Integer(bioSequence.getLength()));
        getDirectDataLoader().store(bioSequence);
        getDirectDataLoader().store(feature);
    }
    
    /**
     * Set the attributes and references common to all LIS BioEntity objects.
     */
    void setCommonAttributes(BioEntity feature) throws ObjectStoreException {
        feature.addDataSets(getDataSet());
        feature.setOrganism(getOrganism());
        if (strainIdentifier!=null) feature.setStrain(getStrain());
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
     * For the given BioJava Sequence object, return the name attribute.
     * @param bioJavaSequence the Sequence
     * @return a name
     */
    protected String getName(AbstractSequence bioJavaSequence) {
        String header = bioJavaSequence.getAccession().getID();
        String[] bits = header.split(" ");
        if (bits.length==1) {
            return null;
        } else {
            String name = "";
            for (int i=1; i<bits.length; i++) {
                if (bits[i].contains("[")) {
                    break;
                } else {
                    name += bits[i]+" ";
                }
            }
            name = name.trim();
            if (name.length()==0) {
                return null;
            } else {
                return name;
            }
        }
    }
}
