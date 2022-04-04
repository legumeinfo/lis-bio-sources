package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.Reader;
import java.io.IOException;
import java.io.InputStream;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.HashSet;
import java.util.Properties;
import static java.util.Map.entry;

import org.apache.log4j.Logger;

import org.biojava.nbio.core.exceptions.ParserException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.Util;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.genome.parsers.gff.Location;

/**
 * Loads data from an LIS datastore annotations collection.
 * Files types loaded are:
 *   README
 *   gene_models_main.gff3
 *   protein.faa
 *   cds.fna
 *   mrna.fna
 *   gfa.tsv
 *   pathway.tsv
 *
 * @author Sam Hokin
 */
public class AnnotationFileConverter extends DatastoreFileConverter {
	
    static final String DEFAULT_VERSION = "legfed_v1_0";

    private static final Logger LOG = Logger.getLogger(AnnotationFileConverter.class);

    // GFF sourced
    Map<String,Item> chromosomes = new HashMap<>();
    Map<String,Item> supercontigs = new HashMap<>();
    Map<String,Item> features = new HashMap<>();       // features other than genes and mRNAs
    Map<String,Item> genes = new HashMap<>();          // also GFA, pathway, etc.
    Map<String,Item> mRNAs = new HashMap<>();          // mRNA sequences,length are FASTA-sourced
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> proteinDomains = new HashMap<>();
    List<Item> ontologyAnnotations = new ArrayList<>();
    List<Item> locations = new ArrayList<>();

    // TSV sourced
    Map<String,Item> geneFamilies = new HashMap<>();
    Map<String,Item> pathways = new HashMap<>();

    // FASTA sourced
    List<Item> sequences = new ArrayList<>();
    Map<String,Item> proteins = new HashMap<>(); // also GFA, etc.
    Map<String,Item> cdses = new HashMap<>();

    // for distinguishing chromosomes from supercontigs
    DatastoreUtils dsu;

    // map GFF types to InterMine classes; be sure to include extras in the additions file!
    Map<String,String> featureClasses = Map.ofEntries(entry("gene", "Gene"),
                                                      entry("mRNA", "MRNA"),
                                                      entry("CDS", "CDSRegion"),
                                                      entry("exon", "Exon"),
                                                      entry("three_prime_UTR", "ThreePrimeUTR"),
                                                      entry("five_prime_UTR", "FivePrimeUTR"),
                                                      entry("lnc_RNA", "LncRNA"),
                                                      entry("transcript", "Transcript"),
                                                      entry("pseudogene", "Pseudogene"),
                                                      entry("primary_transcript", "Transcript"),
                                                      entry("miRNA", "MiRNA"),
                                                      entry("miRNA_primary_transcript", "MiRNA"),
                                                      entry("tRNA", "TRNA"),
                                                      entry("tRNA_primary_transcript", "TRNA"),
                                                      entry("snoRNA", "SnoRNA"),
                                                      entry("snRNA", "SnRNA"),
                                                      entry("rRNA", "RRNA"),
                                                      entry("rRNA_primary_transcript", "RRNA"));

    /**
     * Create a new AnnotationFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public AnnotationFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
        dsu = new DatastoreUtils();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
            setStrain();
        } else if (getCurrentFile().getName().endsWith(".gene_models_main.gff3")) {
            System.out.println("Processing "+getCurrentFile().getName());
            processGFF3File();
        } else if (getCurrentFile().getName().endsWith(".gfa.tsv")) {
            System.out.println("Processing "+getCurrentFile().getName());
            processGFAFile(reader);
        } else if (getCurrentFile().getName().endsWith(".pathway.tsv")) {
            System.out.println("Processing "+getCurrentFile().getName());
            processPathwayFile(reader);
	} else if (getCurrentFile().getName().endsWith(".protein.faa") || getCurrentFile().getName().endsWith(".protein_primary.faa")) {
            System.out.println("Processing "+getCurrentFile().getName());
            processProteinFasta();
        } else if (getCurrentFile().getName().endsWith(".cds.fna") || getCurrentFile().getName().endsWith(".cds_primary.fna")) {
            System.out.println("Processing "+getCurrentFile().getName());
            processCDSFasta();
        } else if (getCurrentFile().getName().endsWith(".mrna.fna") || getCurrentFile().getName().endsWith(".mrna_primary.fna")) {
            System.out.println("Processing "+getCurrentFile().getName());
            processMRNAFasta();
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException, RuntimeException {
        if (readme==null) {
            throw new RuntimeException("README file not read. Aborting.");
        }
        // set Datastore collection attributes and references for annotation-specific objects
        for (Item chromosome : chromosomes.values()) {
            chromosome.setAttribute("assemblyVersion", assemblyVersion);
            chromosome.setReference("organism", organism);
            chromosome.setReference("strain", strain);
        }
        for (Item supercontig : supercontigs.values()) {
            supercontig.setAttribute("assemblyVersion", assemblyVersion);
            supercontig.setReference("organism", organism);
            supercontig.setReference("strain", strain);
        }
        for (Item feature: features.values()) {
            feature.setAttribute("assemblyVersion", assemblyVersion);
            feature.setAttribute("annotationVersion", annotationVersion);
            feature.setReference("organism", organism);
            feature.setReference("strain", strain);
        }            
        for (Item cds : cdses.values()) {
            cds.setAttribute("assemblyVersion", assemblyVersion);
            cds.setAttribute("annotationVersion", annotationVersion);
            cds.setReference("organism", organism);
            cds.setReference("strain", strain);
        }
        for (Item gene : genes.values()) {
            gene.setAttribute("assemblyVersion", assemblyVersion);
            gene.setAttribute("annotationVersion", annotationVersion);
            gene.setReference("organism", organism);
            gene.setReference("strain", strain);
        }
        for (Item mrna : mRNAs.values()) {
            mrna.setAttribute("assemblyVersion", assemblyVersion);
            mrna.setAttribute("annotationVersion", annotationVersion);
            mrna.setReference("organism", organism);
            mrna.setReference("strain", strain);
        }
        for (Item protein : proteins.values()) {
            protein.setAttribute("assemblyVersion", assemblyVersion);
            protein.setAttribute("annotationVersion", annotationVersion);
            protein.setReference("organism", organism);
            protein.setReference("strain", strain);
        }
        // set references and collections for objects loaded from FASTAs based on matching identifiers
        // Note: protein<->transcript is not in SO but set here given general practice
        for (String primaryIdentifier : cdses.keySet()) {
            Item cds = cdses.get(primaryIdentifier);
            if (mRNAs.containsKey(primaryIdentifier)) cds.setReference("transcript", mRNAs.get(primaryIdentifier));
            if (proteins.containsKey(primaryIdentifier)) cds.setReference("protein", proteins.get(primaryIdentifier));
            String geneId = getGeneIdFromDotId(primaryIdentifier);
            if (genes.containsKey(geneId)) cds.setReference("gene", genes.get(geneId));
        }
        for (String primaryIdentifier : mRNAs.keySet()) {
            Item mRNA = mRNAs.get(primaryIdentifier);
            if (cdses.containsKey(primaryIdentifier)) mRNA.addToCollection("CDSs", cdses.get(primaryIdentifier));
            if (proteins.containsKey(primaryIdentifier)) mRNA.setReference("protein", proteins.get(primaryIdentifier));
            String geneId = getGeneIdFromDotId(primaryIdentifier);
            if (genes.containsKey(geneId)) mRNA.setReference("gene", genes.get(geneId));
        }
        for (String primaryIdentifier : proteins.keySet()) {
            Item protein = proteins.get(primaryIdentifier);
            if (cdses.containsKey(primaryIdentifier)) protein.setReference("CDS", cdses.get(primaryIdentifier));
            if (mRNAs.containsKey(primaryIdentifier)) protein.setReference("transcript", mRNAs.get(primaryIdentifier));
            String geneId = getGeneIdFromDotId(primaryIdentifier);
            if (genes.containsKey(geneId)) protein.addToCollection("genes", genes.get(geneId));
        }
        // store
        storeCollectionItems();
        store(chromosomes.values());
        store(supercontigs.values());
        store(features.values());
        store(cdses.values());
        store(genes.values());
        store(mRNAs.values());
        store(proteins.values());
        store(geneFamilies.values());
        store(ontologyAnnotations);
        store(ontologyTerms.values());
        store(pathways.values());
        store(proteinDomains.values());
        store(locations);
        store(sequences);
    }

    /**
     * Process a protein FASTA faa file.
     */
    void processProteinFasta() throws IOException {
        boolean isPrimary = getCurrentFile().getName().endsWith(".protein_primary.faa");
        LinkedHashMap<String,ProteinSequence> sequenceMap = FastaReaderHelper.readFastaProteinSequence(getCurrentFile());
        for (ProteinSequence bioJavaSequence : sequenceMap.values()) {
            String residues = bioJavaSequence.getSequenceAsString();
            String md5checksum = Util.getMd5checksum(residues);
            Item sequence = createItem("Sequence");
            sequence.setAttribute("residues", residues);
            sequence.setAttribute("length", String.valueOf(residues.length()));
            sequence.setAttribute("md5checksum", md5checksum);
            sequences.add(sequence);
            String identifier = getFastaIdentifier(bioJavaSequence);
            // HACK: don't allow spaces or tabs in primary identifiers; set symbol=extra part if present
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
            // Protein Item
            Item protein = getProtein(identifier);
            if (isPrimary) protein.setAttribute("isPrimary", "true");
            protein.setReference("sequence", sequence);
            protein.setAttribute("length", String.valueOf(residues.length()));
            protein.setAttribute("md5checksum", md5checksum);
            if (symbol!=null) protein.setAttribute("symbol", symbol);
        }
    }

    /**
     * Process a CDS nucleotide FASTA fna file
     */
    void processCDSFasta() throws IOException {
        boolean isPrimary = getCurrentFile().getName().endsWith(".cds_primary.fna");
        LinkedHashMap<String,DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(getCurrentFile());
        for (DNASequence bioJavaSequence : sequenceMap.values()) {
            String residues = bioJavaSequence.getSequenceAsString();
            String md5checksum = Util.getMd5checksum(residues);
            Item sequence = createItem("Sequence");
            sequence.setAttribute("residues", residues);
            sequence.setAttribute("length", String.valueOf(residues.length()));
            sequence.setAttribute("md5checksum", md5checksum);
            sequences.add(sequence);
            String identifier = getFastaIdentifier(bioJavaSequence);
            // HACK: don't allow spaces or tabs in primary identifiers; set symbol=extra part if present
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
            // CDS Item
            Item cds = getCDS(identifier);
            if (isPrimary) cds.setAttribute("isPrimary", "true");
            cds.setReference("sequence", sequence);
            cds.setAttribute("length", String.valueOf(residues.length()));
            if (symbol!=null) cds.setAttribute("symbol", symbol);
        }
    }

    /**
     * Process an mRNA nucleotide FASTA fna file
     */
    void processMRNAFasta() throws IOException {
        boolean isPrimary = getCurrentFile().getName().endsWith(".mrna_primary.fna");
        LinkedHashMap<String,DNASequence> sequenceMap = FastaReaderHelper.readFastaDNASequence(getCurrentFile());
        for (DNASequence bioJavaSequence : sequenceMap.values()) {
            String residues = bioJavaSequence.getSequenceAsString();
            String md5checksum = Util.getMd5checksum(residues);
            Item sequence = createItem("Sequence");
            sequence.setAttribute("residues", residues);
            sequence.setAttribute("length", String.valueOf(residues.length()));
            sequence.setAttribute("md5checksum", md5checksum);
            sequences.add(sequence);
            String identifier = getFastaIdentifier(bioJavaSequence);
            // some mRNA FASTAs contain non-mRNAs, hopefully identified by their identifier
            if (identifier.contains("rRNA")) continue;
            if (identifier.contains("tRNA")) continue;
            if (identifier.contains("snRNA")) continue;
            if (identifier.contains("snoRNA")) continue;
            // HACK: don't allow spaces or tabs in primary identifiers; set symbol=extra part if present
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
            // MRNA Item
            Item mRNA = getMRNA(identifier);
            if (isPrimary) mRNA.setAttribute("isPrimary", "true");
            mRNA.setReference("sequence", sequence);
            mRNA.setAttribute("length", String.valueOf(residues.length()));
            if (symbol!=null) mRNA.setAttribute("symbol", symbol);
        }
    }
    
    /**
     * Process a gfa.tsv file which contains relationships between gene families, genes and proteins, along with a score value.
     *
     * 0     1      2    3    4    5           6    7   8
     * phalu.G27455.gnm1.ann1.JD7C.legfed_v1_0.M65K.gfa.tsv
     * 0                1
     * ScoreMeaning	e-value
     *
     * 0=gene                                 1=gene family        2=protein                                3=score
     * phalu.G27455.gnm1.ann1.tig000546640010 legfed_v1_0.L_00CL8T phalu.G27455.gnm1.ann1.tig000546640010.1 2.4e-68
     */
    void processGFAFile(Reader reader) throws IOException, RuntimeException {
        String[] fileParts = getCurrentFile().getName().split("\\.");
        if (fileParts.length!=9) {
            System.err.println("WARNING: GFA file does not have the required 9 dot-separated parts: "+getCurrentFile().getName());
        }
        String version = DEFAULT_VERSION;
        if (fileParts.length>5) version = fileParts[5];
        String scoreMeaning = null;
        // spin through the file
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        int linenumber = 0;
        while ((line=br.readLine())!=null) {
            linenumber++;
            if (line.startsWith("#")) continue;
            String[] fields = line.split("\t");
            if (fields.length==2) {
                // header
                if (fields[0].equals("ScoreMeaning")) scoreMeaning = fields[1];
            } else if (fields.length>2) {
                String geneIdentifier = fields[0];
                String geneFamilyIdentifier = fields[1];
                String proteinIdentifier = fields[2];
                double score = 0.0;
                boolean hasScore = false;
                if (fields.length>3) {
                    hasScore = true;
                    score = Double.parseDouble(fields[3]);
                }
                // validation
                if (geneIdentifier==null || geneIdentifier.trim().length()==0) {
                    throw new RuntimeException("ERROR: Gene.primaryIdentifier="+geneIdentifier+" at line "+linenumber);
                }
                if (proteinIdentifier==null || proteinIdentifier.trim().length()==0) {
                    throw new RuntimeException("ERROR: Protein.primaryIdentifier="+proteinIdentifier+" at line "+linenumber);
                }
                // Gene Family
                Item geneFamily = getGeneFamily(geneFamilyIdentifier);
                geneFamily.setAttribute("version", version);
                // Gene
                Item gene = getGene(geneIdentifier);
                gene.setReference("geneFamily", geneFamily);
                if (hasScore) {
                    if (scoreMeaning!=null) gene.setAttribute("geneFamilyScoreMeaning", scoreMeaning);
                    gene.setAttribute("geneFamilyScore", String.valueOf(score));
                }
                // Protein - this is just the primary transcript; there are other proteins associated with this gene that don't get the geneFamily reference
                Item protein = getProtein(proteinIdentifier);
                protein.setReference("geneFamily", geneFamily);
                if (hasScore) {
                    if (scoreMeaning!=null) protein.setAttribute("geneFamilyScoreMeaning", scoreMeaning);
                    protein.setAttribute("geneFamilyScore", String.valueOf(score));
                }
            }
        }
        br.close();
    }

    /**
     * Process a pathway.tsv file, associating genes with pathways.
     */
    void processPathwayFile(Reader reader) throws IOException {
        // spin through the file
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#")) continue;
            String[] fields = line.split("\t");
            if (fields.length==3) {
                String pathwayIdentifier = fields[0];
                String pathwayName = fields[1];
                String geneIdentifier = fields[2];
                Item gene = getGene(geneIdentifier);
                Item pathway = getPathway(pathwayIdentifier);
                pathway.setAttribute("name", pathwayName);
                gene.addToCollection("pathways", pathway);
            }
        }
    }


    /**
     * Process a GFF3 file, referenced by filename because GFFReader doesn't have a method to parse a Reader.
     * Assumes that ID=full-yuck-LIS-identifier and Name=name
     */
    void processGFF3File() throws IOException, RuntimeException {
        if (readme==null) {
            throw new RuntimeException("README not read before "+getCurrentFile().getName()+". Aborting.");
        }
        FeatureList featureList = GFF3Reader.read(getCurrentFile().getPath());
        for (FeatureI featureI : featureList) {
            String seqname = featureI.seqname();
            Location location = featureI.location();
            String type = featureI.type();
            // attributes
            String id = featureI.getAttribute("ID");
            String name = featureI.getAttribute("Name");
            String parent = featureI.getAttribute("Parent");
            String note = featureI.getAttribute("Note");
            String dbxref = featureI.getAttribute("Dbxref");
            String ontology_term = featureI.getAttribute("Ontology_term");
            String alleles = featureI.getAttribute("alleles");
            // check that id exists and matches collection
            if (id==null) {
                throw new RuntimeException("GFF line does not include ID: "+featureI.toString());
            }
            if (!matchesCollection(id)) {
                throw new RuntimeException("ID "+id+" does not match collection "+readme.identifier);
            }
            // get associated class
            String featureClass = featureClasses.get(type);
            if (featureClass==null) {
                throw new RuntimeException("GFF3 type "+type+" is not associated with a class in the featureClasses Map.");
            }
            // get the feature, e.g. ID=glyma.Lee.gnm1.ann1.GlymaLee.02G198600;
            // Gene and MRNA Items are stored separately
            Item feature = null;
            if (featureClass.equals("Gene")) {
                feature = getGene(id);
                placeFeatureOnSequence(feature, seqname, location);
                feature.setAttribute("length", String.valueOf(location.length()));
            } else if (featureClass.equals("MRNA")) {
                // don't set MRNA.length since that comes from FASTA
                feature = getMRNA(id);
                placeFeatureOnSequence(feature, seqname, location);
            } else {
                feature = getFeature(id, featureClass);
                placeFeatureOnSequence(feature, seqname, location);
                feature.setAttribute("length", String.valueOf(location.length()));
            }
            // Name=GlymaLee.02G198600;
            if (name!=null) {
                feature.setAttribute("name", name);
            }
            // Dbxref=Gene3D:G3DSA:1.10.630.10,InterPro:IPR001128,InterPro:IPR002401,InterPro:IPR017972,PANTHER:PTHR24298,...
            if (dbxref!=null) {
                String[] terms = dbxref.split(",");
                for (String term : terms) {
                    if (term.startsWith("InterPro:")) {
                        // InterPro are protein domains, not ontology annotations, and we only associate them with genes
                        if (type.equals("gene")) {
                            String[] parts = term.split(":");
                            String identifier = parts[1]; // IPR001128
                            Item proteinDomain = getProteinDomain(identifier);
                            feature.addToCollection("proteinDomains", proteinDomain);
                        }
                    } else {
                        // create/add ontology annotation for this feature
                        createOntologyAnnotation(feature, term);
                    }
                }
            }
            // Note=Cytochrome P450 superfamily protein%3B IPR001128 (Cytochrome P450)%3B GO:0005506 (iron ion binding)%2C GO:0020037 (heme binding)%2C ...
            if (note!=null) {
                feature.setAttribute("description", note);
            }
            // Parent is not required but must already be loaded if present
            if (parent!=null) {
                Item parentItem = features.get(parent);
                if (parentItem==null) parentItem = genes.get(parent);
                if (parentItem==null) parentItem = mRNAs.get(parent);
                if (parentItem==null) {
                    throw new RuntimeException("Parent "+parent+" not loaded before child "+id+". Is the GFF sorted?");
                }
                parentItem.addToCollection("childFeatures", feature);
            }
            // Ontology_term=GO:0005506,GO:0016705,GO:0020037,GO:0055114;
            if (ontology_term!=null) {
                String[] terms = ontology_term.split(",");
                for (String term : terms) {
                    createOntologyAnnotation(feature, term);
                }
            }
        }
    }

    /**
     * Add an OntologyAnnotation with the given identifier to the given feature's collection
     * NOTE: GO terms are GOTerm objects.
     */
    void createOntologyAnnotation(Item feature, String identifier) {
        Item ontologyTerm = ontologyTerms.get(identifier);
        if (ontologyTerm==null) {
            if (identifier.startsWith("GO:")) {
                ontologyTerm = createItem("GOTerm");
            } else {
                ontologyTerm = createItem("OntologyTerm");
            }
            ontologyTerm.setAttribute("identifier", identifier);
            ontologyTerms.put(identifier, ontologyTerm);
        }
        Item annotation = createItem("OntologyAnnotation");
        annotation.setReference("subject", feature);
        annotation.setReference("ontologyTerm", ontologyTerm);
        ontologyAnnotations.add(annotation);
    }

    /**
     * Add to a gene's proteinDomains collection if the 
     * Right now we only support InterPro terms (protein domains) associated with genes.
     */
    void addProteinDomain(Item gene, String term) {
        if (term.startsWith("InterPro:")) {
        }
    }

    /**
     * Get/add a feature of the given class in the features map.
     */
    Item getFeature(String primaryIdentifier, String className) {
        if (features.containsKey(primaryIdentifier)) {
            return features.get(primaryIdentifier);
        } else {
            Item feature = createItem(className);
            feature.setAttribute("primaryIdentifier", primaryIdentifier);
            feature.setAttribute("secondaryIdentifier", getSecondaryIdentifier(primaryIdentifier));
            features.put(primaryIdentifier, feature);
            return feature;
        }
    }

    /**
     * Place a feature on a sequence, determining whether it's a Chromosome or Supercontig from its name.
     */
    void placeFeatureOnSequence(Item feature, String seqname, Location location) {
        if (dsu.isSupercontig(seqname)) {
            Item supercontig = getSupercontig(seqname);
            // reference feature on supercontig
            feature.setReference("supercontig", supercontig);
            // reference feature on new IM Location
            Item supercontigLocation = createItem("Location");
            supercontigLocation.setReference("feature", feature);
            if (location.isNegative()) {
                supercontigLocation.setAttribute("strand", "-1");
            } else {
                supercontigLocation.setAttribute("strand", "1");
            }
            supercontigLocation.setAttribute("start", String.valueOf(location.bioStart()));
            supercontigLocation.setAttribute("end", String.valueOf(location.bioEnd()));
            supercontigLocation.setReference("locatedOn", supercontig);
            locations.add(supercontigLocation);
            feature.setReference("supercontigLocation", supercontigLocation);
        } else {
            Item chromosome = getChromosome(seqname);
            // reference feature on chromosome
            feature.setReference("chromosome", chromosome);
            // reference feature on new IM Location
            Item chromosomeLocation = createItem("Location");
            chromosomeLocation.setReference("feature", feature);
            if (location.isNegative()) {
                chromosomeLocation.setAttribute("strand", "-1");
            } else {
                chromosomeLocation.setAttribute("strand", "1");
            }
            chromosomeLocation.setAttribute("start", String.valueOf(location.bioStart()));
            chromosomeLocation.setAttribute("end", String.valueOf(location.bioEnd()));
            chromosomeLocation.setReference("locatedOn", chromosome);
            locations.add(chromosomeLocation);
            feature.setReference("chromosomeLocation", chromosomeLocation);
        }
    }

    /**
     * Get/add a Gene Item, keyed by primaryIdentifier
     */
    Item getGene(String primaryIdentifier) {
        if (genes.containsKey(primaryIdentifier)) {
            return genes.get(primaryIdentifier);
        } else {
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", primaryIdentifier);
            gene.setAttribute("secondaryIdentifier", getSecondaryIdentifier(primaryIdentifier));
            genes.put(primaryIdentifier, gene);
            return gene;
        }
    }

    /**
     * Get/add a Pathway Item, keyed by primaryIdentifier
     */
    Item getPathway(String identifier) {
        if (pathways.containsKey(identifier)) {
            return pathways.get(identifier);
        } else {
            Item pathway = createItem("Pathway");
            pathway.setAttribute("identifier", identifier);
            pathways.put(identifier, pathway);
            return pathway;
        }
    }

    /**
     * Extract the secondaryIdentifier from an annotation feature primaryIdentifier.
     * Do not use this for chromosomes or supercontigs.
     */
    String getSecondaryIdentifier(String primaryIdentifier) throws RuntimeException {
        String secondaryIdentifier = null;
        secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, true);
        if (secondaryIdentifier==null) {
            throw new RuntimeException("Could not get secondaryIdentifier for feature:"+primaryIdentifier+" in file:"+getCurrentFile().getName());
        }
        return secondaryIdentifier;
    }

    /**
     * Get/add a Chromosome Item, keyed by primaryIdentifier with secondaryIdentifier from primaryIdentifier.
     */
    Item getChromosome(String primaryIdentifier) throws RuntimeException {
        if (chromosomes.containsKey(primaryIdentifier)) {
            return chromosomes.get(primaryIdentifier);
        } else {
            Item chromosome = createItem("Chromosome");
            chromosome.setAttribute("primaryIdentifier", primaryIdentifier);
            String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, false);
            if (secondaryIdentifier==null) {
                throw new RuntimeException("Could not get secondaryIdentifier for chromosome:"+primaryIdentifier);
            }
            chromosomes.put(primaryIdentifier, chromosome);
            return chromosome;
        }
    }

    /**
     * Get/add a Supercontig Item, keyed by primaryIdentifier with secondaryIdentifier from primaryIdentifier.
     */
    Item getSupercontig(String primaryIdentifier) throws RuntimeException {
        if (supercontigs.containsKey(primaryIdentifier)) {
            return supercontigs.get(primaryIdentifier);
        } else {
            Item supercontig = createItem("Supercontig");
            supercontig.setAttribute("primaryIdentifier", primaryIdentifier);
            String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, false);
            if (secondaryIdentifier==null) {
                throw new RuntimeException("Could not get secondaryIdentifier for supercontig:"+primaryIdentifier);
            }
            supercontigs.put(primaryIdentifier, supercontig);
            return supercontig;
        }
    }

    /**
     * Get/add a Protein Item, keyed by primaryIdentifier, with secondaryIdentifier from primaryIdentifier.
     */
    Item getProtein(String primaryIdentifier) {
        if (proteins.containsKey(primaryIdentifier)) {
            return proteins.get(primaryIdentifier);
        } else {
            Item protein = createItem("Protein");
            protein.setAttribute("primaryIdentifier", primaryIdentifier);
            protein.setAttribute("secondaryIdentifier", getSecondaryIdentifier(primaryIdentifier));
            proteins.put(primaryIdentifier, protein);
            return protein;
        }
    }

    /**
     * Get/add a ProteinDomain item, keyed by primaryIdentifier.
     */
    Item getProteinDomain(String primaryIdentifier) {
        if (proteinDomains.containsKey(primaryIdentifier)) {
            return proteinDomains.get(primaryIdentifier);
        } else {
            Item proteinDomain = createItem("ProteinDomain");
            proteinDomain.setAttribute("primaryIdentifier", primaryIdentifier);
            proteinDomains.put(primaryIdentifier, proteinDomain);
            return proteinDomain;
        }
    }
    
    /**
     * Get/add a CDS Item, keyed by primaryIdentifier, with secondaryIdentifier from primaryIdentifier
     */
    Item getCDS(String primaryIdentifier) {
        if (cdses.containsKey(primaryIdentifier)) {
            return cdses.get(primaryIdentifier);
        } else {
            Item cds = createItem("CDS");
            cds.setAttribute("primaryIdentifier", primaryIdentifier);
            cds.setAttribute("secondaryIdentifier", getSecondaryIdentifier(primaryIdentifier));
            cdses.put(primaryIdentifier, cds);
            return cds;
        }
    }
    
    /**
     * Get/add a MRNA Item, keyed by primaryIdentifier, with secondaryIdentifier from primaryIdentifier
     */
    Item getMRNA(String primaryIdentifier) {
        if (mRNAs.containsKey(primaryIdentifier)) {
            return mRNAs.get(primaryIdentifier);
        } else {
            Item mRNA = createItem("MRNA");
            mRNA.setAttribute("primaryIdentifier", primaryIdentifier);
            mRNA.setAttribute("secondaryIdentifier", getSecondaryIdentifier(primaryIdentifier));
            mRNAs.put(primaryIdentifier, mRNA);
            return mRNA;
        }
    }
    
    /**
     * Get/add a GeneFamily, keyed by identifier
     */
    Item getGeneFamily(String identifier) {
        if (geneFamilies.containsKey(identifier)) {
            return geneFamilies.get(identifier);
        } else {
            Item geneFamily = createItem("GeneFamily");
            geneFamily.setAttribute("identifier", identifier);
            geneFamilies.put(identifier, geneFamily);
            return geneFamily;
        }
    }

    /**
     * Get/add an OntologyTerm Item, associated with the given ontology and keyed by identifier.
     *
     * @param identifier the ontology term identifier
     * @param ontology the associated ontology Item
     * @return the OntologyTerm item
     */
    Item getOntologyTerm(String identifier, Item ontology) {
        if (ontologyTerms.containsKey(identifier)) {
            return ontologyTerms.get(identifier);
        } else {
            Item ontologyTerm = createItem("OntologyTerm");
            ontologyTerm.setAttribute("identifier", identifier);
            ontologyTerm.setReference("ontology", ontology);
            ontologyTerms.put(identifier, ontologyTerm);
            return ontologyTerm;
        }
    }

    /**
     * For the given BioJava Sequence object, return an identifier to be used when creating the corresponding BioEntity.
     *
     * @param bioJavaSequence the Sequence
     * @return an identifier
     */
    protected String getFastaIdentifier(AbstractSequence bioJavaSequence) {
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
    static String getFastaAttribute(AbstractSequence bioJavaSequence, String field) {
        String header = bioJavaSequence.getAccession().getID();
        if (!header.contains(field)) return null;
        String[] split = header.split(" "+field+"=");
        if (split.length==1) return null;
        String secondHalf = split[1];
        String[] parts = secondHalf.split(" ");
        return parts[0];
    }

    /**
     * Return a gene ID as the portion of an ID with the .1, .2, etc. removed
     *
     * @param dotId the dot.number ID
     * @return the corresponding gene ID
     */
    static String getGeneIdFromDotId(String dotId) {
        String[] parts = dotId.split("\\.");
        String geneId = parts[0];
        for (int i=1; i<parts.length-1; i++) geneId += "."+parts[i];
        return geneId;
    }
}
