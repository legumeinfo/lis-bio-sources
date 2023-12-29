package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
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
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.genome.parsers.gff.Location;

import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.Util;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.ncgr.datastore.Readme;
import org.ncgr.datastore.validation.AnnotationCollectionValidator;
import org.ncgr.zip.GZIPFastaReader;
import org.ncgr.zip.GZIPBufferedReader;

/**
 * Loads data from an LIS datastore annotations collection.
 * Files types loaded are:
 *   README.yaml
 *   gene_models_main.gff3.gz
 *   protein.faa.gz
 *   protein_primary.faa.gz
 *   cds.fna.gz
 *   cds_primary.fna.gz
 *   mrna.fna.gz
 *   mrna_primary.fna.gz
 *   gfa.tsv.gz
 *   iprscan.gff3.gz
 *
 * @author Sam Hokin
 */
public class AnnotationFileConverter extends DatastoreFileConverter {

    // spit out debug lines if not null and ID starts with this
    private static final String DEBUG_ID = null;

    private static final Logger LOG = Logger.getLogger(AnnotationFileConverter.class);
    private static final String TEMPGENEFILE = "/tmp/gene_models_main.gff3";
    private static final String TEMPIPRSCANFILE = "/tmp/iprscan.gff3";

    // GFF sourced
    Map<String,Item> chromosomes = new HashMap<>();
    Map<String,Item> supercontigs = new HashMap<>();
    Map<String,Item> features = new HashMap<>();       // features other than genes and mRNAs
    Map<String,Item> genes = new HashMap<>();          // also GFA
    Map<String,Item> mRNAs = new HashMap<>();          // mRNA sequences,length are FASTA-sourced
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> proteinDomains = new HashMap<>();
    Map<String,Item> featureLocations = new HashMap<>(); // temp storage keyed by feature.primaryIdentifier
    List<Item> ontologyAnnotations = new ArrayList<>();
    List<Item> locations = new ArrayList<>();
    Map<String,Item> cdses = new HashMap<>();              // also FASTA sourced
    Map<String,Item> cdsContigLocations = new HashMap<>(); // store the CDS chromosomeLocation/supercontigLocation for updating

    // TSV sourced
    Map<String,Item> geneFamilies = new HashMap<>();

    // GFA sourced
    List<Item> geneFamilyAssignments = new ArrayList<>();

    // FASTA sourced
    List<Item> sequences = new ArrayList<>();
    Map<String,Item> proteins = new HashMap<>(); // also GFA, etc.

    // validate the collection first and only once by storing a flag
    boolean collectionValidated = false;

    // map GFF types to InterMine classes; be sure to include extras in the additions file!
    Map<String,String> featureClasses = Map.ofEntries(entry("gene", "Gene"),
                                                      entry("mRNA", "MRNA"),
                                                      entry("CDS", "CDS"),
                                                      entry("exon", "Exon"),
                                                      entry("intron", "Intron"),
                                                      entry("three_prime_UTR", "ThreePrimeUTR"),
                                                      entry("five_prime_UTR", "FivePrimeUTR"),
                                                      entry("lnc_RNA", "LncRNA"),
                                                      entry("transcript", "Transcript"),
                                                      entry("pseudogene", "Pseudogene"),
                                                      entry("primary_transcript", "Transcript"),
                                                      entry("miRNA", "MiRNA"),
                                                      entry("miRNA_primary_transcript", "MiRNA"),
                                                      entry("ncRNA", "NcRNA"),
                                                      entry("pre_miRNA", "PreMiRNA"),
                                                      entry("transposable_element_gene", "Gene"),
                                                      entry("tRNA", "TRNA"),
                                                      entry("tRNA_primary_transcript", "TRNA"),
                                                      entry("pseudogenic_tRNA", "TRNA"),
                                                      entry("snoRNA", "SnoRNA"),
                                                      entry("snRNA", "SnRNA"),
                                                      entry("rRNA", "RRNA"),
                                                      entry("rRNA_primary_transcript", "RRNA"),
                                                      entry("inverted_repeat", "InvertedRepeat"),
                                                      entry("region", "Region"),
                                                      entry("repeat_region", "RepeatRegion"));

    /**
     * Create a new AnnotationFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public AnnotationFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (!collectionValidated) {
            AnnotationCollectionValidator validator = new AnnotationCollectionValidator(getCurrentFile().getParent());
            validator.validate();
            if (!validator.isValid()) {
                throw new RuntimeException("Collection " + getCurrentFile().getParent() + " does not pass validation.");
            }
            collectionValidated = true;
        }
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
            setStrain();
            processGenomeReadme(getCurrentFile());
        } else if (getCurrentFile().getName().endsWith(".gene_models_main.gff3.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processGeneModelsMainGFF3File();
        } else if (getCurrentFile().getName().endsWith(".gfa.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processGFAFile();
	} else if (getCurrentFile().getName().endsWith(".protein.faa.gz") || getCurrentFile().getName().endsWith(".protein_primary.faa.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processProteinFasta();
        } else if (getCurrentFile().getName().endsWith(".cds.fna.gz") || getCurrentFile().getName().endsWith(".cds_primary.fna.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processCDSFasta();
        } else if (getCurrentFile().getName().endsWith(".mrna.fna.gz") || getCurrentFile().getName().endsWith(".mrna_primary.fna.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processMRNAFasta();
        } else if (getCurrentFile().getName().endsWith(".iprscan.gff3.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processIPRScanGFF3();
        } else if (getCurrentFile().getName().endsWith(".gz")) {
            System.out.println(" x skipping "+getCurrentFile().getName());
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException, RuntimeException {
        // set missing CDS.transcript reference and add to MRNA.CDSs collection based on matching identifiers
        // (doing this here allows for CDS in GFF to have different parent mRNA identifier)
        for (String primaryIdentifier : cdses.keySet()) {
            Item cds = cdses.get(primaryIdentifier);
            if (cds.getReference("transcript").getRefId() == null && mRNAs.containsKey(primaryIdentifier)) {
                Item mRNA = mRNAs.get(primaryIdentifier);
                cds.setReference("transcript", mRNA);
                mRNA.addToCollection("CDSs", cds);
            }
            // if (proteins.containsKey(primaryIdentifier)) cds.setReference("protein", proteins.get(primaryIdentifier));
        }

        // add publication to all Annotatables
        if (publication != null) {
            for (Item chromosome : chromosomes.values()) {
                chromosome.addToCollection("publications", publication);
            }
            for (Item supercontig : supercontigs.values()) {
                supercontig.addToCollection("publications", publication);
            }
            for (Item feature : features.values()) {
                feature.addToCollection("publications", publication);
            }
            for (Item cds : cdses.values()) {
                cds.addToCollection("publications", publication);
            }
            for (Item gene : genes.values()) {
                gene.addToCollection("publications", publication);
            }
            for (Item mRNA : mRNAs.values()) {
                mRNA.addToCollection("publications", publication);
            }
            for (Item protein : proteins.values()) {
                protein.addToCollection("publications", publication);
            }
            for (Item proteinDomain : proteinDomains.values()) {
                proteinDomain.addToCollection("publications", publication);
            }
        }

        // store our Items
        storeCollectionItems();
        store(chromosomes.values());
        store(supercontigs.values());
        store(features.values());
        store(cdses.values());
        store(cdsContigLocations.values());
        store(genes.values());
        store(mRNAs.values());
        store(proteins.values());
        store(geneFamilies.values());
        store(ontologyAnnotations);
        store(ontologyTerms.values());
        store(proteinDomains.values());
        store(locations);
        store(sequences);
        store(geneFamilyAssignments);
    }

    /**
     * Process a protein FASTA faa file.
     */
    void processProteinFasta() throws IOException {
        LinkedHashMap<String,ProteinSequence> sequenceMap = GZIPFastaReader.readFastaProteinSequence(getCurrentFile());
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
            // we don't parse a name from the FASTA header, so...
            protein.setAttribute("name", getSecondaryIdentifier(identifier));
            if (isPrimary()) protein.setAttribute("isPrimary", "true");
            protein.setReference("sequence", sequence);
            protein.setAttribute("length", String.valueOf(residues.length()));
            protein.setAttribute("md5checksum", md5checksum);
            if (symbol != null) protein.setAttribute("symbol", symbol);
            // associated transcript
            // Item mRNA = getMRNA(identifier);
            // protein.setReference("transcript", mRNA);
            // associated CDS
            // Item cds = getCDS(identifier);
            // protein.setReference("CDS", cds);
        }
    }

    /**
     * Process a CDS nucleotide FASTA fna file
     */
    void processCDSFasta() throws IOException {
        LinkedHashMap<String,DNASequence> sequenceMap = GZIPFastaReader.readFastaDNASequence(getCurrentFile());
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
            if (isPrimary()) cds.setAttribute("isPrimary", "true");
            cds.setReference("sequence", sequence);
            cds.setAttribute("length", String.valueOf(residues.length()));
            if (symbol != null) cds.setAttribute("symbol", symbol);
        }
    }

    /**
     * Process an mRNA nucleotide FASTA fna file
     */
    void processMRNAFasta() throws IOException {
        LinkedHashMap<String,DNASequence> sequenceMap = GZIPFastaReader.readFastaDNASequence(getCurrentFile());
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
            if (isPrimary()) mRNA.setAttribute("isPrimary", "true");
            mRNA.setReference("sequence", sequence);
            mRNA.setAttribute("length", String.valueOf(residues.length()));
            if (symbol != null) mRNA.setAttribute("symbol", symbol);
        }
    }
    
    /**
     * Process a gfa.tsv file which contains relationships between gene families, genes and proteins, along with an e-value, score, and best-domain score.
     * 0     1      2    3    4    5           6    7   8   9
     * phavu.G19833.gnm2.ann1.PB8d.legfed_v1_0.M65K.gfa.tsv.gz
     *
     * gene                                    family                  protein                                   e-value   score  best-domain-score
     * 0                                       1                       2                                         3         4      5
     * vigun.CB5-2.gnm1.ann1.VuCB5-2.03G238200 legfed_v1_0.L_JSS52Z    vigun.CB5-2.gnm1.ann1.VuCB5-2.03G238200.1 5.3e-199  660.9  660.8
     */
    void processGFAFile() throws IOException, RuntimeException {
        String[] fileParts = getCurrentFile().getName().split("\\.");
        if (fileParts.length!=10) {
            System.err.println("WARNING: GFA file name does not have the required 10 dot-separated parts: "+getCurrentFile().getName());
        }
        // spin through the file
        BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        int linenumber = 0;
        while ((line=br.readLine()) != null) {
            linenumber++;
            if (line.startsWith("#")) continue;
            String[] fields = line.split("\t");
            if (fields.length==2) continue; // probably a second header with some item of non-interest
            String geneIdentifier = fields[0];
            String geneFamilyIdentifier = fields[1];
            String proteinIdentifier = fields[2];
            // data columns are optional
            double evalue = 0.0;
            double score = 0.0;
            double bestDomainScore = 0.0;
            if (fields.length>3) evalue = Double.parseDouble(fields[3]);
            if (fields.length>4) score = Double.parseDouble(fields[4]);
            if (fields.length>5) bestDomainScore = Double.parseDouble(fields[5]);
            // Gene Family
            Item geneFamily = getGeneFamily(geneFamilyIdentifier);
            // Gene
            Item gene = getGene(geneIdentifier);
            // Protein - this is just the primary transcript; there are other proteins associated with this gene that don't get the geneFamily reference
            Item protein = getProtein(proteinIdentifier);
            // GeneFamily collections
            geneFamily.addToCollection("genes", gene);
            geneFamily.addToCollection("proteins", protein);
            // scores go into GeneFamilyAssignment along with gene and protein
            createGeneFamilyAssignment(geneFamily, gene, protein, evalue, score, bestDomainScore);
        }
    }

    /**
     * Process a gzipped GFF3 file.
     * Because GFFReader does not handle gzipped files, we have to dump the file to /tmp first and then read that.
     * Assumes that ID=gensp.strain.gnm.ann.identifier and Name=name.
     */
    void processGeneModelsMainGFF3File() throws IOException, RuntimeException {
        if (readme == null) {
            throw new RuntimeException("README not read before "+getCurrentFile().getName()+". Aborting.");
        }
        // uncompress the gff3.gz file to a temp file
        File tempfile = new File(TEMPGENEFILE);
        tempfile.delete();
        BufferedWriter writer = new BufferedWriter(new FileWriter(tempfile));
        BufferedReader reader = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        while ( (line=reader.readLine()) != null ) {
            writer.write(line);
            writer.newLine();
        }
        writer.close();
        // now load the uncompressed GFF
        FeatureList featureList = GFF3Reader.read(TEMPGENEFILE);
        for (FeatureI featureI : featureList) {
            String seqname = featureI.seqname();
            // check that we're on a recognized Chromosome or Supercontig, else bail with warning
            if (!isChromosome(seqname) && !isSupercontig(seqname)) {
                System.out.println("### feature ignored on sequence " + seqname + " because not recognized as chromosome or supercontig.");
                continue;
            }
            Location location = featureI.location();
            String type = featureI.type();
            // GFF spec attributes
            String id = getAttribute(featureI, "ID");
            String name = getAttribute(featureI, "Name"); // if missing, store derived secondaryIdentifier below
            String alias = getAttribute(featureI, "Alias");
            String parent = getAttribute(featureI, "Parent");
            String target = getAttribute(featureI, "Target");
            String gap = getAttribute(featureI, "Gap");
            String derivesFrom = getAttribute(featureI, "Derives_from");
            String note = getAttribute(featureI, "Note");
            String dbxref = getAttribute(featureI, "Dbxref");
            String ontology_term = getAttribute(featureI, "Ontology_term");
            String isCircular = getAttribute(featureI, "Is_circular");
            // LIS-specific attributes
            String alleles = getAttribute(featureI, "alleles");
            String symbol = getAttribute(featureI, "symbol");
            // check that id exists and matches collection
            if (id == null) {
                throw new RuntimeException("GFF record does not include ID: "+featureI.toString());
            }
            if (!matchesCollection(id)) {
                throw new RuntimeException("ID " + id + " format does not match collection " + readme.identifier);
            }
            // get associated class
            String featureClass = featureClasses.get(type);
            if (featureClass == null) {
                throw new RuntimeException("GFF3 type " + type + " is not associated with a class in the featureClasses Map.");
            }
            // optional DEBUG
            if (DEBUG_ID != null && id.startsWith(DEBUG_ID)) printGFFDebugOutput(featureI);
            // get the feature, e.g. ID=glyma.Lee.gnm1.ann1.GlymaLee.02G198600;
            // Gene and MRNA Items are stored separately for specific treatment
            Item feature = null;
            if (featureClass.equals("Gene")) {
                feature = getGene(id);
                placeFeatureOnContig(feature, seqname, location);
                feature.setAttribute("length", String.valueOf(location.length()));
                // set Gene.ensemblName, if available (else null)
                if (name != null) {
                    String ensemblName = DatastoreUtils.getEnsemblName(name);
                    if (ensemblName != null) feature.setAttribute("ensemblName", ensemblName);
                }
            } else if (featureClass.equals("CDS")) {
                // CDS record must have a parent attribute for CDS.transcript reference
                if (parent == null) {
                    throw new RuntimeException("CDS " + id + " lacks Parent transcript attribute.");
                }
                if (!mRNAs.containsKey(parent)) {
                    throw new RuntimeException("CDS " + id + " parent mRNA " + parent + " has not yet been loaded. Is the GFF sorted?");
                }
                // use the id of the parent for the CDS - this should match the id in the CDS FASTAs
                feature = getCDS(parent);
                // add this location to the CDS.locations collection and update chromosome/supercontigLocation
                placeCDSOnContig(parent, seqname, location);
                // don't set CDS.length since that comes from FASTA or stitching together CDS.locations
                // set the CDS.transcript reference to the parent mRNA
                Item mRNA = mRNAs.get(parent);
                feature.setReference("transcript", mRNA);
                // set name = mRNA.name
                if (mRNA.getAttribute("name") != null) {
                    name = mRNA.getAttribute("name").getValue();
                }
            } else if (featureClass.equals("Exon")) {
                // exon record must have a parent attribute for Exon.transcripts collection
                if (parent == null) {
                    throw new RuntimeException("Exon " + id + " lacks Parent attribute.");
                }
                if (!mRNAs.containsKey(parent)) {
                    throw new RuntimeException("Exon " + id + " parent mRNA " + parent + " <has not yet been loaded. Is the GFF sorted?");
                }
                feature = getFeature(id, featureClass);
                placeFeatureOnContig(feature, seqname, location);
                feature.setAttribute("length", String.valueOf(location.length()));
                // add the parent mRNA to the Exon.transcripts collection
                Item mRNA = mRNAs.get(parent);
                feature.addToCollection("transcripts", mRNA);
            } else if (featureClass.equals("MRNA")) {
                // mRNA record must have a parent attribute for Transcript.gene reference
                if (parent == null) {
                    throw new RuntimeException("mRNA " + id + " lacks Parent attribute.");
                }
                if (!genes.containsKey(parent)) {
                    throw new RuntimeException("mRNA " + id + " parent gene " + parent + " has not yet been loaded. Is the GFF sorted?");
                }
                feature = getMRNA(id);
                placeFeatureOnContig(feature, seqname, location);
                // don't set MRNA.length since that comes from FASTA or stitching together exons
                // set the MRNA.gene reference to the parent Gene
                Item gene = genes.get(parent);
                feature.setReference("gene", gene);
                // set the MRNA.protein reference to the identifier-matching protein (which may not yet be loaded)
                Item protein = getProtein(id);
                feature.setReference("protein", protein);
                protein.setReference("transcript", feature);
                // since we have the gene and the protein, add it to the protein.genes collection
                protein.addToCollection("genes", gene);
            } else if (featureClass.equals("ThreePrimeUTR") || featureClass.equals("FivePrimeUTR")) {
                // UTR record must have a parent attribute for UTR.transcripts collection
                if (parent == null) {
                    throw new RuntimeException("UTR " + id + " lacks Parent attribute.");
                }
                if (!mRNAs.containsKey(parent)) {
                    throw new RuntimeException("UTR " + id + " parent mRNA " + parent + " has not yet been loaded. Is the GFF sorted?");
                }
                feature = getFeature(id, featureClass);
                placeFeatureOnContig(feature, seqname, location);
                feature.setAttribute("length", String.valueOf(location.length()));
                // add the parent mRNA to the UTR.transcripts collection
                Item mRNA = mRNAs.get(parent);
                feature.addToCollection("transcripts", mRNA);
                // add the UTR to the MRNA.UTRs collection
                mRNA.addToCollection("UTRs", feature);
            } else {
                // just a plain vanilla feature
                feature = getFeature(id, featureClass);
                placeFeatureOnContig(feature, seqname, location);
                feature.setAttribute("length", String.valueOf(location.length()));
            }
            // only Region has isCircular attribute
            if (featureClass.equals("Region") && isCircular != null) {
                if (isCircular.equals("true")) {
                    feature.setAttribute("isCircular", "true");
                } else {
                    feature.setAttribute("isCircular", "false");
                }
            }
            // Name=GlymaLee.02G198600; default to secondaryIdentifier formed from ID
            if (name != null && name.trim().length() > 0) {
                feature.setAttribute("name", name);
            } else {
                feature.setAttribute("name", feature.getAttribute("secondaryIdentifier").getValue());
            }
            // Note=Cytochrome P450 superfamily protein%3B IPR001128 (Cytochrome P450)%3B GO:0005506 (iron ion binding)%2C GO:0020037 (heme binding)%2C ...
            if (note != null && note.trim().length() > 0) feature.setAttribute("description", DatastoreUtils.unescape(note));
            // Symbol=RGB4
            if (symbol != null && symbol.trim().length() > 0) feature.setAttribute("symbol", symbol);
            // Dbxref=Gene3D:G3DSA:1.10.630.10,InterPro:IPR001128,InterPro:IPR002401,InterPro:IPR017972,PANTHER:PTHR24298,...
            if (dbxref != null && dbxref.trim().length() > 0) {
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

            // add to childFeatures collection of parent
            // parent must be loaded before this feature
            if (parent != null) {
                List<String> parents = new ArrayList<>();
                if (parent.contains(",")) {
                    for (String p : parent.split(",")) parents.add(p);
                } else {
                    parents.add(parent);
                }
                for (String p : parents) {
                    Item parentItem = null;
                    if (features.containsKey(p)) {
                        // parent is not Gene or MRNA
                        parentItem = features.get(p);
                    } else if (genes.containsKey(p)) {
                        // parent is a gene
                        parentItem = genes.get(p);
                    } else if (mRNAs.containsKey(p)) {
                        // parent is an mRNA
                        parentItem = mRNAs.get(p);
                    } else {
                        throw new RuntimeException("Parent " + p + " not loaded before child " + id + ". Is the GFF sorted?");
                    }
                    // add this feature to the parent's children
                    parentItem.addToCollection("childFeatures", feature);
                }
            }

            // create ontology annotations for this feature from ontology_term attribute
            if (ontology_term != null && ontology_term.trim().length() > 0) {
                String[] terms = ontology_term.split(",");
                for (String term : terms) {
                    createOntologyAnnotation(feature, term);
                }
            }
        }
    }

    /**
     * Process an IPRScan GFF file which has Proteins in the sequence column.
     * 0     1            2    3    4    5       6    7
     * medsa.XinJiangDaYe.gnm1.ann1.RKB9.iprscan.gff3.gz
     *
     * medsa.XinJiangDaYe.gnm1.ann1.MS_gene000000.t1 PANTHER protein_match 1 463 . + . Name=PTHR10178:SF14;status=T;ID=match$1047424_1_463;date=04-02-2021
     * medsa.XinJiangDaYe.gnm1.ann1.MS_gene000008.t1 ProSiteProfiles protein_match  10 264   . + . Name=PS50294;status=T;ID=match$461640_10_264;date=03-02-2021;
     *    signature_desc=Trp-Asp (WD) repeats circular profile.
     * medsa.XinJiangDaYe.gnm1.ann1.MS_gene000000.t1 Pfam    protein_hmm_match 3 88 4.8E-11 + . Name=PF03732;status=T;ID=match$1047423_3_88;date=04-02-2021;
     *    signature_desc=Retrotransposon gag protein;Target=PF03732 6 96;
     *
     * NOTE: protein_match and protein_hmm_match IDs are non-unique between files for different species,
     * therefore we prefix for primaryIdentifier and store plain ID as secondaryIdentifier.
     *
     */
    void processIPRScanGFF3() throws IOException {
        // get prefix for uniquefying ProteinMatch and ProteinHmmMatch primaryIdentifier
        String prefix = DatastoreUtils.extractPrefixFromAnnotationFilename(getCurrentFile().getName());
        // uncompress the gff3.gz file to a temp file
        File tempfile = new File(TEMPIPRSCANFILE);
        tempfile.delete();
        BufferedWriter writer = new BufferedWriter(new FileWriter(tempfile));
        BufferedReader reader = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        while ( (line=reader.readLine()) != null ) {
            writer.write(line);
            writer.newLine();
        }
        writer.close();
        // now load the uncompressed GFF
        FeatureList featureList = GFF3Reader.read(TEMPIPRSCANFILE);
        for (FeatureI featureI : featureList) {
            String seqname = featureI.seqname();
            Location location = featureI.location();
            String type = featureI.type();
            // feature
            Item feature = getFeatureOnProtein(type, seqname, location);
            String id = getAttribute(featureI, "ID");
            if (type.equals("protein_match") || type.equals("protein_hmm_match")) {
                // prefix for unique primaryIdentifier
                feature.setAttribute("name", id);
                id = prefix + "." + id;
            }
            feature.setAttribute("primaryIdentifier", id);
            features.put(id, feature);
            // source isn't supplied by FeatureI
            // feature.setAttribute("source", rec.getSource());
            // attributes
            String name = getAttribute(featureI, "Name");
            String status = getAttribute(featureI, "status");
            String date = getAttribute(featureI, "date");
            String target = getAttribute(featureI, "Target");
            String signatureDesc = getAttribute(featureI, "signature_desc");
            // accession=Name
            feature.setAttribute("accession", name);
            // status
            if (status != null) feature.setAttribute("status", status);
            // date
            if (date != null) feature.setAttribute("date", date);
            // target
            if (target != null) feature.setAttribute("target", target);
            // signatureDesc
            if (signatureDesc != null) feature.setAttribute("signatureDesc", signatureDesc);
        }
    }

    /**
     * Add an OntologyAnnotation with the given identifier to the given feature's collection
     * NOTE: GO terms are GOTerm objects.
     */
    void createOntologyAnnotation(Item feature, String identifier) {
        Item ontologyTerm = ontologyTerms.get(identifier);
        if (ontologyTerm == null) {
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
            features.put(primaryIdentifier, feature);
            feature.setAttribute("primaryIdentifier", primaryIdentifier);
            feature.setAttribute("secondaryIdentifier", getSecondaryIdentifier(primaryIdentifier));
            feature.setAttribute("assemblyVersion", assemblyVersion);
            feature.setAttribute("annotationVersion", annotationVersion);
            feature.setReference("organism", organism);
            feature.setReference("strain", strain);
            return feature;
        }
    }

    /**
     * Place a feature on a contig from its GFF Location, determining whether the contig is a Chromosome or Supercontig from its name.
     */
    void placeFeatureOnContig(Item feature, String seqname, Location location) throws RuntimeException {
        // chromosome or supercontig?
        boolean onChromosome = isChromosome(seqname);
        boolean onSupercontig = isSupercontig(seqname);
        if (!onChromosome && !onSupercontig) {
            throw new RuntimeException("Attempted to place feature on sequence " + seqname + " not recognized as a Chromosome or Supercontig.");
        }
        Item contig = null;
        if (onChromosome) {
            contig = getChromosome(seqname);
            feature.setReference("chromosome", contig);
        } else {
            contig = getSupercontig(seqname);
            feature.setReference("supercontig", contig);
        }
        Item contigLocation = createItem("Location");
        contigLocation.setReference("feature", feature);
        contigLocation.setReference("locatedOn", contig);
        if (location.isNegative()) {
            contigLocation.setAttribute("strand", "-1");
        } else {
            contigLocation.setAttribute("strand", "1");
        }
        contigLocation.setAttribute("start", String.valueOf(location.bioStart()));
        contigLocation.setAttribute("end", String.valueOf(location.bioEnd()));
        locations.add(contigLocation);
        if (onChromosome) {
            feature.setReference("chromosomeLocation", contigLocation);
        } else {
            feature.setReference("supercontigLocation", contigLocation);
        }
    }

    /**
     * Place a CDS on a contig, extending its range on the contig if it's not the first one, and adding to its locations collection.
     * Note: getCDS should have been called before this.
     */
    void placeCDSOnContig(String primaryIdentifier, String seqname, Location location) throws RuntimeException {
        Item cds = cdses.get(primaryIdentifier); // should NOT be null
        // chromosome or supercontig?
        boolean onChromosome = isChromosome(seqname);
        boolean onSupercontig = isSupercontig(seqname);
        if (!onChromosome && !onSupercontig) {
            throw new RuntimeException("Attempted to place CDS on sequence " + seqname + " not recognized as a Chromosome or Supercontig.");
        }
        // get the contig
        Item contig = null;
        if (onChromosome) {
            contig = getChromosome(seqname);
            cds.setReference("chromosome", contig);
        } else {
            contig = getSupercontig(seqname);
            cds.setReference("supercontig", contig);
        }
        // create a new location for CDS.locations collection
        Item thisLocation = createItem("Location");
        thisLocation.setReference("feature", cds);
        thisLocation.setReference("locatedOn", contig);
        if (location.isNegative()) {
            thisLocation.setAttribute("strand", "-1");
        } else {
            thisLocation.setAttribute("strand", "1");
        }
        thisLocation.setAttribute("start", String.valueOf(location.bioStart()));
        thisLocation.setAttribute("end", String.valueOf(location.bioEnd()));
        cds.addToCollection("locations", thisLocation);
        locations.add(thisLocation);
        // update/create the chromosomeLocation/supercontigLocation
        if (cdsContigLocations.containsKey(primaryIdentifier)) {
            // get from storage map
            Item contigLocation = cdsContigLocations.get(primaryIdentifier);
            if (location.isNegative()) {
                contigLocation.setAttribute("strand", "-1");
            } else {
                contigLocation.setAttribute("strand", "1");
            }
            int start = Integer.parseInt(contigLocation.getAttribute("start").getValue());
            int end = Integer.parseInt(contigLocation.getAttribute("end").getValue());
            if (location.bioStart() < start) {
                contigLocation.setAttribute("start", String.valueOf(location.bioStart()));
            }
            if (location.bioEnd() > end) {
                contigLocation.setAttribute("end", String.valueOf(location.bioEnd()));
            }
        } else {
            Item contigLocation = createItem("Location");
            if (location.isNegative()) {
                contigLocation.setAttribute("strand", "-1");
            } else {
                contigLocation.setAttribute("strand", "1");
            }
            contigLocation.setAttribute("start", String.valueOf(location.bioStart()));
            contigLocation.setAttribute("end", String.valueOf(location.bioEnd()));
            // set these only once
            contigLocation.setReference("feature", cds);
            contigLocation.setReference("locatedOn", contig);
            if (onChromosome) {
                cds.setReference("chromosomeLocation", contigLocation);
            } else {
                cds.setReference("supercontigLocation", contigLocation);
            }
            // add to storage map
            cdsContigLocations.put(primaryIdentifier, contigLocation);
        }
    }

    /**
     * Place a feature on a protein.
     */
    Item getFeatureOnProtein(String type, String seqname, Location location) throws RuntimeException {
        Item protein = getProtein(seqname);
        Item feature = null;
        if (type.equals("protein_match")) {
            feature = createItem("ProteinMatch");
        } else if (type.equals("protein_hmm_match")) {
            feature = createItem("ProteinHmmMatch");
        } else {
            throw new RuntimeException("IPRSCAN GFF record type "+type+" is not supported by this loader.");
        }
        feature.setReference("protein", protein);
        // reference feature on new IM Location
        Item proteinLocation = createItem("Location");
        proteinLocation.setReference("feature", feature);
        proteinLocation.setAttribute("start", String.valueOf(location.bioStart()));
        proteinLocation.setAttribute("end", String.valueOf(location.bioEnd()));
        proteinLocation.setReference("locatedOn", protein);
        locations.add(proteinLocation);
        feature.addToCollection("locations", proteinLocation);
        return feature;
    }
    
    /**
     * Get/add a Gene Item, keyed by primaryIdentifier
     */
    Item getGene(String primaryIdentifier) {
        if (genes.containsKey(primaryIdentifier)) {
            return genes.get(primaryIdentifier);
        } else {
            Item gene = createItem("Gene");
            genes.put(primaryIdentifier, gene);
            gene.setAttribute("primaryIdentifier", primaryIdentifier);
            gene.setAttribute("secondaryIdentifier", getSecondaryIdentifier(primaryIdentifier));
            gene.setAttribute("assemblyVersion", assemblyVersion);
            gene.setAttribute("annotationVersion", annotationVersion);
            gene.setReference("organism", organism);
            gene.setReference("strain", strain);
            return gene;
        }
    }

    /**
     * Extract the secondaryIdentifier from an annotation feature primaryIdentifier.
     * Do not use this for chromosomes or supercontigs.
     */
    String getSecondaryIdentifier(String primaryIdentifier) throws RuntimeException {
        String secondaryIdentifier = null;
        secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, true);
        if (secondaryIdentifier == null) {
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
            chromosomes.put(primaryIdentifier, chromosome);
            chromosome.setAttribute("primaryIdentifier", primaryIdentifier);
            String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, false);
            if (secondaryIdentifier == null) {
                throw new RuntimeException("Could not get secondaryIdentifier for chromosome:"+primaryIdentifier);
            }
            chromosome.setAttribute("assemblyVersion", assemblyVersion);
            chromosome.setReference("organism", organism);
            chromosome.setReference("strain", strain);
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
            supercontigs.put(primaryIdentifier, supercontig);
            supercontig.setAttribute("primaryIdentifier", primaryIdentifier);
            String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, false);
            if (secondaryIdentifier == null) {
                throw new RuntimeException("Could not get secondaryIdentifier for supercontig:"+primaryIdentifier);
            }
            supercontig.setAttribute("assemblyVersion", assemblyVersion);
            supercontig.setReference("organism", organism);
            supercontig.setReference("strain", strain);
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
            proteins.put(primaryIdentifier, protein);
            protein.setAttribute("primaryIdentifier", primaryIdentifier);
            protein.setAttribute("secondaryIdentifier", getSecondaryIdentifier(primaryIdentifier));
            protein.setAttribute("assemblyVersion", assemblyVersion);
            protein.setAttribute("annotationVersion", annotationVersion);
            protein.setReference("organism", organism);
            protein.setReference("strain", strain);
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
            proteinDomains.put(primaryIdentifier, proteinDomain);
            proteinDomain.setAttribute("primaryIdentifier", primaryIdentifier);
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
            cds.setAttribute("assemblyVersion", assemblyVersion);
            cds.setAttribute("annotationVersion", annotationVersion);
            cds.setReference("organism", organism);
            cds.setReference("strain", strain);
            // add to cdses map
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
            mRNAs.put(primaryIdentifier, mRNA);
            mRNA.setAttribute("primaryIdentifier", primaryIdentifier);
            mRNA.setAttribute("secondaryIdentifier", getSecondaryIdentifier(primaryIdentifier));
            mRNA.setAttribute("assemblyVersion", assemblyVersion);
            mRNA.setAttribute("annotationVersion", annotationVersion);
            mRNA.setReference("organism", organism);
            mRNA.setReference("strain", strain);
            return mRNA;
        }
    }
    
    /**
     * Get/add a GeneFamily, keyed by primaryIdentifier
     */
    Item getGeneFamily(String primaryIdentifier) {
        if (geneFamilies.containsKey(primaryIdentifier)) {
            return geneFamilies.get(primaryIdentifier);
        } else {
            Item geneFamily = createItem("GeneFamily");
            geneFamilies.put(primaryIdentifier, geneFamily);
            geneFamily.setAttribute("primaryIdentifier", primaryIdentifier);
            return geneFamily;
        }
    }

    /**
     * Create a GeneFamilyAssignment and add it to the list.
     */
    void createGeneFamilyAssignment(Item geneFamily, Item gene, Item protein, double evalue, double score, double bestDomainScore) {
        Item gfa = createItem("GeneFamilyAssignment");
        geneFamilyAssignments.add(gfa);
        if (evalue > 0.0) gfa.setAttribute("evalue", String.valueOf(evalue));
        if (score > 0.0) gfa.setAttribute("score", String.valueOf(score));
        if (bestDomainScore > 0.0) gfa.setAttribute("bestDomainScore", String.valueOf(bestDomainScore));
        gfa.setReference("geneFamily", geneFamily);
        gfa.setReference("gene", gene);
        gfa.setReference("protein", protein);
        gene.addToCollection("geneFamilyAssignments", gfa);
        protein.addToCollection("geneFamilyAssignments", gfa);
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
            ontologyTerms.put(identifier, ontologyTerm);
            ontologyTerm.setAttribute("identifier", identifier);
            ontologyTerm.setReference("ontology", ontology);
            return ontologyTerm;
        }
    }

    /**
     * Return true if the current file contains only primary transcripts.
     */
    boolean isPrimary() {
        return getCurrentFile().getName().contains("_primary");
    }

    /**
     * For the given BioJava Sequence object, return an identifier to be used when creating the corresponding BioEntity.
     *
     * @param bioJavaSequence the Sequence
     * @return an identifier
     */
    static String getFastaIdentifier(AbstractSequence bioJavaSequence) {
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
     * Return an attribute for the given name ignoring case
     */
    static String getAttribute(FeatureI featureI, String name) {
        Map<String,String> attributeMap = featureI.getAttributes();
        for (String attributeName : attributeMap.keySet()) {
            if (attributeName.equalsIgnoreCase(name)) {
                return attributeMap.get(attributeName);
            }
        }
        return null;
    }

    /**
     * Print debug output on a GFF feature
     */
    static void printGFFDebugOutput(FeatureI featureI) {
        System.err.println("#### printGFFDebugOutput ####");
        System.err.println("seqname:\t" + featureI.seqname());
        System.err.println("location:\t" + featureI.location());
        System.err.println("type:\t" + featureI.type());
        System.err.println("ID:\t" + getAttribute(featureI, "ID"));
        System.err.println("Name:\t" + getAttribute(featureI, "Name"));
        System.err.println("Alias:\t" + getAttribute(featureI, "Alias"));
        System.err.println("Parent:\t" + getAttribute(featureI, "Parent"));
        System.err.println("Target:\t" + getAttribute(featureI, "Target"));
        System.err.println("Gap:\t" + getAttribute(featureI, "Gap"));
        System.err.println("Derives_from:\t" + getAttribute(featureI, "Derives_from"));
        System.err.println("Note:\t" + getAttribute(featureI, "Note"));
        System.err.println("Dbxref:\t" + getAttribute(featureI, "Dbxref"));
        System.err.println("Ontology_term:\t" + getAttribute(featureI, "Ontology_term"));
        System.err.println("IsCircular:\t" + getAttribute(featureI, "Is_circular"));
        System.err.println("alleles:\t" + getAttribute(featureI, "alleles"));
        System.err.println("symbol:\t" + getAttribute(featureI, "symbol"));
    }

}
