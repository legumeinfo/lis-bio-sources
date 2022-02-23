package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.io.IOException;
import java.io.InputStream;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.Properties;
import static java.util.Map.entry;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.genome.parsers.gff.Location;

/**
 * Loads data from an LIS GFF3 file, rather than using the core InterMine GFF3 loader, which is incompatible with the LIS Datastore.
 *
 * @author Sam Hokin
 */
public class LISGFF3FileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(LISGFF3FileConverter.class);

    // all Items to be stored
    Map<String,Item> features = new HashMap<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> proteinDomains = new HashMap<>();
    List<Item> ontologyAnnotations = new ArrayList<>();
    List<Item> locations = new ArrayList<>();

    // for distinguishing chromosomes from supercontigs
    DatastoreUtils dsu;

    // map GFF types to InterMine classes; be sure to include extras in the additions file!
    Map<String,String> featureClasses = Map.ofEntries(
                                                      entry("gene", "Gene"),
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
                                                      entry("rRNA_primary_transcript", "RRNA"),
                                                      entry("genetic_marker", "GeneticMarker")
                                                      );

    /**
     * Create a new LISGFF3FileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public LISGFF3FileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
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
        } else if (getCurrentFile().getName().endsWith(".gff3")) {
            processGFF3File(getCurrentFile().getPath());
	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        // collection attributes and references
        for (Item feature : features.values()) {
            if (assemblyVersion!=null) feature.setAttribute("assemblyVersion", assemblyVersion);
            if (annotationVersion!=null) feature.setAttribute("annotationVersion", annotationVersion);
            feature.setReference("organism", organism);
            feature.setReference("strain", strain);
        }
        storeCollectionItems();
        // local items
        store(features.values());
        store(ontologyTerms.values());
        store(proteinDomains.values());
        store(ontologyAnnotations);
        store(locations);
    }

    /**
     * Process a GFF3 file, referenced by filename because GFFReader doesn't have a method to parse a Reader.
     * Assumes that ID=full-yuck-LIS-identifier and Name=name/secondaryIdentifier
     */
    void processGFF3File(String filename) throws IOException {
        System.out.println("Processing "+getCurrentFile().getName());
        FeatureList featureList = GFF3Reader.read(filename);
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
                throw new RuntimeException("GFF3 type "+type+" is not associated with a class in the data model.");
            }
            // ID=glyma.Lee.gnm1.ann1.GlymaLee.02G198600;
            Item feature = getFeature(id, featureClass, location, seqname);
            // Name=GlymaLee.02G198600;
            if (name!=null) {
                feature.setAttribute("name", name);
                feature.setAttribute("secondaryIdentifier", name);
            }
            // Note=Cytochrome P450 superfamily protein%3B IPR001128 (Cytochrome P450)%3B GO:0005506 (iron ion binding)%2C GO:0020037 (heme binding)%2C ...
            if (note!=null) {
                feature.setAttribute("description", note);
            }
            // Parent is not required
            if (parent!=null) {
                Item parentItem = features.get(parent);
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
            // Dbxref=Gene3D:G3DSA:1.10.630.10,InterPro:IPR001128,InterPro:IPR002401,InterPro:IPR017972,PANTHER:PTHR24298,...
            // only genes
            if (featureClass.equals("Gene") && dbxref!=null) {
                String[] terms = dbxref.split(",");
                for (String term : terms) {
                    addProteinDomain(feature, term);
                }
            }
            // GeneticMarker gets type=SNP if length==1 and alleles
            if (featureClass.equals("GeneticMarker")) {
                if (location.length()==1) feature.setAttribute("type", "SNP");
                if (alleles!=null) feature.setAttribute("alleles", alleles);
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
            String[] parts = term.split(":");
            String identifier = parts[1]; // IPR001128
            Item proteinDomain = proteinDomains.get(identifier);
            if (proteinDomain==null) {
                proteinDomain = createItem("ProteinDomain");
                proteinDomain.setAttribute("primaryIdentifier", identifier);
                proteinDomains.put(identifier, proteinDomain);
            }
            gene.addToCollection("proteinDomains", proteinDomain);
        }
    }

    /**
     * Get/add a feature Item of the given class, keyed by primaryIdentifier.
     * Sets reference to chromosome or supercontig given by seqname.
     */
    Item getFeature(String primaryIdentifier, String className, Location location, String seqname) {
        if (features.containsKey(primaryIdentifier)) {
            return features.get(primaryIdentifier);
        } else if (dsu.isSupercontig(seqname)) {
            Item supercontig = features.get(seqname);
            if (supercontig==null) {
                // create new supercontig
                supercontig = createItem("Supercontig");
                supercontig.setAttribute("primaryIdentifier", seqname);
                features.put(seqname, supercontig);
            }
            // create new feature on supercontig
            Item feature = createItem(className);
            feature.setAttribute("primaryIdentifier", primaryIdentifier);
            feature.setReference("supercontig", supercontig);
            feature.setAttribute("length", String.valueOf(location.length()));
            // create new IM Location
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
            features.put(primaryIdentifier, feature);
            return feature;
        } else {
            Item chromosome = features.get(seqname);
            if (chromosome==null) {
                // create new chromosome
                chromosome = createItem("Chromosome");
                chromosome.setAttribute("primaryIdentifier", seqname);
                features.put(seqname, chromosome);
            }
            // create new feature on chromosome
            Item feature = createItem(className);
            feature.setAttribute("primaryIdentifier", primaryIdentifier);
            feature.setReference("chromosome", chromosome);
            feature.setAttribute("length", String.valueOf(location.length()));
            // create new IM Location
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
            features.put(primaryIdentifier, feature);
            return feature;
        }
    }
}
