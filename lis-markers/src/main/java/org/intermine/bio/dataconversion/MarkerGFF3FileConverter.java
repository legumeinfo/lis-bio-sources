package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
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
import java.util.Set;
import java.util.HashSet;
import java.util.Properties;
import static java.util.Map.entry;

import org.apache.commons.text.WordUtils;
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

import org.ncgr.datastore.validation.MarkerCollectionValidator;
import org.ncgr.zip.GZIPBufferedReader;

/**
 * Loads data from an LIS genetic marker GFF3 file.
 *
 * All features will be loaded as GeneticMarker. GeneticMarker.type will be set to the value in the GFF
 * type column (e.g. SNP), and, if the marker is a single base it will be forced to type=SNP.
 *
 * @author Sam Hokin
 */
public class MarkerGFF3FileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(MarkerGFF3FileConverter.class);
    private static final String TEMPFILE = "/tmp/marker.gff3";

    // all Items to be stored
    Map<String,Item> chromosomes = new HashMap<>();
    Map<String,Item> supercontigs = new HashMap<>();
    Map<String,Item> geneticMarkers = new HashMap<>();
    List<Item> locations = new ArrayList<>();

    // validate the collection first by storing a flag
    boolean collectionValidated = false;

    /**
     * Create a new MarkerGFF3FileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public MarkerGFF3FileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (!collectionValidated) {
            MarkerCollectionValidator validator = new MarkerCollectionValidator(getCurrentFile().getParent());
            validator.validate();
            if (!validator.isValid()) {
                throw new RuntimeException("Collection "+getCurrentFile().getParent()+" does not pass validation.");
            }
            collectionValidated = true;
        }
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
            // markers are mapped to a specific strain assembly
            setStrain();
        } else if (getCurrentFile().getName().endsWith(".gff3.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processMarkerGFF3File();
        } else {
            System.out.println("## - Skipping "+getCurrentFile().getName());
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        if (readme==null) {
            throw new RuntimeException("README file not read. Aborting.");
        }
        if (geneticMarkers.size()==0) {
            throw new RuntimeException("No genetic markers loaded. Aborting.");
        }
        // set collection attributes and references
        for (Item geneticMarker : geneticMarkers.values()) {
            geneticMarker.setAttribute("genotypingPlatform", readme.genotyping_platform);
            geneticMarker.setAttribute("assemblyVersion", assemblyVersion);
            geneticMarker.setReference("organism", organism);
            geneticMarker.setReference("strain", strain);
        }
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
        // add publication to Annotatables (but not chromosome/supercontig)
        if (publication!=null) {
            for (Item geneticMarker : geneticMarkers.values()) {
                geneticMarker.addToCollection("publications", publication);
            }
        }
        storeCollectionItems();
        // local items
        store(chromosomes.values());
        store(supercontigs.values());
        store(geneticMarkers.values());
        store(locations);
    }

    /**
     * Process a genetic marker GFF3 file, referenced by filename because GFFReader doesn't have a method to parse a Reader.
     * Assumes that ID=full-yuck-LIS-identifier and Name=name
     */
    void processMarkerGFF3File() throws IOException {
        if (readme==null) {
            throw new RuntimeException("README not read before "+getCurrentFile().getName()+". Aborting.");
        }
        // uncompress the gff3.gz file to a temp file
        File tempfile = new File(TEMPFILE);
        tempfile.delete();
        BufferedWriter writer = new BufferedWriter(new FileWriter(tempfile));
        BufferedReader reader = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        while ( (line=reader.readLine())!=null ) {
            writer.write(line);
            writer.newLine();
        }
        writer.close();
        // now load the uncompressed GFF3
        FeatureList featureList = GFF3Reader.read(TEMPFILE);
        for (FeatureI featureI : featureList) {
            String seqname = featureI.seqname();
            Location location = featureI.location();
            String type = featureI.type();
            // attributes (case-insensitive to initcap)
            String id = getAttribute(featureI, "ID");
            String name = getAttribute(featureI, "Name");
            String note = getAttribute(featureI, "Note");
            String alleles = getAttribute(featureI, "Alleles");
            String alias = getAttribute(featureI, "Alias");
            String symbol = getAttribute(featureI, "Symbol");
            String motif = getAttribute(featureI, "Motif");
            // check that id exists and matches collection
            if (id==null) {
                throw new RuntimeException("GFF line does not include ID: "+featureI.toString());
            }
            // the following presumes that the README has already been read, which should be the case since capital R comes before lower case
            if (!matchesStrainAndAssembly(id)) {
                throw new RuntimeException("ID "+id+" does not match strain.assembly from collection "+readme.identifier);
            }
            // GeneticMarker
            Item geneticMarker = getGeneticMarker(id, location, seqname);
            // type
            if (type.equals("genetic_marker")) {
                // auto-set type to SNP if one base
                if (location.length()==1) geneticMarker.setAttribute("type", "SNP");
            } else {
                // use the value in the type column
                geneticMarker.setAttribute("type", type);
            }
            // name
            if (name!=null) geneticMarker.setAttribute("name", name);
            // symbol
            if (symbol!=null) geneticMarker.setAttribute("symbol", symbol);
            // alias
            if (alias!=null) geneticMarker.setAttribute("alias", alias);
            // note is normally not present
            if (note!=null) geneticMarker.setAttribute("description", note);
            // alleles
            if (alleles!=null) geneticMarker.setAttribute("alleles", alleles);
            // motif
            if (motif!=null) geneticMarker.setAttribute("motif", motif);
        }
    }

    /**
     * Get/add a Chromosome Item, keyed by primaryIdentifier with secondaryIdentifier from primaryIdentifier.
     */
    public Item getChromosome(String primaryIdentifier) throws RuntimeException {
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
    public Item getSupercontig(String primaryIdentifier) throws RuntimeException {
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
     * Get/add a geneticMarker Item of the given class, keyed by primaryIdentifier.
     * Sets reference to chromosome or supercontig given by seqname.
     */
    Item getGeneticMarker(String primaryIdentifier, Location location, String seqname) {
        if (geneticMarkers.containsKey(primaryIdentifier)) {
            return geneticMarkers.get(primaryIdentifier);
        } else {
            Item geneticMarker = createItem("GeneticMarker");
            geneticMarker.setAttribute("primaryIdentifier", primaryIdentifier);
            geneticMarker.setAttribute("secondaryIdentifier", DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, false));
            geneticMarker.setAttribute("length", String.valueOf(location.length()));
            geneticMarkers.put(primaryIdentifier, geneticMarker);
            if (isChromosome(seqname)) {
                Item chromosome = getChromosome(seqname);
                geneticMarker.setReference("chromosome", chromosome);
                Item chromosomeLocation = createItem("Location");
                chromosomeLocation.setReference("feature", geneticMarker);
                if (location.isNegative()) {
                    chromosomeLocation.setAttribute("strand", "-1");
                } else {
                    chromosomeLocation.setAttribute("strand", "1");
                }
                chromosomeLocation.setAttribute("start", String.valueOf(location.bioStart()));
                chromosomeLocation.setAttribute("end", String.valueOf(location.bioEnd()));
                chromosomeLocation.setReference("locatedOn", chromosome);
                locations.add(chromosomeLocation);
                geneticMarker.setReference("chromosomeLocation", chromosomeLocation);
            } else if (isSupercontig(seqname)) {
                Item supercontig = getSupercontig(seqname);
                geneticMarker.setReference("supercontig", supercontig);
                Item supercontigLocation = createItem("Location");
                supercontigLocation.setReference("feature", geneticMarker);
                if (location.isNegative()) {
                    supercontigLocation.setAttribute("strand", "-1");
                } else {
                    supercontigLocation.setAttribute("strand", "1");
                }
                supercontigLocation.setAttribute("start", String.valueOf(location.bioStart()));
                supercontigLocation.setAttribute("end", String.valueOf(location.bioEnd()));
                supercontigLocation.setReference("locatedOn", supercontig);
                locations.add(supercontigLocation);
                geneticMarker.setReference("supercontigLocation", supercontigLocation);
            }
            return geneticMarker;
        }
    }

    /**
     * Return an attribute for either the given name ignoring case
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
}
