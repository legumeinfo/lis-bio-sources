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
 * Loads data from an LIS genetic marker GFF3 file.
 *
 * @author Sam Hokin
 */
public class MarkerGFF3FileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(MarkerGFF3FileConverter.class);

    // all Items to be stored
    Map<String,Item> chromosomes = new HashMap<>();
    Map<String,Item> supercontigs = new HashMap<>();
    Map<String,Item> geneticMarkers = new HashMap<>();
    List<Item> locations = new ArrayList<>();

    // for distinguishing chromosomes from supercontigs
    DatastoreUtils dsu;

    /**
     * Create a new MarkerGFF3FileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public MarkerGFF3FileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
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
            System.out.println("Processing "+getCurrentFile().getName());
            processMarkerGFF3File();
	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        if (geneticMarkers.size()==0) {
            throw new RuntimeException("No genetic markers loaded. Aborting.");
        }
        // set collection attributes and references
        for (Item geneticMarker : geneticMarkers.values()) {
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
        FeatureList featureList = GFF3Reader.read(getCurrentFile().getPath());
        for (FeatureI featureI : featureList) {
            String seqname = featureI.seqname();
            Location location = featureI.location();
            String type = featureI.type();
            // attributes
            String id = featureI.getAttribute("ID");
            String name = featureI.getAttribute("Name");
            String note = featureI.getAttribute("Note");
            String alleles = featureI.getAttribute("alleles");
            // check that id exists and matches collection
            if (id==null) {
                throw new RuntimeException("GFF line does not include ID: "+featureI.toString());
            }
            if (!type.equals("genetic_marker")) continue;
            // the following presumes that the README has already been read, which should be the case since capital R comes before lower case
            if (!matchesStrainAndAssembly(id)) {
                throw new RuntimeException("ID "+id+" does not match strain.assembly from collection "+readme.identifier);
            }
            // get associated class
            Item geneticMarker = getGeneticMarker(id, location, seqname);
            // Name
            if (name!=null) {
                geneticMarker.setAttribute("name", name);
                geneticMarker.setAttribute("secondaryIdentifier", name);
            }
            // Note is normally not present
            if (note!=null) {
                geneticMarker.setAttribute("description", note);
            }
            // GeneticMarker gets type=SNP if length==1 as well as alleles
            if (location.length()==1) geneticMarker.setAttribute("type", "SNP");
            if (alleles!=null) geneticMarker.setAttribute("alleles", alleles);
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
            geneticMarker.setAttribute("length", String.valueOf(location.length()));
            geneticMarkers.put(primaryIdentifier, geneticMarker);
            if (dsu.isSupercontig(seqname)) {
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
            } else {
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
            }
            return geneticMarker;
        }
    }
}
