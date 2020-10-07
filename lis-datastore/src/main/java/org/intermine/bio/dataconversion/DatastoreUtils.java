package org.intermine.bio.dataconversion;

import java.io.InputStream;
import java.io.IOException;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Properties;

/**
 * Utility non-static methods for datastore loaders. 
 *
 * @author Sam Hokin
 */
public class DatastoreUtils {
        
    static final String ORGANISM_PROP_FILE = "organism_config.properties";
    static final String DATASTORE_PROP_FILE = "datastore_config.properties";

    Map<String,String> taxonIdGenus = new HashMap<>();
    Map<String,String> taxonIdSpecies = new HashMap<>();
    Map<String,String> taxonIdGensp = new HashMap<>();
    Map<String,String> genspTaxonId = new HashMap<>();
    Map<String,String> genusSpeciesTaxonId = new HashMap<>();

    Map<String,List<String>> supercontigStrings = new HashMap<>();

    public DatastoreUtils() {
        // map taxonId to genus, species
        // 0     1    2
        // taxon.3702.genus=Arabidopsis
        // taxon.3702.species=thaliana
        // taxon.3702.uniprot=ARATH
        Properties orgProps = new Properties();
        try {
            // DEBUG
            // class_keys.properties
            // datastore_config.properties
            // genomic_keyDefs.properties
            // genomic_precompute.properties
            // genomic_priorities.properties
            // keyword_search.properties
            // objectstoresummary.config.properties
            // organism_config.properties
            // so_terms
            if (getClass().getClassLoader().getResourceAsStream("so_terms")==null) {
                System.err.println("Did not find so_terms");
            }
            if (getClass().getClassLoader().getResourceAsStream("class_keys.properties")==null) {
                System.err.println("Did not find class_keys.properties");
            }
            if (getClass().getClassLoader().getResourceAsStream("datastore_config.properties")==null) {
                System.err.println("Did not find datastore_config.properties");
            }
            //
            InputStream orgPropsResource = getClass().getClassLoader().getResourceAsStream(ORGANISM_PROP_FILE);
            if (orgPropsResource == null) {
                System.err.println("Did not find organism properties file:"+ORGANISM_PROP_FILE);
                System.exit(1);
            }
            orgProps.load(orgPropsResource);
            for (Object obj : orgProps.keySet()) {
                String key = (String) obj;
                String value = orgProps.getProperty(key);
                String[] dotparts = key.split("\\.");
                if (dotparts[0].equals("taxon")) {
                    String taxonId = dotparts[1];
                    if (dotparts[2].equals("genus")) {
                        taxonIdGenus.put(taxonId, value);
                    } else if (dotparts[2].equals("species")) {
                        taxonIdSpecies.put(taxonId, value);
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("Problem loading properties from:"+ORGANISM_PROP_FILE, e);
        }
        // map gensp, Genus_species to taxonId
        for (String taxonId : taxonIdGenus.keySet()) {
            String genus = taxonIdGenus.get(taxonId);
            String species = taxonIdSpecies.get(taxonId);
            String gensp = genus.substring(0,3).toLowerCase()+species.substring(0,2).toLowerCase();
            String genusSpecies = genus+"_"+species;
            genspTaxonId.put(gensp, taxonId);
            taxonIdGensp.put(taxonId, gensp);
            genusSpeciesTaxonId.put(genusSpecies, taxonId);
        }
        
        // get datastore properties, like supercontig-matching strings
        Properties datastoreProps = new Properties();
        try {
            InputStream datastorePropsResource = getClass().getClassLoader().getResourceAsStream(DATASTORE_PROP_FILE);
            if (datastorePropsResource == null) {
                System.err.println("Did not find datastore properties file:"+DATASTORE_PROP_FILE);
                System.exit(1);
            }
            datastoreProps.load(datastorePropsResource);
            for (Object obj : datastoreProps.keySet()) {
                String key = (String) obj;
                String[] dotparts = key.split("\\.");
                if (dotparts[0].equals("supercontig")) {
                    // supercontig.3870.Amiga=Chr00c,Foo
                    String value = datastoreProps.getProperty(key);
                    List<String> matchStrings = Arrays.asList(value.split(","));
                    supercontigStrings.put(dotparts[1]+"."+dotparts[2], matchStrings);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("Problem loading properties from:"+DATASTORE_PROP_FILE, e);
        }
    }

    /**
     * Get the taxon ID for a gensp string like "phavu".
     */
    public String getTaxonId(String gensp) {
        String taxonId = genspTaxonId.get(gensp);
        if (taxonId==null) {
            throw new RuntimeException("Taxon ID not available for "+gensp);
        }
        return taxonId;
    }

    /**
     * Get the gensp (like "phavu") for a taxon ID.
     */
    public String getGensp(String taxonId) {
        String gensp = taxonIdGensp.get(taxonId);
        if (gensp==null) {
            throw new RuntimeException("gensp value not available for Taxon ID "+gensp);
        }
        return gensp;
    }

    /**
     * Get the Genus for a gensp string like "phavu".
     */
    public String getGenus(String gensp) {
        String taxonId = getTaxonId(gensp);
        if (taxonIdGenus.containsKey(taxonId)) {
            return taxonIdGenus.get(taxonId);
        } else {
            throw new RuntimeException("Genus not available for taxon ID "+taxonId);
        }
    }
    
    /**
     * Get the species for a gensp string like "phavu".
     */
    public String getSpecies(String gensp) {
        String taxonId = genspTaxonId.get(gensp);
        if (taxonIdSpecies.containsKey(taxonId)) {
            return taxonIdSpecies.get(taxonId);
        } else {
            throw new RuntimeException("Species not available for taxon ID "+taxonId);
        }
    }

    /**
     * Determine whether the given primaryIdentifier is for a Supercontig (based on the given matching strings).
     * TODO: THIS IS HORRIBLY SIMPLIFIED, BUT A NEW STARTING POINT.
     * NOTE: non-static since we need supercontigStrings!
     */
    public boolean isSupercontig(String gensp, String strainIdentifier, String primaryIdentifier) {
        String key = gensp+"."+strainIdentifier;
	List<String> matchStrings = supercontigStrings.get(key);
	if (matchStrings==null) {
	    throw new RuntimeException("You must add a supercontig matching entry for "+key+" in datastore_config.properties.");
	}
	for (String matchString : matchStrings) {
	    if (primaryIdentifier.contains(matchString)) return true;
	}
	return false;
    }
}
