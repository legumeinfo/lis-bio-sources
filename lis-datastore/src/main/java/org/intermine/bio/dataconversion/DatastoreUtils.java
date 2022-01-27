package org.intermine.bio.dataconversion;

import java.io.File;
import java.io.InputStream;
import java.io.IOException;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Properties;

/**
 * Utility static and non-static methods for datastore loaders. 
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

    /**
     * Constructor for methods that require loading properties files.
     */
    public DatastoreUtils() {
        Properties orgProps = new Properties();
        try {
            // load the organism properties into maps
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
        
        // load datastore properties, like supercontig-matching strings
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
            throw new RuntimeException("Taxon ID not available for gensp="+gensp);
        }
        return taxonId;
    }

    /**
     * Get the gensp (like "phavu") for a taxon ID.
     */
    // public String getGensp(String taxonId) {
    //     String gensp = taxonIdGensp.get(taxonId);
    //     if (gensp==null) {
    //         throw new RuntimeException("gensp value not available for Taxon ID "+gensp);
    //     }
    //     return gensp;
    // }

    /**
     * Get the Genus for a gensp string like "phavu".
     */
    // public String getGenus(String gensp) {
    //     String taxonId = getTaxonId(gensp);
    //     if (taxonIdGenus.containsKey(taxonId)) {
    //         return taxonIdGenus.get(taxonId);
    //     } else {
    //         throw new RuntimeException("Genus not available for taxon ID "+taxonId);
    //     }
    // }
    
    /**
     * Get the species for a gensp string like "phavu".
     */
    // public String getSpecies(String gensp) {
    //     String taxonId = genspTaxonId.get(gensp);
    //     if (taxonIdSpecies.containsKey(taxonId)) {
    //         return taxonIdSpecies.get(taxonId);
    //     } else {
    //         throw new RuntimeException("Species not available for taxon ID "+taxonId);
    //     }
    // }

    /**
     * Determine whether the given primaryIdentifier is for a Supercontig (based on the given matching strings).
     * 0     1    2    3
     * glyma.Wm82.gnm4.Gm01 is a chromosome, it's name does not contain a supercontig-identifying string
     * glyma.Wm82.gnm4.scaffold_99 is a supercontig because its name contains 'scaffold'
     */
    public boolean isSupercontig(String primaryIdentifier) {
        String[] fields = primaryIdentifier.split("\\.");
        String gensp = fields[0];
        String strainIdentifier = fields[1];
        String assy = fields[2];
        String name = fields[3];
        String key = gensp+"."+strainIdentifier;
	List<String> matchStrings = supercontigStrings.get(key);
	if (matchStrings==null) {
	    throw new RuntimeException("You must add a supercontig matching entry for "+key+" in datastore_config.properties.");
	}
	for (String matchString : matchStrings) {
	    if (name.contains(matchString)) return true;
	}
	return false;
    }

    /**
     * Extract the gensp string from the given identifier.
     * glyma.Wm82.gnm1.Chr04
     * glyma.Wm82.gnm1.ann1.Gene01
     */
    public static String extractGensp(String identifier) {
        String[] fields = identifier.split("\\.");
        if (fields.length>=4) {
            return fields[0];
        } else {
            return null;
        }
    }

    /**
     * Extract the collection identifier from a README filename.
     */
    public static String extractCollectionFromReadme(File readme) {
        String[] fields = readme.getName().split("\\.");
        String collection = fields[1];
        for (int i=2; i<fields.length; i++) collection += "."+fields[i];
        return collection;
    }

    /**
     * Extract the Strain identifier from the given collection identifier.
     * strain.assy.anno.KEY4
     */
    public static String extractStrainIdentifierFromCollection(String identifier) {
	String[] fields = identifier.split("\\.");
        return fields[0];
    }

    /**
     * Extract the Strain identifier from the given feature identifier.
     * gensp.strain.assy.secondaryIdentifier
     * gensp.strain.assy.anno.secondaryIdentifier
     */
    public static String extractStrainIdentifierFromFeature(String identifier) {
	String[] fields = identifier.split("\\.");
        return fields[1];
    }

    /**
     * Extract the assembly version from the given collection identifier.
     * 0      1    2    3
     * strain.assy.anno.key4
     * 0       1    2
     * strain.assy.key4
     */
    public static String extractAssemblyVersionFromCollection(String identifier) {
        String[] fields = identifier.split("\\.");
        if (fields.length==4 && fields[3].length()==4) {
            return fields[1];
        } else if (fields.length==3 && fields[2].length()==4) {
            return fields[1];
        } else {
            return null;
        }
    }

    /**
     * Extract the assembly version from the given feature identifier.
     * 0     1      2    3
     * gensp.strain.assy.secondaryIdentifier
     * 0     1      2    3    4
     * gensp.strain.assy.anno.secondaryIdentifier
     */
    public static String extractAssemblyVersionFromFeature(String identifier) {
	String[] fields = identifier.split("\\.");
        if (fields.length==4 || fields.length==5) {
            return fields[2];
        } else {
            return null;
        }
    }

    /**
     * Extract the annotation version from the given collection identifier.
     * 0      1    2    3
     * strain.assy.anno.key4
     */
    public static String extractAnnotationVersionFromCollection(String identifier) {
        String[] fields = identifier.split("\\.");
        if (fields.length==4 && fields[3].length()==4) {
            return fields[2];
        } else {
            // are there other scenarios?
            return null;
        }
    }

    /**
     * Extract the annotation version from the given feature identifier.
     * 0     1      2    3    4
     * gensp.strain.assy.anno.secondaryIdentifier
     */
    public static String extractAnnotationVersionFromFeature(String identifier) {
	String[] fields = identifier.split("\\.");
        if (fields.length==4) {
            return fields[3];
        } else {
            return null;
        }
    }

    /**
     * Extract the KEY4 portion of the given collection identifier.
     * strain.assy.xxx.key4, etc.
     */
    public static String extractKEY4(String identifier) {
        String[] fields = identifier.split("\\.");
        int i = fields.length -1;
        if (fields[i].length()==4) {
            return fields[i];
        } else {
            return null;
        }
    }

    /**
     * Extract the secondaryIdentifier from a full-yuck LIS identifier.
     * Set isAnnotationFeature=true if this is an annotation feature (at least five dot-separated parts) as opposed to an assembly feature (at least four dot-separated parts).
     *
     * isAnnotationFeature==true:
     * 0      1      2   3 
     * genesp.strain.gnm.secondaryIdentifier
     *
     * isAnnotationFeature==false:
     * 0      1      2   3   4
     * genesp.strain.gnm.ann.secondaryIdentifier
     *
     * @param   lisIdentifier the LIS full-yuck identifier
     * @param   isAnnotationFeature true if this is an annotation feature with five or more dot-separated parts
     * @returns the secondaryIdentifier
     */
    public static String extractSecondaryIdentifier(String lisIdentifier, boolean isAnnotationFeature) {
        String[] fields = lisIdentifier.split("\\.");
        if (isAnnotationFeature && fields.length>=5) {
            String secondaryIdentifier = fields[4];
            for (int i=5; i<fields.length; i++) {
        	secondaryIdentifier += "."+fields[i];
            }
            return secondaryIdentifier;
        } else if (!isAnnotationFeature && fields.length>=4) {
            String secondaryIdentifier = fields[3];
            for (int i=4; i<fields.length; i++) {
        	secondaryIdentifier += "."+fields[i];
            }
            return secondaryIdentifier;
        } else {
            return null;
        }
    }

    /**
     * Extract the gene identifier from a protein identifier, which is assumed to be [geneIdentifier].n.
     * 
     * @param proteinIdentifier the protein LIS identifier
     * @returns the corresponding gene LIS identifier
     */
    public static String extractGeneIdentifierFromProteinIdentifier(String proteinIdentifier) {
	String[] fields = proteinIdentifier.split("\\.");
	String geneIdentifier = fields[0];
	for (int i=1; i<(fields.length-1); i++) {
	    geneIdentifier += "."+fields[i];
	}
	return geneIdentifier;
    }
}
