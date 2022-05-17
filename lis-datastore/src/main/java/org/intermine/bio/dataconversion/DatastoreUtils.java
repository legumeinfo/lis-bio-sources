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

    Map<String,List<String>> chromosomePrefixes = new HashMap<>();
    Map<String,List<String>> supercontigPrefixes = new HashMap<>();

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
        
        // load datastore properties, like chromosome and supercontig prefixes
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
                if (dotparts[0].equals("chromosome")) {
                    // chromosome.medtr.A17=Chr,CP,MT
                    String value = datastoreProps.getProperty(key);
                    List<String> prefixes = Arrays.asList(value.split(","));
                    chromosomePrefixes.put(dotparts[1]+"."+dotparts[2], prefixes);
                } else if (dotparts[0].equals("supercontig")) {
                    // supercontig.medtr.A17=Chr0
                    String value = datastoreProps.getProperty(key);
                    List<String> prefixes = Arrays.asList(value.split(","));
                    supercontigPrefixes.put(dotparts[1]+"."+dotparts[2], prefixes);
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
     * Determine whether the given primaryIdentifier is for a Chromosome (based on the given matching strings).
     *
     * datastore_config.properties:
     * supercontig.medtr.A17=MtrunA17Chr0c
     * chromosome.medtr.A17=MtrunA17Chr
     *
     * 0     1   2    3=name
     * medtr.A17.gnm1.MtrunA17Chr0c03 is a supercontig, its name "MtrunA17Chr0c03" starts with "MtrunA17Chr0c"
     * medtr.A17.gnm5.MtrunA17Chr1 is a chromosome, its name "MtrunA17Chr1" starts with "MtrunA17Chr" and does NOT start with "MtrunA17Chr0"
     *
     * Default is to return false.
     */
    public boolean isChromosome(String primaryIdentifier) {
        String[] fields = primaryIdentifier.split("\\.");
        try {
            String gensp = fields[0];
            String strainIdentifier = fields[1];
            String assy = fields[2];
            String name = fields[3];
            String key = gensp+"."+strainIdentifier;
            List<String> chrPrefixes = chromosomePrefixes.get(key);
            List<String> scPrefixes = supercontigPrefixes.get(key);
            // default
            if (chrPrefixes==null) {
                return false;
            } else {
                // supercontig prefix match takes precedence
                if (scPrefixes!=null) {
                    for (String prefix : scPrefixes) {
                        if (name.startsWith(prefix)) return false;
                    }
                }
                // chromosome prefix match is required to return true
                for (String prefix : chrPrefixes) {
                    if (name.startsWith(prefix)) return true;
                }
            }
        } catch (ArrayIndexOutOfBoundsException ex) {
            throw new RuntimeException(primaryIdentifier+" does not have enough dot-delimited parts!");
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
     * 0      1    2   3
     * strain.assy.mrk.markerset
     * 0      1    2
     * strain.assy.key4
     *
     * Disregarded:
     * -.qtl.-
     * -.gwas.-
     * -.gen.-
     */
    public static String extractAssemblyVersionFromCollection(String identifier) {
        String[] fields = identifier.split("\\.");
        if (fields[1].equals("gen")) {
            return null;
        } else if (fields[1].equals("gwas")) {
            return null;
        } else if (fields[1].equals("qtl")) {
            return null;
        } else {
            return fields[1];
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
     *
     * Disregard:
     * 0      1    2   3
     * strain.assy.mrk.markerset
     */
    public static String extractAnnotationVersionFromCollection(String identifier) {
        String[] fields = identifier.split("\\.");
        if (fields[2].equals("mrk")) {
            return null;
        } else {
            return fields[2];
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
