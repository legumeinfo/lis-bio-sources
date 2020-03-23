package org.intermine.bio.dataconversion;

import java.io.InputStream;
import java.io.IOException;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Properties;

/**
 * Utility methods for data store.
 *
 * @author Sam Hokin
 */
public class DatastoreUtils {

    // defaults for LIS datasource
    public static final String DEFAULT_DATASOURCE_NAME = "Legume Information System";
    public static final String DEFAULT_DATASOURCE_URL = "https://legumeinfo.org/";
    public static final String DEFAULT_DATASOURCE_DESCRIPTION =
        "A collaborative, community resource to facilitate crop improvement by integrating genetic, genomic, and trait data across legume species.";
        
    static final String ORGANISM_PROP_FILE = "organism_config.properties";
    static final String DATASTORE_PROP_FILE = "datastore_config.properties";

    Map<String,String> taxonIdGenus = new HashMap<>();
    Map<String,String> taxonIdSpecies = new HashMap<>();
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
            InputStream orgPropsResource = DatastoreUtils.class.getClassLoader().getResourceAsStream(ORGANISM_PROP_FILE);
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
            genusSpeciesTaxonId.put(genusSpecies, taxonId);
        }
        
        // get datastore properties, like supercontig-matching strings
        Properties datastoreProps = new Properties();
        try {
            InputStream datastorePropsResource = DatastoreUtils.class.getClassLoader().getResourceAsStream(DATASTORE_PROP_FILE);
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
    public boolean isSupercontig(String taxonId, String strainIdentifier, String primaryIdentifier) {
        String key = taxonId+"."+strainIdentifier;
        for (String matchString : supercontigStrings.get(key)) {
            if (primaryIdentifier.contains(matchString)) return true;
        }
        return false;
        // String lc = primaryIdentifier.toLowerCase();
	// if (lc.contains("scaffold")
        //     || lc.contains("contig")
        //     || lc.contains("pilon")
        //     || primaryIdentifier.contains("Aipa")
        //     || primaryIdentifier.contains("Adur")) {
	//     return true;
        // }
	// // tricky ones
	// String[] parts = primaryIdentifier.split("\\.");
        // if (parts.length>=4) {
        //     // 0     1           2    3
        //     // cicar.CDCFrontier.gnm1.C11044140
        //     if (parts[3].length()==9 && parts[3].charAt(0)=='C') return true;
        //     // 0     1   2    3
        //     // glyma.Lee.gnm1.sc119
        //     if (parts[1].equals("Lee") && parts[3].startsWith("sc")) return true;
        //     // 0     1        2    3
        //     // glyso.PI483463.gnm1.sc255
        //     if (parts[0].equals("glyso") && parts[3].startsWith("sc")) return true;
        //     // 0     1            2    3
        //     // medtr.jemalong_A17.gnm5.MtrunA17Chr0c01
        //     if (parts[3].contains("Chr0c")) return true;
        // }
	// // it's a chromosome!
	// return false;
    }

    /**
     * Extract the annotation version from the given filename.
     *
     * gensp.strain------.assy.anno.keyy.content.ext
     * 0----.1-----------.2---.3---.4---.5------.6--
     * vigun.IT97K-499-35.gnm1.ann1.zb5D.protein.faa
     * vigun.IT97K-499-35.gnm1.ann1.zb5D.gene_models_main.gff3
     */
    public static String extractAnnotationVersion(String filename) {
        String[] fields = filename.split("\\.");
        if (fields[3].startsWith("ann")) {
            return fields[3];
        } else {
            return null;
        }
    }

    /**
     * Extract the assembly version from the given filename.
     *
     * gensp.strain------.assy.anno.keyy.content.ext
     * 0----.1-----------.2---.3---.4---.5------.6--
     * vigun.IT97K-499-35.gnm1.ann1.zb5D.protein.faa
     * vigun.IT97K-499-35.gnm1.ann1.zb5D.gene_models_main.gff3
     *
     * gensp.strain------.assy.keyy.content.ext
     * 0----.1-----------.2---.3---.4------.5--
     * vigun.IT97K-499-35.gnm1.QnBw.genome_main.fna
     */
    public static String extractAssemblyVersion(String filename) {
        String[] fields = filename.split("\\.");
        if (fields[2].startsWith("gnm")) {
            return fields[2];
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
     * Form a primaryIdentifier from the secondaryIdentifier, gensp, strain, assembly and annotation
     */
    public static String formPrimaryIdentifier(String gensp, String strain, String assembly, String annotation, String secondaryIdentifier) {
        return gensp+"."+strain+"."+assembly+"."+annotation+"."+secondaryIdentifier;
    }
}
