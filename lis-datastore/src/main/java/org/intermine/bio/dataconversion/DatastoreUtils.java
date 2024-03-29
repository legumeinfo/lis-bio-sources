package org.intermine.bio.dataconversion;

import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import org.ncgr.datastore.Readme;
import org.ncgr.pubmed.PubMedSummary;

/**
 * Static utility methods for datastore loaders. 
 *
 * @author Sam Hokin
 */
public class DatastoreUtils {

    private static final String PUBMED_API_KEY = "48cb39fb23bf1190394ccbae4e1d35c5c809";
    
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
     * strain.gnmN.key4
     * strain.gnmN.annN.key4.other
     * strain.gnmN.syn.key4
     * strain.gnmN.mrk.markerset
     * strain.gnmN.annN.expr.samplestrain.author1_author2_year
     *
     * Disregard any of these:
     * strain.gwas.author1_author2_year
     * strains.map.author1_author2_year
     * strains.qtl.author1_author2_year
     */
    public static String extractAssemblyVersionFromCollection(String identifier) {
        String[] fields = identifier.split("\\.");
        if (fields[1].equals("gwas") || fields[1].equals("map") || fields[1].equals("qtl")) {
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
     * strain.assy.annN.key4.other
     * G19833.gnm1.ann1.pScz
     * G19833.gnm1.ann1.expr.Negro_jamapa.ORourke_Iniguez_2014
     *
     * Disregard any with less than four fields:
     * strain.gnmN.key4
     * strain.gwas.author1_author2_year
     * strains.map.author1_author2_year
     * strains.qtl.author1_author2_year
     *
     * Disregard these with 4 fields:
     * strain.gnmN.syn.key4
     * strain.gnmN.mrk.markerset
     */
    public static String extractAnnotationVersionFromCollection(String identifier) {
        String[] fields = identifier.split("\\.");
        if (fields.length<4) {
            return null;
        } else if (fields[2].equals("syn") || fields[2].equals("mrk")) {
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
     * Extract the LIS prefix from an annotation filename.
     * 0     1            2    3    4    5       6    7
     * medsa.XinJiangDaYe.gnm1.ann1.RKB9.iprscan.gff3.gz
     */
    public static String extractPrefixFromAnnotationFilename(String filename) {
        String[] pieces = filename.split("\\.");
        return pieces[0]+"."+pieces[1]+"."+pieces[2]+"."+pieces[3];
    }
    
    /**
     * Unescape some URL escaped characters used in GFF notes, etc.
     */
    public static String unescape(String s) {
        return s.
            replaceAll("%09","\t").replaceAll("%26","&").
            replaceAll("%2509","\t").replaceAll("%2526","&").
            replaceAll("%2C",",").replaceAll("%3B",";").replaceAll("%3D","=").
            replaceAll("%252C",",").replaceAll("%253B",";").replaceAll("%253D","=").
            replaceAll("%2c",",").replaceAll("%3b",";").replaceAll("%3d","=").
            replaceAll("%252c",",").replaceAll("%253b",";").replaceAll("%253d","=");
    }

    /**
     * Retrieve the PubMed ID from PubMed for a given DOI.
     *
     * @param doi the DOI of the publication
     * @return the PubMed ID, or 0 if not found 
     */
    public static int getPubMedId(String doi) throws IOException, ParserConfigurationException, SAXException {
        PubMedSummary summary = new PubMedSummary();
        summary.searchDOI(doi, PUBMED_API_KEY);
        if (summary.id!=0) {
            return summary.id;
        } else {
            return 0;
        }
    }
    
}
