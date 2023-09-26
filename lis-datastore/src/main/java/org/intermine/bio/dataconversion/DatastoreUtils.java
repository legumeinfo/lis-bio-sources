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

    /**
     * Return the Ensembl Plants / Plant Reactome name for a gene, given its LIS name
     *
     * R-ADU-1119289.1	Arginine degradation    Arachis duranensis      Aradu.16RQU.
     * R-AIP-9030908.1	Underwater shoot and internode elongation	Arachis ipaensis	Araip.QX18T.
     * R-CCA-1119495.1	Citrulline biosynthesis	Cajanus cajan	C.cajan_36043.1
     * R-CAR-1119394.1	Pantothenate and coenzyme A biosynthesis III	Cicer arietinum	Ca_10073
     * R-GMA-1119402.1	Phospholipid biosynthesis I	Glycine max	GLYMA_08G085800
     * R-LAN-1119342.1	Gamma-glutamyl cycle	Lupinus angustifolius	TanjilG_00390
     * R-MTR-1119438.1	Secologanin and strictosidine biosynthesis	Medicago truncatula	MTR_7g090310
     * R-PVU-1119370.1	Sterol biosynthesis	Phaseolus vulgaris	PHAVU_004G160000g
     * R-PSA-1119273.1	Lysine biosynthesis I	Pisum sativum	Psat4g004640
     * R-TPR-1119583.1	Phytocassane biosynthesis	Trifolium pratense	Tp57577_TGAC_v2_gene3601
     * R-VAN-1119374.1	Abscisic acid biosynthesis	Vigna angularis	LR48_Vigan05g014300
     * R-VRA-1119351.1	Mitochondrial pyruvate metabolism	Vigna radiata	Vradi07g24310
     * R-VUN-1119312.1	Photorespiration	Vigna unguiculata	Vigun04g097100.v1.2
     *
     * @param lisName the gene name from the LIS annotation GFF Name attribute
     * @return the corresponding Ensembl gene name, or null
     */
    public static String getEnsemblName(String lisName) {
        // Aradu.16RQU
        if (lisName.startsWith("Aradu")) {
            return lisName + ".";
        }
        // Araip.QX18T
        if (lisName.startsWith("Araip")) {
            return lisName + ".";
        }
        // cajca.C.cajan_36043 FIX!
        if (lisName.startsWith("cajca.C.cajan")) {
            return lisName.replace("cajca.", "") + ".1";
        }
        // C.cajan_36043
        if (lisName.startsWith("C.cajan")) {
            return lisName + ".1";
        }
        // cicar.CDCFrontier.Ca_06796 cicar.ICC4958.Ca_06796  FIX!
        if (lisName.startsWith("cicar")) {
            String[] dotparts = lisName.split("\\.");
            return (dotparts[2]);
        }
        // Ca_06796
        if (lisName.startsWith("Ca_")) {
            return lisName;
        }
        // Glyma.08G085800
        if (lisName.startsWith("Glyma")) {
            return lisName.replace("Glyma.", "GLYMA_");
        }
        // lupan.Lup003900
        if (lisName.startsWith("lupan.Lup")) {
            String number = lisName.replace("lupan.Lup", "");
            return "TanjilG_" + number;
        }
        // Medtr7g090310
        if (lisName.startsWith("Medtr")) {
            return lisName.replace("Medtr", "MTR_");
        }
        // Phvul.009G204800
        if (lisName.startsWith("Phvul")) {
            return lisName.replace("Phvul.", "PHAVU_") + "g";
        }
        // Psat4g004640
        if (lisName.startsWith("Psat")) {
            return lisName;
        }
        // tripr.gene36010 FIX?
        if (lisName.startsWith("tripr.gene")) {
            return "Tp57577_TGAC_V2_" + lisName.replace("tripr.", "");
        }
        // Vigan.05G014300
        if (lisName.startsWith("Vigan.")) {
            return "LR48_Vigan" + lisName.replace("Vigan.", "");
        }
        // Vradi07g24310
        if (lisName.startsWith("Vradi")) {
            return lisName;
        }
        // Vigun04g097100
        if (lisName.startsWith("Vigun")) {
            return lisName + ".v1.2";
        }
        // vigun.IT97K-499-35.Vigun04g097100 FIX!
        if (lisName.startsWith("vigun.")) {
            String[] dotparts = lisName.split("\\.");
            return dotparts[2] + ".v1.2";
        }
        // default
        return null;
    }
    
}
