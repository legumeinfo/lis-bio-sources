package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2015-2016 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

/**
 * Encapsulates a single tab-delimited GWAS experiment file record.
 *
 * 0               1          2            3
 * #identifier     phenotype  marker       pvalue
 * Seed oil 4-g14  Seed oil   ss715591641  3.16E-09
 * Seed oil 4-g15  Seed oil   ss715591642  3.16E-08
 *
 * NOTE: identifier and phenotype are converted to capitalized lower-case.
 *
 * @author Sam Hokin, NCGR
 */
public class GWASFileRecord {

    // file record values
    String identifier;             // 0
    String phenotype;              // 1
    String marker;                 // 2
    double pvalue;                 // 3
    double lod = Double.MAX_VALUE; // 4, optional
    
    /**
     * Instantiate from a line in the GWAS file.
     */
    public GWASFileRecord(String line) {
        String[] parts = line.split("\t");
        if (parts.length<4) {
            System.err.println("ERROR: GWASFileRecord input has fewer than 4 fields:");
            System.err.println(line);
            System.exit(1);
        }
	identifier = parts[0].trim().toLowerCase();
        phenotype = parts[1].trim().toLowerCase();
        identifier = identifier.substring(0,1).toUpperCase() + identifier.substring(1);
        phenotype = phenotype.substring(0,1).toUpperCase() + phenotype.substring(1);
	marker = parts[2].trim();
        try {
            pvalue = Double.parseDouble(parts[3].trim());
        } catch (NumberFormatException ex) {
            System.err.println("Error parsing p value in line:");
            System.err.println(line);
            System.err.println(ex);
            System.exit(1);
        }
        if (parts.length>4 && !parts[4].toLowerCase().equals("na")) {
            try {
                lod = Double.parseDouble(parts[4].trim());
            } catch (NumberFormatException ex) {
                System.err.println("Error parsing lod value in line:");
                System.err.println(line);
                System.err.println(ex);
                System.exit(1);
            }
        }
    }

    /**
     * For diagnostics
     */
    public String toString() {
        return "[identifier,phenotype,marker,pvalue,lod]=["+identifier+","+phenotype+","+marker+","+pvalue+","+lod+"]";
    }
}
