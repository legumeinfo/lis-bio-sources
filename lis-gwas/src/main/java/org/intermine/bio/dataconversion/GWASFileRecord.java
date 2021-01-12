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
            System.err.println("ERROR: GWASFileRecord input has fewer than 4 fields:"+line);
            System.exit(1);
        }
	identifier = parts[0].trim();
        phenotype = parts[1].trim();
	marker = parts[2].trim();
	pvalue = Double.parseDouble(parts[3].trim());
        if (parts.length>4 && !parts[4].toLowerCase().equals("na")) {
            lod = Double.parseDouble(parts[4].trim());
        }
    }

    /**
     * For diagnostics
     */
    public String toString() {
        String str = "";
	str += " identifier="+identifier;
        str += " phenotype="+phenotype;
        str += " marker="+marker;
        str += " pvalue="+pvalue;
        str += " lod="+lod;
        return str;
    }
}
