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

import org.apache.log4j.Logger;

import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.ncgr.datastore.Readme;

/**
 * Loads data from an LIS datastore gfa.tsv file, e.g. phalu.G27455.gnm1.ann1.JD7C.legfed_v1_0.M65K.gfa.tsv.
 *
 * @author Sam Hokin
 */
public class GFAFileConverter extends DatastoreFileConverter {
	
    static final String DEFAULT_VERSION = "legfed_v1_0";

    private static final Logger LOG = Logger.getLogger(GFAFileConverter.class);

    // local items
    Map<String,Item> geneFamilies = new HashMap<>();
    Map<String,Item> genes = new HashMap<>();
    Map<String,Item> proteins = new HashMap<>();

    /**
     * Create a new GFAFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GFAFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
            setStrain();
        } else if (getCurrentFile().getName().endsWith(".gfa.tsv")) {
            processGFAFile(reader);
	}
    }


    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        if (readme==null) {
            throw new RuntimeException("README file missing. Aborting.");
        }
        if (geneFamilies.size()==0) {
            throw new RuntimeException("No GFA file found. Aborting.");
        }
        // set collection references
        for (Item gene : genes.values()) {
            gene.setReference("organism", organism);
            gene.setReference("strain", strain);
        }
        for (Item protein : proteins.values()) {
            protein.setReference("organism", organism);
            protein.setReference("strain", strain);
        }
        // store
        storeCollectionItems();
        store(geneFamilies.values());
        store(genes.values());
        store(proteins.values());
    }

    /**
     * Process a gfa.tsv file which contains relationships between gene families, genes and proteins, along with a score value.
     *
     * 0     1      2    3    4    5           6    7   8
     * phalu.G27455.gnm1.ann1.JD7C.legfed_v1_0.M65K.gfa.tsv
     * 0                1
     * ScoreMeaning	e-value
     *
     * 0=gene                                 1=gene family        2=protein                                3=score
     * phalu.G27455.gnm1.ann1.tig000546640010 legfed_v1_0.L_00CL8T phalu.G27455.gnm1.ann1.tig000546640010.1 2.4e-68
     */
    void processGFAFile(Reader reader) throws IOException {
        System.out.println("Processing "+getCurrentFile().getName());
        String[] fileParts = getCurrentFile().getName().split("\\.");
        if (fileParts.length!=9) {
            System.err.println("WARNING: GFA file does not have the required 9 dot-separated parts: "+getCurrentFile().getName());
        }
        String version = DEFAULT_VERSION;
        if (fileParts.length>5) version = fileParts[5];
        String scoreMeaning = null;
        // spin through the file
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        int linenumber = 0;
        while ((line=br.readLine())!=null) {
            linenumber++;
            if (line.startsWith("#")) continue;
            String[] fields = line.split("\t");
            if (fields.length==2) {
                // header
                if (fields[0].equals("ScoreMeaning")) scoreMeaning = fields[1];
            } else if (fields.length>2) {
                String geneIdentifier = fields[0];
                String geneFamilyIdentifier = fields[1];
                String proteinIdentifier = fields[2];
                double score = 0.0;
                boolean hasScore = false;
                if (fields.length>3) {
                    hasScore = true;
                    score = Double.parseDouble(fields[3]);
                }
                // validation
                if (geneIdentifier==null || geneIdentifier.trim().length()==0) {
                    throw new RuntimeException("ERROR: Gene.primaryIdentifier="+geneIdentifier+" at line "+linenumber);
                }
                if (proteinIdentifier==null || proteinIdentifier.trim().length()==0) {
                    throw new RuntimeException("ERROR: Protein.primaryIdentifier="+proteinIdentifier+" at line "+linenumber);
                }
                // Gene Family
                Item geneFamily = getGeneFamily(geneFamilyIdentifier);
                geneFamily.setAttribute("version", version);
                // Gene
                Item gene = getGene(geneIdentifier);
                gene.setReference("geneFamily", geneFamily);
                if (hasScore) {
                    if (scoreMeaning!=null) gene.setAttribute("geneFamilyScoreMeaning", scoreMeaning);
                    gene.setAttribute("geneFamilyScore", String.valueOf(score));
                }
                // Protein
                Item protein = getProtein(proteinIdentifier);
                protein.setReference("geneFamily", geneFamily);
                if (hasScore) {
                    if (scoreMeaning!=null) protein.setAttribute("geneFamilyScoreMeaning", scoreMeaning);
                    protein.setAttribute("geneFamilyScore", String.valueOf(score));
                }
                // Gene.proteins collection
                gene.addToCollection("proteins", protein);
            }
        }
        br.close();
    }

    /**
     * Get/add a Gene Item, keyed by primaryIdentifier
     */
    public Item getGene(String primaryIdentifier) {
        if (genes.containsKey(primaryIdentifier)) {
            return genes.get(primaryIdentifier);
        } else {
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", primaryIdentifier);
            genes.put(primaryIdentifier, gene);
            return gene;
        }
    }

    /**
     * Get/add a Protein Item, keyed by primaryIdentifier
     */
    public Item getProtein(String primaryIdentifier) {
        if (proteins.containsKey(primaryIdentifier)) {
            return proteins.get(primaryIdentifier);
        } else {
            Item protein = createItem("Protein");
            protein.setAttribute("primaryIdentifier", primaryIdentifier);
            proteins.put(primaryIdentifier, protein);
            return protein;
        }
    }

    /**
     * Get/add a GeneFamily, keyed by identifier
     */
    public Item getGeneFamily(String identifier) {
        if (geneFamilies.containsKey(identifier)) {
            return geneFamilies.get(identifier);
        } else {
            Item geneFamily = createItem("GeneFamily");
            geneFamily.setAttribute("identifier", identifier);
            geneFamilies.put(identifier, geneFamily);
            return geneFamily;
        }
    }
}
