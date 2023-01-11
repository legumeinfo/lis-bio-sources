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

/**
 * Loads data from an LIS datastore hsh.tsv file, e.g. glysp.mixed.pan2.TV81.hsh.tsv.
 *
 * glysp.mixed.pan2.SoyPan000001   glyma.Lee.gnm1.ann1.GlymaLee.16G153800.1
 * glysp.mixed.pan2.SoyPan000001   glyma.Wm82.gnm1.ann1.Glyma16g31020.2
 * glysp.mixed.pan2.SoyPan000001   glyma.Lee.gnm1.ann1.GlymaLee.16G154600.1
 *
 * @author Sam Hokin
 */
public class HSHFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(HSHFileConverter.class);

    // local things to store
    Map<String,Item> panGeneSets = new HashMap<>();
    Map<String,Item> genes = new HashMap<>();
    Map<String,Item> proteins = new HashMap<>();

    /**
     * Create a new HSHFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public HSHFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
        } else if (getCurrentFile().getName().endsWith(".hsh.tsv")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processHSHFile(reader);
	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        // standard collection items
        storeCollectionItems();
        // add publication to Annotatables
        if (publication!=null) {
            for (Item gene : genes.values()) {
                gene.addToCollection("publications", publication);
            }
            for (Item protein : proteins.values()) {
                protein.addToCollection("publications", publication);
            }
        }
        // local items
        store(panGeneSets.values());
        store(genes.values());
        store(proteins.values());
    }

    /**
     * Process a hsh.tsv file which contains relationships between pan-gene sets and proteins.
     *
     * 0     1     2    3    4   5
     * glysp.mixed.pan2.TV81.hsh.tsv.
     *
     * 0=PanGeneSet.identifier         1=Protein.primaryIdentifier
     * glysp.mixed.pan2.SoyPan000001   glyma.Lee.gnm1.ann1.GlymaLee.16G153800.1
     */
    void processHSHFile(Reader reader) throws IOException {
        String[] fileParts = getCurrentFile().getName().split("\\.");
        if (fileParts.length!=6) {
            throw new RuntimeException("HSH file does not have the required 6 dot-separated parts: "+getCurrentFile().getName());
        }
        String version = fileParts[2];
        // spin through the file
        String line = null;
        BufferedReader br = new BufferedReader(reader);
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#")) continue;
            String[] fields = line.split("\t");
            if (fields.length!=2) continue;
            String panGeneSetIdentifier = fields[0];
            String proteinIdentifier = fields[1];
            String[] proteinParts = proteinIdentifier.split("\\.");
            String geneIdentifier = proteinParts[0];
            for (int i=1; i<(proteinParts.length-1); i++) {
                geneIdentifier += "." + proteinParts[i];
            }
            // PanGeneSet
            Item panGeneSet = getPanGeneSet(panGeneSetIdentifier);
            panGeneSet.setAttribute("version", version);
            // Protein
            Item protein = getProtein(proteinIdentifier);
            protein.setReference("panGeneSet", panGeneSet);
            // Gene
            Item gene = getGene(geneIdentifier);
            gene.setReference("panGeneSet", panGeneSet);
            gene.addToCollection("proteins", protein);
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
            // phavu.Phvul.002G040500
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", primaryIdentifier);
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, true);
	    if (secondaryIdentifier!=null) gene.setAttribute("secondaryIdentifier", secondaryIdentifier);
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
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, true);
	    if (secondaryIdentifier!=null) protein.setAttribute("secondaryIdentifier", secondaryIdentifier);
            proteins.put(primaryIdentifier, protein);
            return protein;
        }
    }

    /**
     * Get/add a PanGeneSet, keyed by identifier
     */
    public Item getPanGeneSet(String identifier) {
        if (panGeneSets.containsKey(identifier)) {
            return panGeneSets.get(identifier);
        } else {
            Item panGeneSet = createItem("PanGeneSet");
            panGeneSet.setAttribute("identifier", identifier);
            panGeneSets.put(identifier, panGeneSet);
            return panGeneSet;
        }
    }
}
