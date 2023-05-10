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

import org.ncgr.datastore.validation.PangeneCollectionValidator;
import org.ncgr.zip.GZIPBufferedReader;


/**
 * Loads data from an LIS datastore clister.tsv file, e.g. Cicer.pan1.SV8C.clust.tsv.gz.
 * Each row contains the PanGeneSet identifier and tab-separated protein identifiers.
 *
 * Cicer.pan1.pan00001 cicar.CDCFrontier.gnm3.ann1.Ca2g036100.1 cicar.CDCFrontier.gnm3.ann1.Ca2g036500.1 ...
 * Cicer.pan1.pan00002 cicar.CDCFrontier.gnm3.ann1.Ca3g160700.1 cicar.ICC4958.gnm2.ann1.Ca_06696.1 ...
 *
 * @author Sam Hokin
 */
public class PangeneFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(PangeneFileConverter.class);

    // local things to store
    Map<String,Item> pangeneSets = new HashMap<>();
    Map<String,Item> genes = new HashMap<>();
    Map<String,Item> proteins = new HashMap<>();

    // toggle for validator
    boolean collectionValidated = false;

    /**
     * Create a new PangeneFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public PangeneFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (!collectionValidated) {
            PangeneCollectionValidator validator = new PangeneCollectionValidator(getCurrentFile().getParent());
            validator.validate();
            if (!validator.isValid()) {
                throw new RuntimeException("Collection "+getCurrentFile().getParent()+" does not pass validation.");
            }
            collectionValidated = true;
        }
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
        } else if (getCurrentFile().getName().endsWith(".clust.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processPangeneFile();
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
            for (Item pangeneSet : pangeneSets.values()) {
                pangeneSet.addToCollection("publications", publication);
            }
        }
        // local items
        store(pangeneSets.values());
        store(genes.values());
        store(proteins.values());
    }

    /**
     * Process a hsh.tsv file which contains relationships between pan-gene sets and proteins.
     *
     * 0                   1                                        2                                        ...
     * identifier          protein 1                                protein 2                                ...
     * Cicer.pan1.pan00001 cicar.CDCFrontier.gnm3.ann1.Ca2g036100.1 cicar.CDCFrontier.gnm3.ann1.Ca2g036500.1 ...
     */
    void processPangeneFile() throws IOException {
        // spin through the file
        BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        int linenumber = 0;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#")) continue;
            linenumber++;
            String[] fields = line.split("\t");
            Item pangeneSet = getPanGeneSet(fields[0]);
            for (int i=1; i<fields.length; i++) {
                Item protein = getProtein(fields[i]);
                pangeneSet.addToCollection("proteins", protein);
                Item gene = getGene(fields[i]);
                pangeneSet.addToCollection("genes", gene);
            }
        }
        br.close();
    }

    /**
     * Get/add a Gene Item, keyed by identifier, given a protein identifier
     *
     * @param proteinIdentifier the identifier of the corresponding protein
     */
    public Item getGene(String proteinIdentifier) {
        String[] fields = proteinIdentifier.split("\\.");
        String geneIdentifier = fields[0];
        for (int i=1; i<(fields.length-2); i++) {
            geneIdentifier += "." + fields[i];
        }
        if (genes.containsKey(geneIdentifier)) {
            return genes.get(geneIdentifier);
        } else {
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", geneIdentifier);
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(geneIdentifier, true);
	    if (secondaryIdentifier!=null) {
                gene.setAttribute("secondaryIdentifier", secondaryIdentifier);
            }
            genes.put(geneIdentifier, gene);
            return gene;
        }
    }

    /**
     * Get/add a Protein Item, keyed by identifier
     */
    public Item getProtein(String identifier) {
        if (proteins.containsKey(identifier)) {
            return proteins.get(identifier);
        } else {
            Item protein = createItem("Protein");
            protein.setAttribute("primaryIdentifier", identifier);
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(identifier, true);
	    if (secondaryIdentifier!=null) protein.setAttribute("secondaryIdentifier", secondaryIdentifier);
            proteins.put(identifier, protein);
            return protein;
        }
    }

    /**
     * Get/add a PanGeneSet, keyed by identifier
     */
    public Item getPanGeneSet(String identifier) {
        if (pangeneSets.containsKey(identifier)) {
            return pangeneSets.get(identifier);
        } else {
            Item pangeneSet = createItem("PanGeneSet");
            pangeneSet.setAttribute("primaryIdentifier", identifier);
            pangeneSets.put(identifier, pangeneSet);
            return pangeneSet;
        }
    }
}
