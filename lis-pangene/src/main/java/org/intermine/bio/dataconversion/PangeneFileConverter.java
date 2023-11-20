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
 * Loads data from an LIS datastore hsh.tsv file, e.g. Phaseolus.pan1.X2PC.hsh.tsv.gz.
 * Each row contains the PanGeneSet identifier and a transcript (mRNA) identifier.
 *
 * Phaseolus.pan1.pan00001	phaac.Frijol_Bayo.gnm1.ann1.Phacu.CVR.011G222700.1
 * Phaseolus.pan1.pan00001	phaac.W6_15578.gnm2.ann1.Phacu.WLD.011G221300.1
 * Phaseolus.pan1.pan00001	phalu.G27455.gnm1.ann1.Pl11G0000365800.1
 *
 * @author Sam Hokin
 */
public class PangeneFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(PangeneFileConverter.class);

    // local things to store
    Map<String,Item> pangeneSets = new HashMap<>();
    Map<String,Item> transcripts = new HashMap<>();

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
        } else if (getCurrentFile().getName().endsWith(".hsh.tsv.gz")) {
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
        if (publication != null) {
            for (Item transcript : transcripts.values()) {
                transcript.addToCollection("publications", publication);
            }
            for (Item pangeneSet : pangeneSets.values()) {
                pangeneSet.addToCollection("publications", publication);
            }
        }
        // local items
        store(pangeneSets.values());
        store(transcripts.values());
    }

    /**
     * Process a hsh.tsv file which lists pan-gene set and transcript (mRNA) identifiers.
     *
     * 0                        1
     * identifier               transcript
     * Phaseolus.pan1.pan00001	phaac.Frijol_Bayo.gnm1.ann1.Phacu.CVR.011G222700.1
     */
    void processPangeneFile() throws IOException {
        // spin through the file
        BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        int linenumber = 0;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#")) continue;
            String[] fields = line.split("\t");
            Item pangeneSet = getPanGeneSet(fields[0]);
            Item transcript = getTranscript(fields[1]);
            pangeneSet.addToCollection("transcripts", transcript);
        }
        br.close();
    }

    /**
     * Get/add a transcript (mRNA) Item, keyed by identifier
     */
    public Item getTranscript(String identifier) {
        if (transcripts.containsKey(identifier)) {
            return transcripts.get(identifier);
        } else {
            // we store MRNA to avoid conflict with matching mRNAs loaded from GFF
            Item transcript = createItem("MRNA");
            transcript.setAttribute("primaryIdentifier", identifier);
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(identifier, true);
	    if (secondaryIdentifier!=null) transcript.setAttribute("secondaryIdentifier", secondaryIdentifier);
            transcripts.put(identifier, transcript);
            return transcript;
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
