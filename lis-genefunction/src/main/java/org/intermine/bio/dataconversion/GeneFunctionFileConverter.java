package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.Reader;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.HashSet;
import java.util.Properties;
import static java.util.Map.entry;
    
import org.apache.log4j.Logger;

import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.Util;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.ncgr.datastore.Readme;
import org.ncgr.zip.GZIPFastaReader;
import org.ncgr.zip.GZIPBufferedReader;

/**
 * Loads data from an LIS datastore gene_functions collection
 * Files types loaded are:
 *   README.yaml
 *
 * @author Andrew Farmer
 */
public class GeneFunctionFileConverter extends DatastoreFileConverter {

    // spit out debug lines if not null and ID starts with this
    private static final String DEBUG_ID = null;

    private static final Logger LOG = Logger.getLogger(GeneFunctionFileConverter.class);

    Map<String,Item> publications = new HashMap<>();

    /**
     * Create a new GeneFunctionFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GeneFunctionFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
        } else if (getCurrentFile().getName().endsWith(".citations.txt")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processCitations();
        } 
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException, RuntimeException {
        // store our Items
        //storeCollectionItems();
        store(publications.values());
        store(authors.values());
    }

    /**
     * Process a citations file
     */
    void processCitations() throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(getCurrentFile()));
        String line;
        while ((line = reader.readLine()) != null) {
            String[] fields = line.split("\t");
            String doi = fields[0];
            String citation = fields[2];
	    if (publications.get(doi) == null) {
                try {
                    Item publication = createPublication(doi);
                    publication.setAttribute("citation", citation);
                    publications.put(doi, publication);
                } catch (Exception ex) {
                    throw new RuntimeException(ex);
                }
	    }
        }
        reader.close();
    }

}
