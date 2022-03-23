package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2020 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 */

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;
import java.util.Map;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.apache.log4j.Logger;

/**
 * Load gene-pathway associations from three-column LIS files.
 *
 * @author Sam Hokin
 */
public class PathwayFileConverter extends DatastoreFileConverter {
    private static final Logger LOG = Logger.getLogger(PathwayFileConverter.class);
    
    private Map<String,Item> pathways = new HashMap<>();      // keyed by Pathway.identifier
    private Map<String,Item> genes = new HashMap<>();         // keyed by Gene.name

    /**
     * Constructor
     *
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public PathwayFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
            setStrain();
        } else if (getCurrentFile().getName().endsWith("pathway.tsv")) {
            processPathwayFile(reader);
        }
    }

    /**
     * Process the pathway.tsv file
     */
    void processPathwayFile(Reader reader) throws IOException {
        System.out.println("Processing "+getCurrentFile().getName());
        // spin through the file
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#")) continue;
            String[] fields = line.split("\t");
            if (fields.length==2) {
                // not doing anything with header
            } else if (fields.length==3) {
                String pathwayIdentifier = fields[0];
                String pathwayName = fields[1];
                String geneIdentifier = fields[2];
                Item gene = genes.get(geneIdentifier);
                if (gene==null) {
                    gene = createItem("Gene");
                    gene.setAttribute("primaryIdentifier", geneIdentifier);
                    genes.put(geneIdentifier, gene);
                }
                Item pathway = pathways.get(pathwayIdentifier);
                if (pathway==null) {
                    pathway = createItem("Pathway");
                    pathway.setAttribute("identifier", pathwayIdentifier);
                    pathway.setAttribute("name", pathwayName);
                    pathways.put(pathwayIdentifier, pathway);
                }
                gene.addToCollection("pathways", pathway);
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    public void close() throws ObjectStoreException {
        if (readme==null) {
            throw new RuntimeException("README file missing. Aborting.");
        }
        if (pathways.size()==0) {
            throw new RuntimeException("No pathway file found. Aborting.");
        }
        // augment our Items with README stuff
        for (Item gene : genes.values()) {
            gene.setReference("organism", organism);
            gene.setReference("strain", strain);
        }
        storeCollectionItems();
        store(genes.values());
        store(pathways.values());
    }
}
