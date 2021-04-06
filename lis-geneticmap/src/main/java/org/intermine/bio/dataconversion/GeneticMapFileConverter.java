package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Store GeneticMap/marker/linkage group data from tab-delimited files.
 * 
 * @author Sam Hokin
 */
public class GeneticMapFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GeneticMapFileConverter.class);

    // local things to store
    Map<String,Item> linkageGroups = new HashMap<>();
    Map<String,Item> markers = new HashMap<>();
    List<Item> linkageGroupPositions = new LinkedList<>();

    // global singleton Items
    Item geneticMap = createItem("GeneticMap");
    Item publication = createItem("Publication");
    
    /**
     * Create a new GeneticMapFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GeneticMapFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        // arahy.mixed.map.TT_Tifrunner_x_GT-C20_c.lg.tsv
        // arahy.mixed.map.TT_Tifrunner_x_GT-C20_c.mrk.tsv
        // README.TT_Tifrunner_x_GT-C20_c.yml
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
        } else if (getCurrentFile().getName().endsWith("lg.tsv")) {
            processLgFile(reader);
        } else if (getCurrentFile().getName().endsWith("mrk.tsv")) {
            processMrkFile(reader);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        // DatastoreFileConverter
        store(dataSource);
        store(dataSets.values());
        store(organisms.values());
        store(strains.values());
        // local
        store(geneticMap);
        store(publication);
        store(markers.values());
        store(linkageGroupPositions);
        store(linkageGroups.values());
    }
    
    /**
     * Process the README, which contains the GeneticMap details.
     *
     * README.TT_Tifrunner_x_GT-C20_c.yml
     * ---------------
     * identifier: TT_Tifrunner_x_GT-C20_c
     * subject: "Genetic map of Tifrunner x GT-C20 for the study of early..."
     * taxid: 3818
     * genotype: 
     * - Tifrunner
     * - GT-C20
     * description: "Whole-genome resequencing (WGRS) of a recombinant inbred line (RIL) mapping..."
     * publication_doi: 10.1111/pbi.12930
     * publication_title: "High-Density Genetic Map Using Whole-Genome..."
     */
    void processReadme(Reader reader) throws IOException {
        Readme readme = Readme.getReadme(reader);
        // check required stuff
        if (readme.identifier==null ||
            readme.taxid==null ||
            readme.subject==null ||
            readme.description==null ||
            readme.genotype==null ||
            readme.publication_doi==null ||
            readme.publication_title==null) {
            throw new RuntimeException("ERROR: a required field is missing from README. "+
                                       "Required fields are: identifier, taxid, subject, description, genotype, publication_doi, publication_title");
        }
        // Organism
        Item organism = getOrganism(Integer.parseInt(readme.taxid));
        // GeneticMap
        geneticMap.setReference("organism", organism);
        geneticMap.setAttribute("primaryIdentifier", readme.identifier);
        geneticMap.setAttribute("subject", readme.subject);
        geneticMap.setAttribute("mappingDescription", readme.description);
        // Strain = mapping parents
        for (String mappingParent : readme.genotype) {
            Item strain = getStrain(mappingParent, organism);
            geneticMap.addToCollection("mappingParents", strain);
        }
        // create and store Publication here
        publication = createItem("Publication");
        publication.setAttribute("doi", readme.publication_doi);
        publication.setAttribute("title", readme.publication_title);
        geneticMap.addToCollection("publications", publication);
        // override DataSet.description from README
        Item dataSet = getDataSet();
        dataSet.setAttribute("description", readme.description);
        dataSet.setReference("publication", publication);
    }

    /**
     * Process an lg.tsv file
     * 0                                1
     * #uniquename                      length
     * TT_Tifrunner_x_GT-C20_c-A01      176.02
     * TT_Tifrunner_x_GT-C20_c-A02      185.7
     * TT_Tifrunner_x_GT-C20_c-A03      192.46
     */
    void processLgFile(Reader reader) throws IOException {
        BufferedReader br = new BufferedReader(reader);
	Item dataSet = getDataSet();
        Item organism = getOrganism();
        int number = 0; // convenience numbering of LGs
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            number++;
            String[] fields = line.split("\t");
            String identifier = fields[0].trim();
            double length = Double.parseDouble(fields[1]);
            Item linkageGroup = createItem("LinkageGroup");
            linkageGroup.setAttribute("identifier", identifier);
            linkageGroup.setAttribute("number", String.valueOf(number));
            linkageGroup.setAttribute("length", String.valueOf(length));
            linkageGroup.setReference("organism", organism);
            linkageGroup.setReference("geneticMap", geneticMap);
            linkageGroup.addToCollection("dataSets", dataSet);
            linkageGroups.put(identifier, linkageGroup);
        }
        br.close();
    }

    /**
     * Process a mrk.tsv file
     * 0                1       2
     * #marker          mappos  lg
     * A01_859822	0	TT_Tifrunner_x_GT-C20_c-A01
     * B01_15102376	0.75	TT_Tifrunner_x_GT-C20_c-A01
     * A01_304818	2.18	TT_Tifrunner_x_GT-C20_c-A01
     */
    void processMrkFile(Reader reader) throws IOException {
        BufferedReader br = new BufferedReader(reader);
        Item dataSet = getDataSet();
        Item organism = getOrganism();
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
            if (fields.length<3) {
                throw new RuntimeException("File "+getCurrentFile().getName()+" data line does not contain three required fields: marker, linkagegroup, position:"+line);
            }
            String markerId = fields[0].trim();
            Double position = Double.parseDouble(fields[1]);
            String lgId = fields[2].trim();
            // linkage group
            Item linkageGroup = linkageGroups.get(lgId);
            if (linkageGroup==null) {
                linkageGroup = createItem("LinkageGroup");
                linkageGroup.setAttribute("identifier", lgId);
                linkageGroups.put(lgId, linkageGroup);
            }
            // marker
            Item marker = markers.get(markerId);
            if (marker==null) {
                marker = createItem("GeneticMarker");
                marker.setReference("organism", organism);
                marker.setAttribute("secondaryIdentifier", markerId);
                marker.addToCollection("dataSets", dataSet);
                markers.put(markerId, marker);
            }
            // marker-linkage group
            linkageGroup.addToCollection("markers", marker);
            // linkage group position
            Item lgPosition = createItem("LinkageGroupPosition");
            lgPosition.setReference("marker", marker);
            lgPosition.setReference("linkageGroup", linkageGroup);
            lgPosition.setAttribute("position", String.valueOf(position));
            marker.addToCollection("linkageGroupPositions", lgPosition);
            linkageGroupPositions.add(lgPosition);
        }
        br.close();
    }
}
