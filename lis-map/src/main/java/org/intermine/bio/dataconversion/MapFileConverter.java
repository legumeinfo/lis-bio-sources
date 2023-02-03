package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.ncgr.datastore.Readme;
import org.ncgr.datastore.validation.MapCollectionValidator;
import org.ncgr.zip.GZIPBufferedReader;

/**
 * Store genetic maps from tab-delimited files in LIS Datastore /maps/ collections.
 *
 * Actual GeneticMarker objects are NOT created here; rather, their names are stored and the
 * corresponding GeneticMarker objects (SequenceFeatures) are related by a post-processor.
 * 
 * @author Sam Hokin
 */
public class MapFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(MapFileConverter.class);

    // local Items to store
    Item geneticMap;
    List<Item> linkageGroupPositions = new LinkedList<>();
    Map<String,Item> linkageGroups = new HashMap<>();
 
    // validate the collection first by storing a flag
    boolean collectionValidated = false;

    /**
     * Create a new MapFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public MapFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     *
     * mixed.map.GmComposite2003/
     * ├── glyma.mixed.map.GmComposite2003.lg.tsv.gz
     * ├── glyma.mixed.map.GmComposite2003.mrk.tsv.gz
     * └── README.mixed.map.GmComposite2003.yml
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (!collectionValidated) {
            MapCollectionValidator validator = new MapCollectionValidator(getCurrentFile().getParent());
            validator.validate();
            if (!validator.isValid()) {
                throw new RuntimeException("Collection "+getCurrentFile().getParent()+" does not pass validation.");
            }
            collectionValidated = true;
        }
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
            // GeneticMap
            geneticMap = createItem("GeneticMap");
            geneticMap.setReference("organism", organism);
            geneticMap.setAttribute("primaryIdentifier", readme.genetic_map);
            geneticMap.setAttribute("synopsis", readme.synopsis);
            geneticMap.setAttribute("description", readme.description);
            if (readme.genotyping_platform!=null) geneticMap.setAttribute("genotypingPlatform", readme.genotyping_platform);
            if (readme.genotyping_method!=null) geneticMap.setAttribute("genotypingMethod", readme.genotyping_method);
            // store |-delimited list of genotypes
            String genotypes = "";
            for (String genotype : readme.genotype) {
                if (genotypes.length()>0) genotypes += "|";
                genotypes += genotype;
            }
            geneticMap.setAttribute("genotypes", genotypes);
        } else if (getCurrentFile().getName().endsWith("lg.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processLgFile();
        } else if (getCurrentFile().getName().endsWith("mrk.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processMrkFile();
        } else {
            System.out.println(" x skipping "+getCurrentFile().getName());
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        // standard collection items
        if (readme==null) {
            throw new RuntimeException("README not read. Aborting.");
        }
        storeCollectionItems();
        // add publication to Annotatables
        geneticMap.addToCollection("publications", publication);
        // reference genetic map (in case README not read first)
        for (Item linkageGroup : linkageGroups.values()) {
            linkageGroup.setReference("geneticMap", geneticMap);
        }
        // store 'em
        store(geneticMap);
        store(linkageGroups.values());
        store(linkageGroupPositions);
    }
    
    /**
     * Process an lg.tsv.gz file
     * 0                            1       2
     * #linkage_group               length  chromosome number (optional)
     * TT_Tifrunner_x_GT-C20_c-A01  176.02  8
     * TT_Tifrunner_x_GT-C20_c-A02  185.70  10
     * TT_Tifrunner_x_GT-C20_c-A03  192.46  2
     */
    void processLgFile() throws IOException {
        BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
        Set<Integer> numbers = new HashSet<>(); // numbers already used
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
            if (fields.length<2) {
                throw new RuntimeException("File "+getCurrentFile().getName()+" data line does not contain two required fields: linkage_group, position:"+line);
            }
            
            String identifier = fields[0].trim();
            double length = Double.parseDouble(fields[1]);
            Item linkageGroup = getLinkageGroup(identifier);
            linkageGroup.setAttribute("length", String.valueOf(length));
            if (fields.length>2) {
                int number = Integer.parseInt(fields[2]);
                numbers.add(number);
                linkageGroup.setAttribute("number", String.valueOf(number));
            } else {
                // find the first available number starting with 1
                // this presumes that non-numbered LGs are at the bottom of a file where others have numbers!
                for (int n=1; n<1000; n++) {
                    if (!numbers.contains(n)) {
                        numbers.add(n);
                        linkageGroup.setAttribute("number", String.valueOf(n));
                        break;
                    }
                }
            }
        }
        br.close();
    }

    /**
     * Process a mrk.tsv.gz file which gives marker positions on linkage groups.
     * #marker linkage_group       position
     * A053_2  GmComposite2003_A1  34.55
     * A064_3  GmComposite2003_A1  42.08
     * A082_1  GmComposite2003_A1  102.30
     */
    void processMrkFile() throws IOException {
        BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
            if (fields.length<3) {
                throw new RuntimeException("File "+getCurrentFile().getName()+" data line does not contain three required fields: marker, linkage_group, position:"+line);
            }
            String markerName = fields[0].trim();
            String lgId = fields[1].trim();
            Double position = Double.parseDouble(fields[2]);
            // linkage group
            Item linkageGroup = getLinkageGroup(lgId);
            // linkage group position for this marker
            Item lgPosition = createItem("LinkageGroupPosition");
            linkageGroupPositions.add(lgPosition);
            lgPosition.setAttribute("markerName", markerName);
            lgPosition.setAttribute("position", String.valueOf(position));
            lgPosition.setReference("linkageGroup", linkageGroup);
        }
        br.close();
    }

    /**
     * Return a new or existing LinkageGroup
     */
    Item getLinkageGroup(String identifier) {
        if (linkageGroups.containsKey(identifier)) {
            return linkageGroups.get(identifier);
        } else {
            Item linkageGroup = createItem("LinkageGroup");
            linkageGroups.put(identifier, linkageGroup);
            linkageGroup.setAttribute("identifier", identifier);
            return linkageGroup;
        }
    }
}
