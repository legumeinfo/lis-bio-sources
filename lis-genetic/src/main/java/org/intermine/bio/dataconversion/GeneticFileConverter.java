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

import org.ncgr.datastore.Readme;

/**
 * Store GeneticMap/marker/linkage group data from tab-delimited files.
 *
 * Actual GeneticMarker objects are NOT created here; rather, their names are stored and the
 * corresponding GeneticMarker objects (SequenceFeatures) are loaded by a post-processor.
 * 
 * @author Sam Hokin
 */
public class GeneticFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GeneticFileConverter.class);

    // local Items to store
    Item geneticMap;
    List<Item> ontologyAnnotations = new LinkedList<>();
    List<Item> linkageGroupPositions = new LinkedList<>();
    List<Item> populations = new LinkedList<>();

    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> traits = new HashMap<>();
    Map<String,Item> qtls = new HashMap<>();
    Map<String,Item> linkageGroups = new HashMap<>();
    Map<String,String> qtlMarkers = new HashMap<>();
    Map<String,String> linkageGroupMarkers = new HashMap<>();

    /**
     * Create a new GeneticFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GeneticFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     *
     * MAGIC-2017.gen.Huynh_Ehlers_2018/
     * ├── README.MAGIC-2017.gen.Huynh_Ehlers_2018.yml
     * ├── vigun.MAGIC-2017.gen.Huynh_Ehlers_2018.lg.tsv
     * ├── vigun.MAGIC-2017.gen.Huynh_Ehlers_2018.mrk.tsv
     * ├── vigun.MAGIC-2017.gen.Huynh_Ehlers_2018.obo.tsv
     * ├── vigun.MAGIC-2017.gen.Huynh_Ehlers_2018.qtlmrk.tsv
     * └── vigun.MAGIC-2017.gen.Huynh_Ehlers_2018.qtl.tsv
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
            if (readme.genotype==null) {
                throw new RuntimeException("ERROR: README file must contain genotype.");
            }
            // GeneticMap
            geneticMap = createItem("GeneticMap");
            geneticMap.setReference("organism", organism);
            geneticMap.setAttribute("primaryIdentifier", readme.identifier);
            geneticMap.setAttribute("synopsis", readme.synopsis);
            geneticMap.setAttribute("description", readme.description);
            if (readme.genotyping_platform!=null) geneticMap.setAttribute("genotypingPlatform", readme.genotyping_platform);
            
            // Populations (from genotype)
            for (String genotype : readme.genotype) {
                Item population = createItem("Population");
                population.setAttribute("identifier", genotype);
                populations.add(population);
                geneticMap.addToCollection("populations", population);
            }
            if (publication!=null) geneticMap.addToCollection("publications", publication);
        } else if (getCurrentFile().getName().endsWith("lg.tsv")) {
            processLgFile(reader);
	} else if (getCurrentFile().getName().endsWith("qtlmrk.tsv")) {
	    processQTLMrkFile(reader);
        } else if (getCurrentFile().getName().endsWith("mrk.tsv")) {
            processMrkFile(reader);
        } else if (getCurrentFile().getName().endsWith("obo.tsv")) {
            processOboFile(reader);
	} else if (getCurrentFile().getName().endsWith("qtl.tsv")) {
            processQTLFile(reader);
	} else if (getCurrentFile().getName().endsWith("trait.tsv")) {
            processTraitFile(reader);
	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        // set references to collection items
        for (Item linkageGroup : linkageGroups.values()) {
            linkageGroup.setReference("organism", organism);
            linkageGroup.setReference("geneticMap", geneticMap);
        }
        for (Item qtl : qtls.values()) {
            qtl.setReference("organism", organism);
            qtl.setReference("geneticMap", geneticMap);
        }
        // store collection items
        storeCollectionItems();
        // store local Items
        store(geneticMap);
        store(populations);
        store(ontologyAnnotations);
        store(linkageGroupPositions);
        store(ontologyTerms.values());
        store(traits.values());
        store(qtls.values());
        store(linkageGroups.values());
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
        System.out.println("Processing "+getCurrentFile().getName());
        BufferedReader br = new BufferedReader(reader);
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
            linkageGroups.put(identifier, linkageGroup);
        }
        br.close();
    }

    /**
     * Process a mrk.tsv file which gives marker positions on linkage groups.
     * 0                1       2
     * A01_859822	A01     0.0
     * B01_15102376     A01     0.75
     * A01_304818	A01     2.18
     */
    void processMrkFile(Reader reader) throws IOException {
        System.out.println("Processing "+getCurrentFile().getName());
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
            if (fields.length<3) {
                throw new RuntimeException("File "+getCurrentFile().getName()+" data line does not contain three required fields: marker, linkagegroup, position:"+line);
            }
            String markerName = fields[0].trim();
            String lgId = fields[1].trim();
            Double position = Double.parseDouble(fields[2]);
            // linkage group
            Item linkageGroup = linkageGroups.get(lgId);
            if (linkageGroup==null) {
                linkageGroup = createItem("LinkageGroup");
                linkageGroup.setAttribute("identifier", lgId);
                linkageGroups.put(lgId, linkageGroup);
            }
            // marker
            if (linkageGroupMarkers.containsKey(lgId)) {
                String markerNames = linkageGroupMarkers.get(lgId) + "|" + markerName;
                linkageGroupMarkers.put(lgId, markerNames);
                linkageGroup.setAttribute("markerNames", markerNames);
            } else {
                linkageGroupMarkers.put(lgId, markerName);
                linkageGroup.setAttribute("markerNames", markerName);
            }
            // linkage group position
            Item lgPosition = createItem("LinkageGroupPosition");
            lgPosition.setAttribute("markerName", markerName);
            lgPosition.setAttribute("position", String.valueOf(position));
            lgPosition.setReference("linkageGroup", linkageGroup);
            linkageGroupPositions.add(lgPosition);
        }
        br.close();
    }
    
    /**
     * Process a qtlmrk.tsv file. distinction column is optional.
     * 0                    1             2
     * Early leaf spot 1-1  A08_35596996  flanking
     * Early leaf spot 1-1  A08_35776787  nearest
     */
    void processQTLMrkFile(Reader reader) throws IOException {
        System.out.println("Processing "+getCurrentFile().getName());
        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
            // required fields
	    String qtlId = fields[0];
	    String markerName = fields[1];
            //  optional fields
	    String distinction = null;
	    if (fields.length>2) distinction = fields[2];
            // QTL
	    Item qtl = qtls.get(qtlId);
	    if (qtl==null) {
                qtl = createItem("QTL");
                qtl.setAttribute("identifier", qtlId);
                qtls.put(qtlId, qtl);
            }
            // marker
            if (qtlMarkers.containsKey(qtlId)) {
                String markerNames = qtlMarkers.get(qtlId) + "|" + markerName;
                qtlMarkers.put(qtlId, markerNames);
                qtl.setAttribute("markerNames", markerNames);
            } else {
                qtlMarkers.put(qtlId, markerName);
                qtl.setAttribute("markerNames", markerName);
            }
	}
        br.close();
    }

    /**
     * Process a qtl.tsv file. The first two columns are required.
     * 0                    1                 2                            3         4          5         6                        7    8                 9          10        11
     * #qtl_name            trait_id          lg                           left_end  right_end  qtl_peak  favorable_allele_source  lod  likelihood_ratio  marker_r2  total_r2  additivity
     * Early leaf spot 1-1  Early leaf spot   TT_Tifrunner_x_GT-C20_c-A08  100.7     102.9      102                                3.02 12.42             0.56                             
     */
    void processQTLFile(Reader reader) throws IOException {
        System.out.println("Processing "+getCurrentFile().getName());
        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
            // required columns
            String qtlId = fields[0];
            String traitId = fields[1];
            // optional columns
            String lgId = null;
            double leftEnd = Double.MAX_VALUE;
            double rightEnd = Double.MAX_VALUE;
            double peak = Double.MAX_VALUE;
            String favorableAlleleSource = null;
            double lod = Double.MAX_VALUE;
            double likelihoodRatio = Double.MAX_VALUE;
            double markerR2 = Double.MAX_VALUE;
            double totalR2 = Double.MAX_VALUE;
            double additivity = Double.MAX_VALUE;
            if (fields.length>2) lgId = fields[2];
            if (fields.length>3) leftEnd = doubleOrZero(fields[3]);
            if (fields.length>4) rightEnd = doubleOrZero(fields[4]);
            if (fields.length>5) peak = doubleOrZero(fields[5]);
            if (fields.length>6) favorableAlleleSource = fields[6];
            if (fields.length>7) lod = doubleOrZero(fields[7]);
            if (fields.length>8) likelihoodRatio = doubleOrZero(fields[8]);
            if (fields.length>9) markerR2 = doubleOrZero(fields[9]);
            if (fields.length>10) totalR2 = doubleOrZero(fields[10]);
            if (fields.length>11) additivity = doubleOrZero(fields[11]);
            // Trait
            Item trait = traits.get(traitId);
            if (trait==null) {
                trait = createItem("Trait");
                trait.setAttribute("primaryIdentifier", traitId);
                traits.put(traitId, trait);
            }
            // QTL
            Item qtl = qtls.get(qtlId);
            if (qtl==null) {
                qtl = createItem("QTL");
                qtl.setAttribute("identifier", qtlId);
                qtls.put(qtlId, qtl);
            }
            qtl.setReference("trait", trait);
            if (leftEnd<Double.MAX_VALUE) qtl.setAttribute("start", String.valueOf(leftEnd));
            if (rightEnd<Double.MAX_VALUE) qtl.setAttribute("end", String.valueOf(rightEnd));
            if (peak<Double.MAX_VALUE) qtl.setAttribute("peak", String.valueOf(peak));
            if (lod<Double.MAX_VALUE) qtl.setAttribute("lod", String.valueOf(lod));
            if (likelihoodRatio<Double.MAX_VALUE) qtl.setAttribute("likelihoodRatio", String.valueOf(likelihoodRatio));
            if (markerR2<Double.MAX_VALUE) qtl.setAttribute("markerR2", String.valueOf(markerR2));
            // favorableAlleleSource 
            // totalR2
            // additivity = Double.MAX_VALUE;
            // LinkageGroup
            if (lgId!=null) {
                Item linkageGroup = linkageGroups.get(lgId);
                if (linkageGroup==null) {
                    linkageGroup = createItem("LinkageGroup");
                    linkageGroup.setAttribute("identifier", lgId);
                    linkageGroups.put(lgId, linkageGroup);
                }
                qtl.setReference("linkageGroup", linkageGroup);
            }
            // GeneticMap
        }
        br.close();
    }

    /**
     * Process a trait.tsv file, creating Trait records.
     * 0                           2
     * #trait_name                 description/method/etc.
     * Early leaf spot resistance  Leafs were photographed and spots were counted...
     */
    void processTraitFile(Reader reader) throws IOException {
        System.out.println("Processing "+getCurrentFile().getName());
        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
            String identifier = fields[0];
            String description = fields[1];
            // Trait
            Item trait = traits.get(identifier);
            if (trait==null) {
                trait = createItem("Trait");
                trait.setAttribute("primaryIdentifier", identifier);
                traits.put(identifier, trait);
            }
            trait.setAttribute("description", description);
        }
        br.close();
    }

    /**
     * Process an obo.tsv file.
     * 0            1
     * #trait       obo_term
     * Seed weight  TO:0000181
     */
    void processOboFile(Reader reader) throws IOException {
        System.out.println("Processing "+getCurrentFile().getName());
        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
	    if (fields.length<2) continue;
            if (fields[1]==null || fields[1].trim().length()==0) continue;
            String traitId = fields[0];
            String ontologyId = fields[1];
            // OntologyTerm
            Item ontologyTerm = ontologyTerms.get(ontologyId);
            if (ontologyTerm==null) {
                ontologyTerm = createItem("OntologyTerm");
                ontologyTerm.setAttribute("identifier", ontologyId);
                ontologyTerms.put(ontologyId, ontologyTerm);
            }
            // Trait
            Item trait = traits.get(traitId);
            if (trait==null) {
                trait = createItem("Trait");
                trait.setAttribute("primaryIdentifier", traitId);
                traits.put(traitId, trait);
            }
            // OntologyAnnotation
            Item ontologyAnnotation = createItem("OntologyAnnotation");
            ontologyAnnotation.setReference("subject", trait);
            ontologyAnnotation.setReference("ontologyTerm", ontologyTerm);
            ontologyAnnotations.add(ontologyAnnotation);
        }
        br.close();
    }

    /**
     * Return double if length>0 else 0.00.
     */
    static double doubleOrZero(String field) throws NumberFormatException {
        if (field.length()>0) {
            return Double.parseDouble(field);
        } else {
            return 0.00;
        }
    }
}
