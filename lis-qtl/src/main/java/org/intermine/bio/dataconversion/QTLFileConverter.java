package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.ncgr.datastore.Readme;
import org.ncgr.datastore.validation.QTLCollectionValidator;
import org.ncgr.zip.GZIPBufferedReader;

/**
 * Store QTL study data from tab-delimited files in LIS Datastore /qtl/ collections.
 *
 * GeneticMarker objects are NOT created here; rather, their names are stored along with the genetic map,
 * and the corresponding GeneticMarker SequenceFeature object are related by a post-processor.
 * 
 * @author Sam Hokin
 */
public class QTLFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(QTLFileConverter.class);

    // local Items to store
    Item qtlStudy;
    List<Item> ontologyAnnotations = new ArrayList<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> qtls = new HashMap<>();
    Map<String,Item> traits = new HashMap<>();
    Map<String,Item> geneticMaps = new HashMap<>();
    Map<String,Item> linkageGroups = new HashMap<>(); // keyed by geneticMapIdentifier_linkageGroupIdentifier

    // utility
    Map<String,String> qtlMarkers = new HashMap<>();

    // validate the collection first by storing a flag
    boolean collectionValidated = false;

    /**
     * Create a new QTLFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public QTLFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     *
     * MAGIC-2017.qtl.Huynh_Ehlers_2018/
     * ├── README.MAGIC-2017.qtl.Huynh_Ehlers_2018.yml
     * ├── vigun.MAGIC-2017.qtl.Huynh_Ehlers_2018.obo.tsv
     * ├── vigun.MAGIC-2017.qtl.Huynh_Ehlers_2018.qtlmrk.tsv
     * └── vigun.MAGIC-2017.qtl.Huynh_Ehlers_2018.qtl.tsv
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (!collectionValidated) {
            QTLCollectionValidator validator = new QTLCollectionValidator(getCurrentFile().getParent());
            validator.validate();
            if (!validator.isValid()) {
                throw new RuntimeException("Collection "+getCurrentFile().getParent()+" does not pass validation.");
            }
            collectionValidated = true;
        }
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
            // QTLStudy
            qtlStudy = createItem("QTLStudy");
            qtlStudy.setReference("organism", organism);
            qtlStudy.setAttribute("primaryIdentifier", readme.identifier);
            qtlStudy.setAttribute("synopsis", readme.synopsis);
            qtlStudy.setAttribute("description", readme.description);
	    // NEEDS MODEL UPDATE:
            // Parse comma-delimited genetic maps, for each:
            // qtlStudy.addToCollection("geneticMaps", geneticMap);
	    // store |-delimited list of genotypes
            String genotypes = "";
            for (String genotype : readme.genotype) {
                if (genotypes.length()>0) genotypes += "|";
                genotypes += genotype;
            }
            qtlStudy.setAttribute("genotypes", genotypes);
	} else if (getCurrentFile().getName().endsWith("qtlmrk.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
	    processQTLMrkFile();
        } else if (getCurrentFile().getName().endsWith("obo.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processOboFile();
	} else if (getCurrentFile().getName().endsWith("qtl.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processQTLFile();
	} else if (getCurrentFile().getName().endsWith("trait.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processTraitFile();
        } else {
            System.out.println("## - Skipping "+getCurrentFile().getName());
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
        qtlStudy.addToCollection("publications", publication);
        for (Item geneticMap : geneticMaps.values()) {
            geneticMap.addToCollection("publications", publication);
        }
        // associate QTLs with QTLStudy (in case README not read first)
        for (Item qtl : qtls.values()) {
            qtl.setReference("qtlStudy", qtlStudy);
        }
        // associate Traits with QTLStudy (in case README not read first)
        for (Item trait : traits.values()) {
            trait.setReference("qtlStudy", qtlStudy);
        }
        // store 'em
        store(qtlStudy);
        store(geneticMaps.values());
        store(qtls.values());
        store(linkageGroups.values());
        store(traits.values());
        store(ontologyTerms.values());
        store(ontologyAnnotations);
    }
    
    /**
     * Process a qtl.tsv file. The first five columns are required.
     * The qtl_name may not be unique across the mine. Same for trait. We index identifier,qtlStudy for that reason.
     *
     * 0                    1                2                        3              4      5     6         7                        8    9                 10         11        12
     * #qtl_identifier      trait_name       genetic_map              linkage_group  start  end   qtl_peak  favorable_allele_source  lod  likelihood_ratio  marker_r2  total_r2  additivity
     * Early leaf spot 1-1  Early leaf spot  TT_Tifrunner_x_GT-C20_c  A08            100.7  102.9 102                                3.02 12.42             0.56
     */
    void processQTLFile() throws IOException {
        BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
            if (fields.length<5) {
                throw new RuntimeException("ERROR: file "+getCurrentFile().getName()+" does not have five required fields: qtl_name, trait_name, lg, left_end, right_end here:"+line);
            }
            // required columns
            String qtlIdentifier = fields[0];
            String traitName = fields[1];
            String geneticMapIdentifier = fields[2];
            String linkageGroupIdentifier = fields[3];
            double left = doubleOrZero(fields[4]);
            double right = doubleOrZero(fields[5]);
            // QTL.name
            Item qtl = getQTL(qtlIdentifier);
            // QTL.trait
            qtl.setReference("trait", getTrait(traitName));
            // QTL.linkageGroup (with geneticMap)
            qtl.setReference("linkageGroup", getLinkageGroup(linkageGroupIdentifier, geneticMapIdentifier));
            // QTL.start
            qtl.setAttribute("start", String.valueOf(left));
            // QTL.end
            qtl.setAttribute("end", String.valueOf(right));
            // optional columns (skip 6, 10, 11)
            if (fields.length>5) qtl.setAttribute("peak", String.valueOf(doubleOrZero(fields[5])));
            if (fields.length>7) qtl.setAttribute("lod", String.valueOf(doubleOrZero(fields[7])));
            if (fields.length>8) qtl.setAttribute("likelihoodRatio", String.valueOf(doubleOrZero(fields[8])));
            if (fields.length>9) qtl.setAttribute("markerR2", String.valueOf(doubleOrZero(fields[9])));
        }
        br.close();
    }

    /**
     * Process a qtlmrk.tsv file. We ignore the distinction column.
     * 0                  1               2                3       4
     * #qtl_identifier    trait_name      genetic_map      marker  distinction
     * Seed linoleic 1-1  Seed linoleic   GmComposite1999  A082_1  flanking
     */
    void processQTLMrkFile() throws IOException {
        BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
            String qtlIdentifier = fields[0];
            String traitName = fields[1];
            String geneticMapIdentifier = fields[2];
            String markerName = fields[3];
            // QTL
	    Item qtl = getQTL(qtlIdentifier);
            // QTL.markerNames (append)
            if (qtlMarkers.containsKey(qtlIdentifier)) {
                String markerNames = qtlMarkers.get(qtlIdentifier) + "|" + markerName;
                qtlMarkers.put(qtlIdentifier, markerNames);
                qtl.setAttribute("markerNames", markerNames);
            } else {
                qtlMarkers.put(qtlIdentifier, markerName);
                qtl.setAttribute("markerNames", markerName);
            }
	}
        br.close();
    }

    /**
     * Process a trait.tsv file, creating Trait records. Both columns are required.
     * Trait name is likely not unique across mine, so we index identifier,qtlStudy.
     *
     * 0                           2
     * #trait_name                 description/method/etc.
     * Early leaf spot resistance  Leafs were photographed and spots were counted...
     */
    void processTraitFile() throws IOException {
        BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
            // required fields
            String traitName = fields[0];
            String description = fields[1];
            // Trait.name
            Item trait = getTrait(traitName);
            // Trait.description (may be blank)
            if (description!=null && description.length()>0) trait.setAttribute("description", description);
        }
        br.close();
    }

    /**
     * Process an obo.tsv file, both fields required.
     *
     * 0            1
     * #trait       obo_term
     * Seed weight  TO:0000181
     */
    void processOboFile() throws IOException {
        BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
	    if (fields.length<2) {
                throw new RuntimeException("File "+getCurrentFile().getName()+" has less than two fields in line:"+line);
            }
            if (fields[1]==null || fields[1].trim().length()==0) {
                throw new RuntimeException("File "+getCurrentFile().getName()+" has empty OBO term in line:"+line);
            }
            // required fields
            String traitName = fields[0];
            String oboTerm = fields[1];
            // Trait
            Item trait = getTrait(traitName);
            // OntologyTerm
            Item ontologyTerm = getOntologyTerm(oboTerm);
            // OntologyAnnotation
            Item ontologyAnnotation = createItem("OntologyAnnotation");
            ontologyAnnotation.setReference("subject", trait);
            ontologyAnnotation.setReference("ontologyTerm", ontologyTerm);
            ontologyAnnotations.add(ontologyAnnotation);
        }
        br.close();
    }

    /**
     * Return a new or existing QTL Item keyed by identifier
     */
    Item getQTL(String identifier) {
        if (qtls.containsKey(identifier)) {
            return qtls.get(identifier);
        } else {
            Item qtl = createItem("QTL");
            qtl.setAttribute("identifier", identifier);
            qtls.put(identifier, qtl);
            return qtl;
        }
    }
    
    /**
     * Return a new or existing LinkageGroup Item
     */
    Item getLinkageGroup(String linkageGroupIdentifier, String geneticMapIdentifier) {
        String key = geneticMapIdentifier+"_"+linkageGroupIdentifier;
        if (linkageGroups.containsKey(key)) {
            return linkageGroups.get(key);
        } else {
            Item geneticMap = getGeneticMap(geneticMapIdentifier);
            Item linkageGroup = createItem("LinkageGroup");
            linkageGroup.setAttribute("identifier", linkageGroupIdentifier);
            linkageGroup.setReference("geneticMap", geneticMap);
            linkageGroups.put(key, linkageGroup);
            return linkageGroup;
        }
    }

    /**
     * Return a new or existing GeneticMap Item
     */
    Item getGeneticMap(String primaryIdentifier) {
        if (geneticMaps.containsKey(primaryIdentifier)) {
            return geneticMaps.get(primaryIdentifier);
        } else {
            Item geneticMap = createItem("GeneticMap");
            geneticMap.setAttribute("primaryIdentifier", primaryIdentifier);
            geneticMaps.put(primaryIdentifier, geneticMap);
            return geneticMap;
        }
    }

    /**
     * Return a new or existing Trait Item keyed by name.
     */
    Item getTrait(String name) {
        if (traits.containsKey(name)) {
            return traits.get(name);
        } else {
            Item trait = createItem("Trait");
            trait.setAttribute("name", name);
            traits.put(name, trait);
            return trait;
        }
    }

    /**
     * Return a new or existing OntologyTerm Item.
     */
    Item getOntologyTerm(String identifier) {
        if (ontologyTerms.containsKey(identifier)) {
            return ontologyTerms.get(identifier);
        } else {
            Item ontologyTerm = createItem("OntologyTerm");
            ontologyTerm.setAttribute("identifier", identifier);
            ontologyTerms.put(identifier, ontologyTerm);
            return ontologyTerm;
        }
    }

    /**
     * Return the collection identifier portion of the current file.
     *     0     1          2   3                 4   5
     * └── vigun.MAGIC-2017.qtl.Huynh_Ehlers_2018.qtl.tsv
     */
    String getCollectionIdentifier() {
        String[] parts = getCurrentFile().getName().split("\\.");
        return parts[1]+"."+parts[2]+"."+parts[3];
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
