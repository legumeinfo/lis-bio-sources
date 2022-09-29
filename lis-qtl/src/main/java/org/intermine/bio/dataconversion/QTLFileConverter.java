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
import org.ncgr.datastore.validation.QTLCollectionValidator;
import org.ncgr.zip.GZIPBufferedReader;

/**
 * Store QTL study data from tab-delimited files in LIS Datastore /qtl/ collections.
 *
 * GeneticMarker objects are NOT created here; rather, their names are stored and the
 * corresponding GeneticMarker objects (SequenceFeatures) are related by a post-processor.
 * 
 * @author Sam Hokin
 */
public class QTLFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(QTLFileConverter.class);

    // local Items to store
    Item qtlStudy;
    Item geneticMap;
    List<Item> ontologyAnnotations = new LinkedList<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> qtls = new HashMap<>();
    Map<String,Item> traits = new HashMap<>();
    Map<String,Item> linkageGroups = new HashMap<>();

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
            // GeneticMap
            geneticMap = createItem("GeneticMap");
            geneticMap.setAttribute("primaryIdentifier", readme.genetic_map);
            // QTLStudy
            qtlStudy = createItem("QTLStudy");
            qtlStudy.setReference("organism", organism);
            qtlStudy.setAttribute("primaryIdentifier", readme.identifier);
            qtlStudy.setAttribute("synopsis", readme.synopsis);
            qtlStudy.setAttribute("description", readme.description);
            qtlStudy.setReference("geneticMap", geneticMap);
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
        geneticMap.addToCollection("publications", publication);
        qtlStudy.addToCollection("publications", publication);
        // associate QTLs with QTLStudy (in case README not read first)
        for (Item qtl : qtls.values()) {
            qtl.setReference("qtlStudy", qtlStudy);
        }
        // associate LinkageGroups with GeneticMap (in case README not read first)
        for (Item linkageGroup : linkageGroups.values()) {
            linkageGroup.setReference("geneticMap", geneticMap);
        }
        // store 'em
        store(qtlStudy);
        store(geneticMap);
        store(qtls.values());
        store(linkageGroups.values());
        store(traits.values());
        store(ontologyTerms.values());
        store(ontologyAnnotations);
    }
    
    /**
     * Process a qtlmrk.tsv file. The first three columns are required. 
     * Distinction column 3 is optional and not loaded into mines.
     * Note: markers are placed on linkage groups with the lis-map loader.
     *
     * qtl_name           trait_name      marker_name distinction
     * 0                  1               2           3
     * Seed linoleic 1-1  Seed linoleic   A082_1      flanking
     * Seed oleic 1-1     Seed oleic      A082_1      peak
     * Sprout yield 1-1   Sprout yield    A089_2      flanking
     */
    void processQTLMrkFile() throws IOException {
        BufferedReader br = GZIPBufferedReader.getReader(getCurrentFile());
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
            if (fields.length<3) {
                throw new RuntimeException(getCurrentFile().getName()+" has record with less than 3 fields:"+line);
            }
            // required columns
            String qtlName = fields[0];
            String traitName = fields[1];
            String markerName = fields[2];
            // 0 QTL.name
	    Item qtl = getQTL(qtlName);
            // 1 QTL.trait
            Item trait = getTrait(traitName);
            qtl.setReference("trait", trait);
            // 2 QTL.markerNames (append)
            if (qtlMarkers.containsKey(qtlName)) {
                String markerNames = qtlMarkers.get(qtlName) + "|" + markerName;
                qtlMarkers.put(qtlName, markerNames);
                qtl.setAttribute("markerNames", markerNames);
            } else {
                qtlMarkers.put(qtlName, markerName);
                qtl.setAttribute("markerNames", markerName);
            }
	}
        br.close();
    }

    /**

     * Process a qtl.tsv file. The first five columns are required.
     * The qtl_name may not be unique across the mine. Same for trait. We index identifier,qtlStudy for that reason.
     *
     * qtl_identifier       trait_name        lg                           left_end  right_end  qtl_peak  favorable_allele_source  lod  likelihood_ratio  marker_r2  total_r2  additivity
     * 0                    1                 2                            3         4          5         6                        7    8                 9          10        11
     * Early leaf spot 1-1  Early leaf spot   TT_Tifrunner_x_GT-C20_c-A08  100.7     102.9      102                                3.02 12.42             0.56                             
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
            String lgName = fields[2];
            double left = doubleOrZero(fields[3]);
            double right = doubleOrZero(fields[4]);
            // 0 QTL.name
            Item qtl = getQTL(qtlIdentifier);
            // 1 QTL.trait
            qtl.setReference("trait", getTrait(traitName));
            // 2 QTL.linkageGroup
            qtl.setReference("linkageGroup", getLinkageGroup(lgName));
            // 3 QTL.start
            qtl.setAttribute("start", String.valueOf(left));
            // 4 QTL.end
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
            // 0:Trait.name
            Item trait = getTrait(traitName);
            // 1:Trait.description (may be blank)
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
            // 0:Trait
            Item trait = getTrait(traitName);
            // 1:OntologyTerm
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
    Item getLinkageGroup(String identifier) {
        if (linkageGroups.containsKey(identifier)) {
            return linkageGroups.get(identifier);
        } else {
            Item linkageGroup = createItem("LinkageGroup");
            linkageGroup.setAttribute("identifier", identifier);
            linkageGroups.put(identifier, linkageGroup);
            return linkageGroup;
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
