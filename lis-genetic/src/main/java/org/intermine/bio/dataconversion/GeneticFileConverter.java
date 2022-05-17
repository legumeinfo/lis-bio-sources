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
 * Store GWAS or genetic map data from tab-delimited files in LIS Datastore /genetic/ collections.
 *
 * Actual GeneticMarker objects are NOT created here; rather, their names are stored and the
 * corresponding GeneticMarker objects (SequenceFeatures) are related by a post-processor.
 * 
 * @author Sam Hokin
 */
public class GeneticFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GeneticFileConverter.class);

    // local Items to store
    Item gwas;
    Item qtlStudy;
    List<Item> gwasResults = new LinkedList<>();
    List<Item> ontologyAnnotations = new LinkedList<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> qtls = new HashMap<>();
    Map<String,Item> traits = new HashMap<>();
    Map<String,Item> linkageGroups = new HashMap<>();

    // utility
    Map<String,String> qtlMarkers = new HashMap<>();

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
            if (getCurrentFile().getName().contains(".gen.")) {
                // QTLStudy
                qtlStudy = createItem("QTLStudy");
                qtlStudy.setReference("organism", organism);
                qtlStudy.setAttribute("primaryIdentifier", readme.identifier);
                qtlStudy.setAttribute("synopsis", readme.synopsis);
                qtlStudy.setAttribute("description", readme.description);
                if (readme.genetic_map!=null) qtlStudy.setAttribute("geneticMap", readme.genetic_map);
                // form |-delimited list of genotypes
                String genotypes = "";
                for (String genotype : readme.genotype) {
                    if (genotypes.length()>0) genotypes += "|";
                    genotypes += genotype;
                }
                qtlStudy.setAttribute("genotypes", genotypes);
            } else if (getCurrentFile().getName().contains(".gwas.")) {
                // GWAS
                if (readme.genotyping_platform==null) {
                    throw new RuntimeException("ERROR: "+getCurrentFile().getName()+" must contain genotyping_platform.");
                }
                gwas = createItem("GWAS");
                gwas.setReference("organism", organism);
                gwas.setAttribute("primaryIdentifier", readme.identifier);
                gwas.setAttribute("synopsis", readme.synopsis);
                gwas.setAttribute("description", readme.description);
                gwas.setAttribute("genotypingPlatform", readme.genotyping_platform);
                if (readme.genotyping_method!=null) gwas.setAttribute("genotypingMethod", readme.genotyping_method);
                // form |-delimited list of genotypes
                String genotypes = "";
                for (String genotype : readme.genotype) {
                    if (genotypes.length()>0) genotypes += "|";
                    genotypes += genotype;
                }
                gwas.setAttribute("genotypes", genotypes);
            }
	} else if (getCurrentFile().getName().endsWith("qtlmrk.tsv")) {
            System.out.println("Processing "+getCurrentFile().getName());
	    processQTLMrkFile(reader);
        } else if (getCurrentFile().getName().endsWith("obo.tsv")) {
            System.out.println("Processing "+getCurrentFile().getName());
            processOboFile(reader);
	} else if (getCurrentFile().getName().endsWith("qtl.tsv")) {
            System.out.println("Processing "+getCurrentFile().getName());
            processQTLFile(reader);
	} else if (getCurrentFile().getName().endsWith("trait.tsv")) {
            System.out.println("Processing "+getCurrentFile().getName());
            processTraitFile(reader);
	} else if (getCurrentFile().getName().endsWith("result.tsv")) {
            System.out.println("Processing "+getCurrentFile().getName());
            processResultFile(reader);
	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        if (readme==null) {
            throw new RuntimeException("README not read. Aborting.");
        }
        storeCollectionItems();
        // QTLStudy
        if (qtlStudy!=null) {
            qtlStudy.addToCollection("publications", publication);
            for (Item qtl : qtls.values()) {
                qtl.setReference("qtlStudy", qtlStudy);
                qtl.addToCollection("publications", publication);
            }
            for (Item linkageGroup : linkageGroups.values()) {
                linkageGroup.addToCollection("publications", publication);
            }
            store(qtlStudy);
            store(qtls.values());
            store(linkageGroups.values());
        }
        // GWAS
	if (gwas!=null) {
            gwas.addToCollection("publications", publication);
            for (Item gwasResult : gwasResults) gwasResult.setReference("gwas", gwas);
            store(gwas);
            store(gwasResults);
        }
        // Both 
        for (Item trait : traits.values()) {
            trait.addToCollection("publications", publication);
        }
        store(traits.values());
        store(ontologyTerms.values());
        store(ontologyAnnotations);
    }
    
    /**
     * Process a qtlmrk.tsv file. linkage group column 3 is optional; distinction column 4 is optional and not read here.
     * 0                  1               2       3                   4
     * Seed linoleic 1-1  Seed linoleic   A082_1  GmComposite1999_A1  flanking
     * Seed oleic 1-1     Seed oleic      A082_1  GmComposite1999_A1  peak
     * Sprout yield 1-1   Sprout yield    A089_2  GmComposite1999_B1  flanking
     */
    void processQTLMrkFile(Reader reader) throws IOException {
        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
            if (fields.length<3) {
                throw new RuntimeException(getCurrentFile().getName()+" has record with less than 3 fields:"+line);
            }
            // 0 QTL
	    String qtlName = fields[0];
	    Item qtl = getQTL(qtlName);
            // 1 QTL.trait
            String traitName = fields[1];
            Item trait = getTrait(traitName);
            qtl.setReference("trait", trait);
            // 2 QTL.markerNames
	    String markerName = fields[2];
            if (qtlMarkers.containsKey(qtlName)) {
                String markerNames = qtlMarkers.get(qtlName) + "|" + markerName;
                qtlMarkers.put(qtlName, markerNames);
                qtl.setAttribute("markerNames", markerNames);
            } else {
                qtlMarkers.put(qtlName, markerName);
                qtl.setAttribute("markerNames", markerName);
            }
            // 3 QTL.linkageGroup
            if (fields.length>3) {
                String lgId = fields[3];
                Item linkageGroup = getLinkageGroup(lgId);
                qtl.setReference("linkageGroup", linkageGroup);
            }
	}
        br.close();
    }

    /**

     * Process a qtl.tsv file. The first two columns are required.
     * └── vigun.MAGIC-2017.gen.Huynh_Ehlers_2018.qtl.tsv
     * The qtl_name may not be unique across the mine, so we'll prefix the collection identifier for the primaryIdentifier. Same for trait.
     * The first five columns are required.
     * 0                    1                 2                            3         4          5         6                        7    8                 9          10        11
     * #qtl_name            trait_name        lg                           left_end  right_end  qtl_peak  favorable_allele_source  lod  likelihood_ratio  marker_r2  total_r2  additivity
     * Early leaf spot 1-1  Early leaf spot   TT_Tifrunner_x_GT-C20_c-A08  100.7     102.9      102                                3.02 12.42             0.56                             
     */
    void processQTLFile(Reader reader) throws IOException {
        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
            if (fields.length<5) {
                throw new RuntimeException("ERROR: file "+getCurrentFile().getName()+" does not have five required fields: qtl_name, trait_name, lg, left_end, right_end here:"+line);
            }
            // 0 QTL
            Item qtl = getQTL(fields[0]);
            // 1 QTL.trait
            qtl.setReference("trait", getTrait(fields[1]));
            // 2 QTL.linkageGroup
            qtl.setReference("linkageGroup", getLinkageGroup(fields[2]));
            // 3 QTL.start
            qtl.setAttribute("start", String.valueOf(doubleOrZero(fields[3])));
            // 4 QTL.end
            qtl.setAttribute("end", String.valueOf(doubleOrZero(fields[4])));
            // optional columns
            if (fields.length>5) qtl.setAttribute("peak", String.valueOf(doubleOrZero(fields[5])));
            // not loading 6:favorableAlleleSource
            if (fields.length>7) qtl.setAttribute("lod", String.valueOf(doubleOrZero(fields[7])));
            if (fields.length>8) qtl.setAttribute("likelihoodRatio", String.valueOf(doubleOrZero(fields[8])));
            if (fields.length>9) qtl.setAttribute("markerR2", String.valueOf(doubleOrZero(fields[9])));
            // not loading 10:totalR2
            // not loading 11:additivity
        }
        br.close();
    }

    /**
     * Process a trait.tsv file, creating Trait records.
     * Trait name is likely not unique across mine, so primaryIdentifier has prefix collection identifier prefix.
     * 0                           2
     * #trait_name                 description/method/etc.
     * Early leaf spot resistance  Leafs were photographed and spots were counted...
     */
    void processTraitFile(Reader reader) throws IOException {
        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
            // 0:Trait
            String name = fields[0];
            Item trait = getTrait(name);
            // 1:Trait.description
            if (fields.length>1) {
                String description = fields[1];
                if (description!=null && description.length()>0) trait.setAttribute("description", description);
            }
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
        BufferedReader br = new BufferedReader(reader);
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
            // 0:Trait
            String traitName = fields[0];
            Item trait = getTrait(traitName);
            // 1:OntologyTerm
            String ontologyId = fields[1];
            Item ontologyTerm = getOntologyTerm(ontologyId);
            // OntologyAnnotation
            Item ontologyAnnotation = createItem("OntologyAnnotation");
            ontologyAnnotation.setReference("subject", trait);
            ontologyAnnotation.setReference("ontologyTerm", ontologyTerm);
            ontologyAnnotations.add(ontologyAnnotation);
        }
        br.close();
    }

    /**
     * Process a GWASResult file with trait-marker associations.
     * Prefix trait_name with collection identifier for primaryIdentifier.
     * 0                                        1               2
     * #trait_name				marker		pvalue
     * 100 Seed weight from Florida-7 NAM	Affx-152042939	9.12e-9
     * 100 Pod weight from Florida-7 NAM	Affx-152042939	9.12e-9
     */
    void processResultFile(Reader reader) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(reader);
	String line;
        while ((line=bufferedReader.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
            if (fields.length<3) {
                throw new RuntimeException("File "+getCurrentFile().getName()+" has less than three fields in line:"+line);
            }
            // 0:Trait
            String traitName = fields[0];
            Item trait = getTrait(traitName);
            // 1-2:GWASResult
            String markerName = fields[1];
            double pValue = Double.parseDouble(fields[2]);
            Item gwasResult = createItem("GWASResult");
            gwasResult.setAttribute("markerName", markerName);
            gwasResult.setAttribute("pValue", String.valueOf(pValue));
            gwasResult.setReference("trait", trait);
            gwasResults.add(gwasResult);
            trait.addToCollection("gwasResults", gwasResult);
        }
        bufferedReader.close();
    }

    /**
     * Return a new or existing QTL Item keyed by name
     */
    Item getQTL(String name) {
        if (qtls.containsKey(name)) {
            return qtls.get(name);
        } else {
            Item qtl = createItem("QTL");
            qtl.setAttribute("name", name);
            qtl.setAttribute("primaryIdentifier", getCollectionIdentifier()+":"+name);
            qtls.put(name, qtl);
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
            linkageGroup.setAttribute("primaryIdentifier", identifier);
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
            trait.setAttribute("primaryIdentifier", getCollectionIdentifier()+":"+name);
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
     * └── vigun.MAGIC-2017.gen.Huynh_Ehlers_2018.qtl.tsv
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
