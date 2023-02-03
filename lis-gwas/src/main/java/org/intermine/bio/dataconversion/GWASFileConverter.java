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
import org.ncgr.datastore.validation.GWASCollectionValidator;
import org.ncgr.zip.GZIPBufferedReader;

/**
 * Store GWAS data from tab-delimited files in LIS Datastore /gwas/ collections.
 *
 * Genetic markers are identified only by name, there is no reference to a specific genome mapping,
 * so just the names are stored and a post-processor creates the GeneticMarker objects for known mappings.
 *
 * @author Sam Hokin
 */
public class GWASFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GWASFileConverter.class);

    // local Items to store
    Item gwas;
    List<Item> gwasResults = new ArrayList<>();
    List<Item> ontologyAnnotations = new ArrayList<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> traits = new HashMap<>();

    // utility
    String collectionIdentifier;
    
    // validate the collection first by storing a flag
    boolean collectionValidated = false;

    /**
     * Create a new GWASFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GWASFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     *
     * MAGIC-2017.gwas.Huynh_Ehlers_2018/
     * ├── README.MAGIC-2017.gwas.Huynh_Ehlers_2018.yml
     * ├── vigun.MAGIC-2017.gwas.Huynh_Ehlers_2018.obo.tsv
     * ├── vigun.MAGIC-2017.gwas.Huynh_Ehlers_2018.trait.tsv
     * └── vigun.MAGIC-2017.gwas.Huynh_Ehlers_2018.result.tsv
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (!collectionValidated) {
            GWASCollectionValidator validator = new GWASCollectionValidator(getCurrentFile().getParent());
            validator.validate();
            if (!validator.isValid()) {
                throw new RuntimeException("Collection "+getCurrentFile().getParent()+" does not pass validation.");
            }
            collectionValidated = true;
        }
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
            collectionIdentifier = readme.identifier;
            gwas = createItem("GWAS");
            gwas.setReference("organism", organism);
            gwas.setAttribute("primaryIdentifier", collectionIdentifier);
            gwas.setAttribute("synopsis", readme.synopsis);
            gwas.setAttribute("description", readme.description);
            gwas.setAttribute("genotypingPlatform", readme.genotyping_platform);
            if (readme.genotyping_method!=null) gwas.setAttribute("genotypingMethod", readme.genotyping_method);
            // store |-delimited list of genotypes
            String genotypes = "";
            for (String genotype : readme.genotype) {
                if (genotypes.length()>0) genotypes += "|";
                genotypes += genotype;
            }
            gwas.setAttribute("genotypes", genotypes);
        }
        if (getCurrentFile().getName().endsWith("obo.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processOboFile();
	} else if (getCurrentFile().getName().endsWith("trait.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processTraitFile();
	} else if (getCurrentFile().getName().endsWith("result.tsv.gz")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processResultFile();
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
        gwas.addToCollection("publications", publication);
        // associate GWAS with results (in case README not read first)
        for (Item gwasResult : gwasResults) {
            gwasResult.setReference("gwas", gwas);
        }
        // associate GWAS with traits (in case README not read first)
        for (Item trait : traits.values()) {
            trait.setReference("gwas", gwas);
        }
        // store 'em
        store(gwas);
        store(gwasResults);
        store(traits.values());
        store(ontologyTerms.values());
        store(ontologyAnnotations);
    }
    
    /**
     * Process a trait.tsv file, creating Trait records. Both columns are required.
     * Trait name is likely not unique across mine, so primaryIdentifier has collection identifier prefix.
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
            // add ontology annotation
            Item ontologyAnnotation = createItem("OntologyAnnotation");
            ontologyAnnotation.setReference("subject", trait);
            ontologyAnnotation.setReference("ontologyTerm", ontologyTerm);
            ontologyAnnotations.add(ontologyAnnotation);
        }
        br.close();
    }

    /**
     * Process a GWASResult file with trait-marker associations. All three fields are required.
     * Prefix trait_name with collection identifier for primaryIdentifier.
     * 0                                        1               2
     * #trait_name				marker		pvalue
     * 100 Seed weight from Florida-7 NAM	Affx-152042939	9.12e-9
     * 100 Pod weight from Florida-7 NAM	Affx-152042939	9.12e-9
     */
    void processResultFile() throws IOException {
        BufferedReader bufferedReader = GZIPBufferedReader.getReader(getCurrentFile());
	String line;
        while ((line=bufferedReader.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
            if (fields.length<3) {
                throw new RuntimeException("File "+getCurrentFile().getName()+" has less than three fields in line:"+line);
            }
            // required fields
            String traitName = fields[0];
            String markerName = fields[1];
            String pValueString = fields[2];
            // throw Exception if p-value is not double-parseable
            double pValue = 0.0;
            try {
                pValue = Double.parseDouble(pValueString);
            } catch (NumberFormatException ex) {
                throw new RuntimeException("File "+getCurrentFile().getName()+" has a malformatted p-value:"+pValueString+" in line:"+line);
            }                
            // 0:Trait.name
            Item trait = getTrait(traitName);
            // 0,1,2:GWASResult
            Item gwasResult = getGWASResult(trait, markerName, pValue);
            trait.addToCollection("gwasResults", gwasResult);
        }
        bufferedReader.close();
    }

    /**
     * Return a new or existing Trait Item keyed by name; primaryIdentifier is concocted from collection identifier and name.
     */
    Item getTrait(String name) {
        if (traits.containsKey(name)) {
            return traits.get(name);
        } else {
            Item trait = createItem("Trait");
            trait.setAttribute("name", name);
            trait.setAttribute("primaryIdentifier", collectionIdentifier+":"+name);
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
     * Return a GWASResult with the given trait, marker name, and p-value.
     */
    Item getGWASResult(Item trait, String markerName, double pValue) {
        Item gwasResult = createItem("GWASResult");
        gwasResult.setReference("trait", trait);
        gwasResult.setAttribute("markerName", markerName);
        gwasResult.setAttribute("pValue", String.valueOf(pValue));
        gwasResults.add(gwasResult);
        return gwasResult;
    }

    /**
     * Return the collection identifier portion of the current file.
     *     0     1          2   3                 4   5
     * └── vigun.MAGIC-2017.gwas.Huynh_Ehlers_2018.qtl.tsv
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
