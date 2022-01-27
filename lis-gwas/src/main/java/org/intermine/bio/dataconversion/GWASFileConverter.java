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

/**
 * Store GWAS result/trait/marker data from a tab-delimited files.
 *
 * @author Sam Hokin
 */
public class GWASFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GWASFileConverter.class);

    // local things to store
    Item gwas = createItem("GWAS");
    List<Item> gwasResults = new LinkedList<>();
    List<Item> ontologyAnnotations = new LinkedList<>();
    Map<String,Item> traits = new HashMap<>();           // keyed by identifier
    Map<String,Item> markers = new HashMap<>();          // keyed by secondaryIdentifier
    Map<String,Item> ontologyTerms = new HashMap<>();    // keyed by identifier

    /**
     * Create a new GWASFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GWASFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().startsWith("README")) {
	    processReadme(reader);
            if (readme.genotyping_platform==null) {
                throw new RuntimeException("ERROR: "+getCurrentFile().getName()+" must contain genotyping_platform.");
            }
            // GWAS
            gwas = createItem("GWAS");
            gwas.setAttribute("primaryIdentifier", readme.identifier);
            gwas.setAttribute("synopsis", readme.synopsis);
            gwas.setAttribute("description", readme.description);
            gwas.setAttribute("genotypingPlatform", readme.genotyping_platform);
            if (readme.genotyping_method!=null) gwas.setAttribute("genotypingMethod", readme.genotyping_method);
            gwas.setAttribute("population", readme.genotype[0]);
            gwas.setReference("organism", organism);
            if (publication!=null) gwas.addToCollection("publications", publication);
        } else if (getCurrentFile().getName().endsWith("obo.tsv")) {
            processOboFile(reader);
	} else if (getCurrentFile().getName().endsWith("result.tsv")) {
            processResultFile(reader);
	} else if (getCurrentFile().getName().endsWith("trait.tsv")) {
            processTraitFile(reader);
	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        // references and collections
        for (Item marker : markers.values()) {
            marker.setReference("organism", organism);
        }
        for (Item gwasResult : gwasResults) {
            gwasResult.setReference("gwas", gwas);
        }
        // collection stuff
        storeCollectionItems();
        // local items
	store(gwas);
	store(gwasResults);
        store(ontologyAnnotations);
        store(ontologyTerms.values());
        store(traits.values());
	store(markers.values());
    }

    /**
     * Process a trait.tsv file, creating Trait records with primaryIdentifer and description.
     * 0                                   1
     * #trait_id                           description
     * 100 Seed weight from Florida-7 NAM  After harvest and drying to less than 10% water content, 100 seeds were picked randomly and weighed.
     * 100 Pod weight from Florida-7 NAM   After harvest and drying to less than 10% water content, 100 pods were picked randomly and weighed.
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
     * 0                                   1
     * #trait_id                           obo_term
     * 100 Seed weight from Florida-7 NAM  TO:0000181
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
     * Process a GWASResult file with trait-marker associations.
     * 0                                        1               2
     * #trait_id				marker		pvalue
     * 100 Seed weight from Florida-7 NAM	Affx-152042939	9.12e-9
     * 100 Pod weight from Florida-7 NAM	Affx-152042939	9.12e-9
     */
    void processResultFile(Reader reader) throws IOException {
        System.out.println("Processing "+getCurrentFile().getName());
        BufferedReader bufferedReader = new BufferedReader(reader);
	String line;
        while ((line=bufferedReader.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
            String traitId = fields[0];
            String markerId = fields[1];
            double pValue = Double.parseDouble(fields[2]);
            // Trait
            Item trait = traits.get(traitId);
            if (trait==null) {
                trait = createItem("Trait");
                trait.setAttribute("primaryIdentifier", traitId);
                traits.put(traitId, trait);
            }
            // GeneticMarker
            Item marker = markers.get(markerId);
            if (marker==null) {
                marker = createItem("GeneticMarker");
                marker.setAttribute("secondaryIdentifier", markerId);
                markers.put(markerId, marker);
            }
            // GWASResult
            Item gwasResult = createItem("GWASResult");
            gwasResult.setReference("trait", trait);
            gwasResult.setReference("marker", marker);
            gwasResult.setAttribute("pValue", String.valueOf(pValue));
            gwasResults.add(gwasResult);
            trait.addToCollection("gwasResults", gwasResult);
        }
        bufferedReader.close();
    }
}
