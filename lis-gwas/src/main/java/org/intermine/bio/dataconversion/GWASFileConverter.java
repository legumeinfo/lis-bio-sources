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

/**
 * Store GWAS result/trait/marker data from a tab-delimited files.
 *
 * @author Sam Hokin
 */
public class GWASFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GWASFileConverter.class);

    // local things to store
    List<Item> gwasResults = new LinkedList<>();
    List<Item> ontologyAnnotations = new LinkedList<>();
    Map<String,Item> traits = new HashMap<>();           // keyed by identifier
    Map<String,Item> markers = new HashMap<>();          // keyed by secondaryIdentifier
    Map<String,Item> ontologyTerms = new HashMap<>();    // keyed by identifier

    // singleton Items
    Item gwas = createItem("GWAS");
    Item publication = createItem("Publication");
                                  
    /**
     * Create a new GWASFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GWASFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
	dataSource = getDataSource();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().startsWith("README")) {
	    processReadme(reader);
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
        // DatastoreFileConverter
	store(dataSource);
	store(dataSets.values());
	store(organisms.values());
        // local
	store(publication);
	store(gwas);
	store(gwasResults);
        store(traits.values());
        store(ontologyAnnotations);
	store(markers.values());
        store(ontologyTerms.values());
    }

    /**
     * Process the README, which contains the GWAS metadata.
     *
     * identifier: 2020NAMFlor7
     * synopsis: GangurdeSS 2020 NAM Florida-7 GWAS study of pod and seed weight
     * taxid: 3818
     * genotype:
     * - Florida-07 NAM population
     * description: "A GWAS dataset from Gangurde SS, et al., 2020 NAM Florida-7 mapping population. Nested-association mapping..."
     * publication_doi: 10.1111/pbi.13311
     * publication_title: "Nested-association mapping (NAM)-based genetic dissection uncovers candidate genes for seed and pod weights in peanut (Arachis hypogaea)"
     * genotyping_platform: "Affymetrix 58K SNP Axiom_Arachis array"
     * genotyping_method: "Blah di blah di blah"
     */
    void processReadme(Reader reader) throws IOException {
        Readme readme = Readme.getReadme(reader);
        // check required stuff
        if (readme.identifier==null ||
            readme.taxid==null ||
            readme.synopsis==null ||
            readme.description==null ||
            readme.genotype==null ||
            readme.publication_doi==null ||
            readme.publication_title==null ||
            readme.genotyping_platform==null) {
            throw new RuntimeException("ERROR: a required field is missing from "+getCurrentFile().getName()+": "+
                                       "Required fields are: identifier, taxid, synopsis, description, genotype, publication_doi, publication_title, genotyping_platform");
        }
        // Organism from README taxid rather than filename
        Item organism = getOrganism(Integer.parseInt(readme.taxid));
        // GWAS
        gwas.setAttribute("primaryIdentifier", readme.identifier);
        gwas.setReference("organism", organism);
        gwas.setAttribute("synopsis", readme.synopsis);
        gwas.setAttribute("description", readme.description);
        gwas.setAttribute("genotypingPlatform", readme.genotyping_platform);
        if (readme.genotyping_method!=null) gwas.setAttribute("genotypingMethod", readme.genotyping_method);
        gwas.setAttribute("population", readme.genotype[0]);
        // Publication
        publication.setAttribute("doi", readme.publication_doi);
        publication.setAttribute("title", readme.publication_title);
        gwas.addToCollection("publications", publication);
        // override DataSet.description from README
        Item dataSet = getDataSet();
        dataSet.setAttribute("description", readme.description);
        dataSet.setReference("publication", publication);
        gwas.addToCollection("dataSets", dataSet);
    }

    /**
     * Process a trait.tsv file, creating Trait records.
     * 0                                   1            2optional
     * #trait_id                           trait_name   description
     * 100 Seed weight from Florida-7 NAM  Seed weight  After harvest and drying to less than 10% water content, 100 seeds were picked randomly and weighed.
     * 100 Pod weight from Florida-7 NAM   Pod weight   After harvest and drying to less than 10% water content, 100 pods were picked randomly and weighed.
     */
    void processTraitFile(Reader reader) throws IOException {
        Item dataSet = getDataSet();
        gwas.addToCollection("dataSets", dataSet);
        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
            // required
            String identifier = fields[0];
            String name = fields[1];
            // optional
            String description = null;
            if (fields.length>2) description = fields[2];
            // Trait
            Item trait = traits.get(identifier);
            if (trait==null) {
                trait = createItem("Trait");
                trait.setAttribute("primaryIdentifier", identifier);
                traits.put(identifier, trait);
            }
            trait.setAttribute("name", name);
            if (description!=null) trait.setAttribute("description", description);
            trait.setReference("gwas", gwas);
            trait.addToCollection("dataSets", dataSet);
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
        Item dataSet = getDataSet();
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
            trait.addToCollection("dataSets", dataSet);
            // OntologyAnnotation
            Item ontologyAnnotation = createItem("OntologyAnnotation");
            ontologyAnnotation.setReference("subject", trait);
            ontologyAnnotation.setReference("ontologyTerm", ontologyTerm);
            ontologyAnnotation.addToCollection("dataSets", dataSet);
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
        Item dataSet = getDataSet();
        Item organism = getOrganism();
        gwas.addToCollection("dataSets", dataSet);
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
            trait.setReference("gwas", gwas);
            trait.addToCollection("dataSets", dataSet);
            // GeneticMarker
            Item marker = markers.get(markerId);
            if (marker==null) {
                marker = createItem("GeneticMarker");
                marker.setAttribute("secondaryIdentifier", markerId);
                marker.setReference("organism", organism);
                markers.put(markerId, marker);
            }
            marker.addToCollection("dataSets", dataSet);
            // GWASResult
            Item gwasResult = createItem("GWASResult");
            gwasResult.setReference("gwas", gwas);
            gwasResult.setReference("trait", trait);
            gwasResult.setReference("marker", marker);
            gwasResult.setAttribute("pValue", String.valueOf(pValue));
            gwasResults.add(gwasResult);
        }
        bufferedReader.close();
    }
}
