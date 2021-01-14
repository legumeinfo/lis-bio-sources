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
 * Store GWAS/phenotype/marker data from a tab-delimited file.
 * 
 * glyma.mixed.gwas1.1W14.KGK20170714-1.gwas.tsv
 * ---------------------------------------------
 * Identifier	KGK20170714.1
 * Name	Bandillo, Jarquin et al. 2015
 * Description  We did a GWAS on 38 accessions of soybean and here are the results.
 * PlatformName	SoySNP50k
 * PlatformDetails   Illumina Infinium Bead Chip
 * GenotypeDataset   LIS DS folder containing a VCF or HMP file underlying analysis
 * PhenotypeDataset  LIS DS folder containing a trait measurement file underlying analysis
 * MarkerSet         SoyBead
 * PMID 123456
 * DOI	10.3835/plantgenome2015.04.0024
 * #identifier	   phenotype  marker	      pvalue
 * Seed oil 4-g14  Seed oil   ss715591641   3.16E-09
 * Seed oil 4-g15  Seed oil   ss715591642   3.16E-08
 * ...
 *
 * @author Sam Hokin
 */
public class GWASFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GWASFileConverter.class);

    // local things to store
    List<Item> gwases = new LinkedList<>();
    List<Item> gwasResults = new LinkedList<>();
    Map<Integer,Item> publicationPMIDMap = new HashMap<>(); // keyed by PMID if supplied
    Map<String,Item> publicationDOIMap = new HashMap<>();   // keyed by DOI if supplied and PMID isn't
    Map<String,Item> ontologyAnnotationMap = new HashMap<>();
    Map<String,Item> phenotypeMap = new HashMap<>();
    Map<String,Item> markerMap = new HashMap<>(); // keyed by secondaryIdentifier
    Map<String,Item> ontologyTermMap = new HashMap<>();

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
        if (!getCurrentFile().getName().endsWith(".gwas.tsv") && !getCurrentFile().getName().endsWith(".phen.tsv")) return;
        if (getCurrentFile().getName().endsWith(".gwas.tsv")) {
            processGWASFile(reader);
        } else if (getCurrentFile().getName().endsWith(".phen.tsv")) {
            processPhenFile(reader);
        }
    }

    /**
     * Process a GWAS file
     */
    void processGWASFile(Reader reader) throws IOException {
        Item dataSet = getDataSet();
        Item organism = getOrganism();
        Item gwas = createItem("GWAS");
	gwas.setReference("dataSet", dataSet);
        gwas.setReference("organism", organism);
        gwases.add(gwas);
        Item publication = null;
        BufferedReader bufferedReader = new BufferedReader(reader);
	String line;
        while ((line=bufferedReader.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] parts = line.split("\t");
	    if (parts.length<2) continue; // entry without value
            String key = parts[0].toLowerCase();
            String value = parts[1];
	    if (value.length()==0) continue; // entry without value
            if (key.equals("identifier")) {
                gwas.setAttribute("primaryIdentifier", value);
	    } else if (key.equals("name")) {
		gwas.setAttribute("name", value);
            } else if (key.equals("description")) {
                gwas.setAttribute("description", value);
            } else if (key.equals("platformname")) {
                gwas.setAttribute("platformName", value);
            } else if (key.equals("platformdetails")) {
                gwas.setAttribute("platformDetails", value);
            } else if (key.equals("markerset")) {
                // gwas.setAttribute("markerSet", value);
            } else if (key.equals("pmid")) {
                int pmid = Integer.parseInt(value); // make sure it's a number
                if (publication==null) {
                    if (publicationPMIDMap.containsKey(pmid)) {
                        publication = publicationPMIDMap.get(pmid);
                    } else {
                        publication = createItem("Publication");
                        publication.setAttribute("pubMedId", String.valueOf(pmid));
                        publicationPMIDMap.put(pmid, publication);
                    }
                    gwas.addToCollection("publications", publication);
                } else {
                    // update DOI pub
                    publication.setAttribute("pubMedId", String.valueOf(pmid));
                }                    
            } else if (key.equals("doi")) {
                String doi = value;
                if (publication==null) {
                    if (publicationDOIMap.containsKey(doi)) {
                        publication = publicationDOIMap.get(doi);
                    } else {
                        publication = createItem("Publication");
                        publication.setAttribute("doi", doi);
                        publicationDOIMap.put(doi, publication);
                    }
                    gwas.addToCollection("publications", publication);
                } else {
                    // update PMID pub
                    publication.setAttribute("doi", doi);
                }
            } else {
                // data record
                GWASFileRecord rec = new GWASFileRecord(line);
                // GeneticMarker
                String secondaryIdentifier = rec.marker;
                Item marker = markerMap.get(secondaryIdentifier);
		if (marker==null) {
                    marker = createItem("GeneticMarker");
                    marker.setAttribute("secondaryIdentifier", secondaryIdentifier);
                    markerMap.put(secondaryIdentifier, marker);
                }
                marker.setReference("organism", organism);
                marker.addToCollection("dataSets", dataSet);
                if (publication!=null) {
                    marker.addToCollection("publications", publication);
                }
                // Phenotype
                Item phenotype = phenotypeMap.get(rec.phenotype);
		if (phenotype==null) {
                    phenotype = createItem("Phenotype");
		    phenotype.setAttribute("primaryIdentifier", rec.phenotype);
                    phenotypeMap.put(rec.phenotype, phenotype);
                }
                phenotype.setAttribute("secondaryIdentifier", rec.phenotype);
                phenotype.setReference("organism", organism);
                // phenotype.addToCollection("publications", publication);
                phenotype.addToCollection("dataSets", dataSet);
                // GWASResult
                Item gwasResult = createItem("GWASResult");
		gwasResult.setAttribute("identifier", rec.identifier);
                gwasResult.setAttribute("pValue", String.valueOf(rec.pvalue));
                if (rec.lod!=Double.MAX_VALUE) {
                    gwasResult.setAttribute("lod", String.valueOf(rec.lod));
                }
                gwasResult.setReference("gwas", gwas);
                gwasResult.setReference("marker", marker);
                gwasResult.setReference("phenotype", phenotype);
		gwasResults.add(gwasResult);
            }
        }
        bufferedReader.close();
    }

    /**
     * Process a Phenotype file
     */
    void processPhenFile(Reader reader) throws IOException {
        Item dataSet = getDataSet();
        Item organism = getOrganism();
        BufferedReader bufferedReader = new BufferedReader(reader);
	String line;
        while ((line=bufferedReader.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] parts = line.split("\t");
	    if (parts.length<2) continue; // entry without value
            String trait = parts[0];
            String ontologyId = parts[1];
            // Phenotype
            Item phenotype = phenotypeMap.get(trait);
            if (phenotype==null) {
                phenotype = createItem("Phenotype");
                phenotype.setAttribute("primaryIdentifier", trait);
                phenotypeMap.put(trait, phenotype);
            }
            phenotype.setAttribute("secondaryIdentifier", trait);
            phenotype.setReference("organism", organism);
            phenotype.addToCollection("dataSets", dataSet);
            // OntologyTerm
            Item ontologyTerm = ontologyTermMap.get(ontologyId);
            if (ontologyTerm==null) {
                ontologyTerm = createItem("OntologyTerm");
                ontologyTerm.setAttribute("identifier", ontologyId);
                ontologyTermMap.put(ontologyId, ontologyTerm);
            }
            // OntologyAnnotation
            String ontologyAnnotationKey = trait+"|"+ontologyId;
            if (!ontologyAnnotationMap.containsKey(ontologyAnnotationKey)) {
                Item ontologyAnnotation = createItem("OntologyAnnotation");
                ontologyAnnotation.setReference("subject", phenotype);
                ontologyAnnotation.setReference("ontologyTerm", ontologyTerm);
                ontologyAnnotation.addToCollection("dataSets", dataSet);
                ontologyAnnotationMap.put(ontologyAnnotationKey,ontologyAnnotation);
            }
        }
        bufferedReader.close();
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
        Set<Item> publications = new HashSet<>();
        publications.addAll(publicationPMIDMap.values());
        publications.addAll(publicationDOIMap.values());
	store(publications);
	store(gwases);
	store(gwasResults);
        store(ontologyAnnotationMap.values());
	store(markerMap.values());
        store(phenotypeMap.values());
        store(ontologyTermMap.values());
    }
}
