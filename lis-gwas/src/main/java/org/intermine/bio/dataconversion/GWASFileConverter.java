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

/**
 * Store GWAS/phenotype/marker data from a tab-delimited file.
 * 
 * glyma.mixed.gwas1.1W14.KGK20170714-1.gwas.tsv
 * ---------------------------------------------
 * TaxonID	3847
 * Identifier	KGK20170714.1
 * Name	Bandillo, Jarquin et al. 2015
 * PlatformName	SoySNP50k
 * PlatformDetails	Illumina Infinium Bead Chip
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
    List<Item> publications = new ArrayList<>();
    List<Item> gwases = new ArrayList<>();
    List<Item> gwasResults = new ArrayList<>();
    Map<String,Item> phenotypeMap = new HashMap<>();
    Map<String,Item> markerMap = new HashMap<>();

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
     * Process the marker-chromosome relationships by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (!getCurrentFile().getName().endsWith(".tsv")) return;
        LOG.info("Processing file "+getCurrentFile().getName()+"...");
	Item dataSet = getDataSet();
	Item organism = getOrganism();
	
	// header items
        Item gwas = createItem("GWAS");
	gwas.setReference("organism", organism);
	gwas.setReference("dataSet", dataSet);
	
        Item publication = createItem("Publication");

        BufferedReader bufferedReader = new BufferedReader(reader);
	String line;
        while ((line=bufferedReader.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] parts = line.split("\t");
	    if (parts.length<2) continue; // entry without value
            String key = parts[0];
            String value = parts[1];
	    if (value.length()==0) continue; // entry without value

	    // TaxonID      3847
	    // Identifier KGK20170714.1
	    // Name	Bandillo, Jarquin et al. 2015
	    // PlatformName SoySNP50k
	    // PlatformDetails      Illumina Infinium Bead Chip
	    // DOI  10.3835/plantgenome2015.04.0024
            if (key.toLowerCase().equals("taxonid")) {
		// skip, we get organism from gensp in filename

            } else if (key.toLowerCase().equals("identifier")) {
                gwas.setAttribute("primaryIdentifier", value);

	    } else if (key.toLowerCase().equals("name")) {
		gwas.setAttribute("name", value);
		
            } else if (key.toLowerCase().equals("platformname")) {
                gwas.setAttribute("platformName", value);
		
            } else if (key.toLowerCase().equals("platformdetails")) {
                gwas.setAttribute("platformDetails", value);

            } else if (key.toLowerCase().equals("pmid")) {
                int pmid = Integer.parseInt(value);
                publication = createItem("Publication");
                publication.setAttribute("pubMedId", String.valueOf(pmid));
		publications.add(publication);
                LOG.info("Stored publication PMID="+pmid);
                gwas.addToCollection("publications", publication);

            } else if (key.toLowerCase().equals("doi")) {
                String doi = value;
                publication = createItem("Publication");
                publication.setAttribute("doi", doi);
		publications.add(publication);
                LOG.info("Stored publication DOI="+doi);
                gwas.addToCollection("publications", publication);

            } else {
                // process a GWASResult record //
                GWASFileRecord rec = new GWASFileRecord(line);
                Item marker = markerMap.get(rec.marker);
		if (marker==null) {
                    marker = createItem("GeneticMarker");
                    marker.setAttribute("secondaryIdentifier", rec.marker); // GWAS are genetic, so don't have full-yuck marker prefix
		    marker.setReference("organism", organism);
                    markerMap.put(rec.marker, marker);
                }
		marker.addToCollection("publications", publication);
                Item phenotype = phenotypeMap.get(rec.phenotype);
		if (phenotype==null) {
                    phenotype = createItem("Phenotype");
		    phenotype.setAttribute("primaryIdentifier", rec.phenotype);
                    phenotypeMap.put(rec.phenotype, phenotype);
                }
                phenotype.addToCollection("publications", publication);
                Item gwasResult = createItem("GWASResult");
		gwasResult.setAttribute("identifier", rec.identifier);
                gwasResult.setAttribute("pValue", String.valueOf(rec.pvalue));
                gwasResult.setReference("gwas", gwas);
                gwasResult.setReference("phenotype", phenotype);
                gwasResult.setReference("marker", marker);
		gwasResults.add(gwasResult);
            }
        }

        // finally store this GWAS experiment
	gwases.add(gwas);
        
	// and close its file
        bufferedReader.close();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
	store(dataSource);
	store(dataSets.values());
	store(organisms.values());
	store(publications);
	store(gwases);
	store(gwasResults);
	store(markerMap.values());
        store(phenotypeMap.values());
    }
}
