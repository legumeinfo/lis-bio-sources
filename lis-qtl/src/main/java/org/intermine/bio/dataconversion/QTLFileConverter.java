package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2015-2016 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

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

/**
 * Store QTL experiments and QTL-marker relationships.
 *
 * gensp.mixed.qtl.KEY4.identifier.expt.tsv
 * ----------------------------------------
 * Identifier          Experiment unique identifier, e.g. 22691139
 * Name                Experiment name, e.g. Jones, Smith, et al. 2005
 * Description         Description of the experiment; often from the pub abstract
 * MappingParent       gensp.Accession = one of the accessions used in the map
 * MappingParent       gensp.Accession = another of the accessions used in the map
 * MappingParent       gensp.Accession = yet another accession, etc.
 * MappingDescription  Description of the development of the mapping population
 * PMID                PubMed ID of the publication; one of PMID, DOI required
 * DOI                 DOI of the publication; one of PMID, DOI required
 * GenotypingPlatform  Description of the genotyping platform
 * GenotypingMethod    Description of the genotyping method
 * PlatformName        Identifies the GFF files linking markers to genomes, e.g. SoySNP50K.
 *                     This name is used in the GFF filename, following .key4.
 * Treatment           Description of treatment applied to detect specific traits (e.g. high and low salt treatments to find salt tolerance QTLs).
 * AnalysisMethod      Method used to detect QTL or marker-trait associations (e.g., SMA, IM, SIM, CIM, ICIM.)
 *
 * gensp.mixed.qtl.KEY4.identifier.markers.tsv
 * -------------------------------------------
 * #Identifier	Marker	Distinction
 * Hls	1_0083	Flanking
 * Hls	1_0349	
 * Hls	1_0417	Flanking
 * Hls	1_0910	Peak
 * Hls	1_0992	Flanking
 * Hls	1_1013	Flanking
 *
 * @author Sam Hokin
 */
public class QTLFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(QTLFileConverter.class);
    
    // local things to store
    List<Item> qtlMarkers = new LinkedList<>();
    List<Item> publications = new LinkedList<>();
    Map<String,Item> qtlMap = new HashMap<>();
    Map<String,Item> experimentMap = new HashMap<>();
    Map<String,Item> phenotypeMap = new HashMap<>();
    Map<String,Item> markerMap = new HashMap<>();
    Map<String,Item> ontologyAnnotationMap = new HashMap<>();
    Map<String,Item> ontologyTermMap = new HashMap<>();
    
    /**
     * Create a new QTLFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public QTLFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().endsWith("expt.tsv")) {
	    processExperiment(reader);
        } else if (getCurrentFile().getName().endsWith("phen.tsv")) {
            processPhenFile(reader);
	} else if (getCurrentFile().getName().endsWith("mrk.tsv")) {
	    processMarkers(reader);
	}
    }

    /**
     * Process a QTL experiment file. Data lines are all tab-separated pairs.
     */
    void processExperiment(Reader reader) throws IOException {
	Item dataSet = getDataSet();
	Item organism = getOrganism();
	Item experiment = createItem("QTLExperiment");
        experiment.setReference("organism", organism);
        experiment.addToCollection("dataSets", dataSet);
	Item publication = null; // based on PMID or DOI or both
	String line;
        BufferedReader bufferedReader = new BufferedReader(reader);
        while ((line=bufferedReader.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
            String key = fields[0].toLowerCase();
            String value = fields[1];
	    if (key.equals("identifier")) {
                // QTLExperiment.primaryIdentifier because it is Annotatable
                experiment.setAttribute("primaryIdentifier", value);
                experimentMap.put(value, experiment);
	    } else if (key.equals("name")) {
		experiment.setAttribute("name", value);
	    } else if (key.equals("description")) {
		experiment.setAttribute("description", value);
	    } else if (key.equals("mappingparent")) {
                String[] parts = value.split("\\.");
                if (parts.length!=2) {
                    throw new RuntimeException("MappingParent must have form gensp.StrainName. Aborting.");
                }
                String gensp = parts[0].trim();
                String strainName = parts[1].trim();
                experiment.addToCollection("mappingParents", getStrain(strainName, getOrganism(gensp)));
	    } else if (key.equals("mappingdescription")) {
		experiment.setAttribute("mappingDescription", value);
	    } else if (key.equals("genotypingplatform")) {
		experiment.setAttribute("genotypingPlatform", value);
	    } else if (key.equals("genotypingmethod")) {
		experiment.setAttribute("genotypingMethod", value);
            } else if (key.equals("platformname")) {
                experiment.setAttribute("platformName", value);
	    } else if (key.equals("pmid")) {
                if (publication==null) {
                    publication = createItem("Publication");
                    publications.add(publication);
                }
                publication.setAttribute("pubMedId", value);
	    } else if (key.equals("doi")) {
		if (publication==null) {
                    publication = createItem("Publication");
                    publications.add(publication);
                }
                publication.setAttribute("doi", value);
            } else if (key.equals("taxonid")) {
                // do nothing
	    } else {
		// data line
		if (publication==null) {
                    throw new RuntimeException("Experiment file "+getCurrentFile().getName()+" lacks a publication PMID or DOI record.");
                }
		// phenotype-QTL record: QTL and trait identifier to Cap lower case for consistency
		String qtlId = fields[0].toLowerCase();
		String trait = fields[1].trim().toLowerCase();
                qtlId = qtlId.substring(0,1).toUpperCase() + qtlId.substring(1);
                trait = trait.substring(0,1).toUpperCase() + trait.substring(1);
		Item qtl = qtlMap.get(qtlId);
		if (qtl==null) {
		    qtl = createItem("QTL");
                    qtl.setAttribute("identifier", qtlId);
		    qtlMap.put(qtlId, qtl);
                }
                qtl.setReference("organism", organism);
                qtl.setReference("experiment", experiment);
		qtl.addToCollection("dataSets", dataSet);
		Item phenotype = phenotypeMap.get(trait);
		if (phenotype==null) {
		    phenotype = createItem("Phenotype");
                    phenotype.setAttribute("primaryIdentifier", trait);
		    phenotypeMap.put(trait, phenotype);
		}
		phenotype.addToCollection("dataSets", dataSet);
		qtl.setReference("phenotype", phenotype);
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
            String trait = parts[0].trim().toLowerCase();
            String ontologyId = parts[1];
            // initcap trait
            trait = trait.substring(0,1).toUpperCase() + trait.substring(1);
            if (ontologyId==null || ontologyId.trim().length()==0) continue; // placeholder
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
     * Process a QTL-markers file. Data lines are tab-separated fields, some empty.
     *
     * #Identifier       Marker       Distinction
     * Pod dehiscence-1  ss715639553  flanking
     * Pod dehiscence-1  ss715639323  flanking
     */
    void processMarkers(Reader reader) throws IOException {
	Item dataSet = getDataSet();
	Item organism = getOrganism();
        BufferedReader bufferedReader = new BufferedReader(reader);
	String line;
        while ((line=bufferedReader.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
	    String qtlId = fields[0];
	    String markerId = fields[1];
	    String distinction = null;
	    if (fields.length>2) distinction = fields[2];
	    Item qtl = qtlMap.get(qtlId);
	    if (qtl==null) {
                qtl = createItem("QTL");
                qtl.setAttribute("identifier", qtlId);
                qtlMap.put(qtlId, qtl);
            }
	    qtl.addToCollection("dataSets", dataSet);
	    Item marker = markerMap.get(markerId);
	    if (marker==null) {
		marker = createItem("GeneticMarker");
		marker.setReference("organism", organism);
		marker.setAttribute("secondaryIdentifier", markerId);
		markerMap.put(markerId, marker);
	    }
	    marker.addToCollection("dataSets", dataSet);
	    Item qtlMarker = createItem("QTLMarker");
	    if (distinction!=null) qtlMarker.setAttribute("distinction", distinction);
	    qtlMarker.setReference("qtl", qtl);
	    qtlMarker.setReference("marker", marker);
	    qtlMarkers.add(qtlMarker);
	}
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
	// from DataStoreFileConverter
	store(dataSource);
	store(dataSets.values());
	store(organisms.values());
	store(strains.values());
	// local
	store(publications);
	store(qtlMarkers);
	store(experimentMap.values());
	store(qtlMap.values());
	store(phenotypeMap.values());
	store(markerMap.values());
        store(ontologyAnnotationMap.values());
        store(ontologyTermMap.values());
    }
}
