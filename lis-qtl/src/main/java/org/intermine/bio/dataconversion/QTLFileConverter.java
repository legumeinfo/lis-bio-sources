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
 * Identifier	22691139
 * Name	Pottorff et al. 2012
 * Description	In this study, we analyzed the genetics of leaf morphology in a segregating cowpea RIL population....
 * TaxonID	3920
 * MappingParent	Sanzi
 * MappingParent	Vita7
 * MappingDescription	The mapping population consisted of 122 RILs which were advanced from the Sanzi x Vita7 cross.
 * GenotypingPlatform	Illumina GoldenGate Array
 * GenotypingMethod	The Sanzi x Vita 7 population was genotyped at the F8 generation using bi-allelic SNP markers from the 1536-marker Illumina GoldenGate Assay.
 * PMID	22691139
 * #Identifier  Trait
 * Hls          Hastate leaf shape	
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
        } else if (getCurrentFile().getName().endsWith(".phen.tsv")) {
            processPhenFile(reader);
	} else if (getCurrentFile().getName().endsWith("markers.tsv")) {
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
        BufferedReader bufferedReader = new BufferedReader(reader);
	String line;
        while ((line=bufferedReader.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
	    if (fields[0].equals("Identifier")) {
                experiment.setAttribute("primaryIdentifier", fields[1]);
                experimentMap.put(fields[1], experiment);
	    } else if (fields[0].equals("TaxonID")) {
                // ignore, we get organism from filename
	    } else if (fields[0].equals("Name")) {
		experiment.setAttribute("name", fields[1]);
	    } else if (fields[0].equals("Description")) {
		experiment.setAttribute("description", fields[1]);
	    } else if (fields[0].equals("MappingParent")) {
                String[] parts = fields[1].split("\\.");
                if (parts.length!=2) {
                    throw new RuntimeException("MappingParent must have form gensp.StrainName. Aborting.");
                }
                String gensp = parts[0].trim();
                String strainName = parts[1].trim();
                experiment.addToCollection("mappingParents", getStrain(strainName, getOrganism(gensp)));
	    } else if (fields[0].equals("MappingDescription")) {
		experiment.setAttribute("mappingDescription", fields[1]);
	    } else if (fields[0].equals("GenotypingPlatform")) {
		experiment.setAttribute("genotypingPlatform", fields[1]);
	    } else if (fields[0].equals("GenotypingMethod")) {
		experiment.setAttribute("genotypingMethod", fields[1]);
	    } else if (fields[0].equals("PMID")) {
                if (publication==null) {
                    publication = createItem("Publication");
                    publications.add(publication);
                }
                publication.setAttribute("pubMedId", fields[1]);
	    } else if (fields[0].equals("DOI")) {
		if (publication==null) {
                    publication = createItem("Publication");
                    publications.add(publication);
                }
                publication.setAttribute("doi", fields[1]);
	    } else {
		// data line
		if (publication==null) {
                    throw new RuntimeException("Experiment file "+getCurrentFile().getName()+" lacks a publication PMID or DOI record.");
                }
		// phenotype-QTL record
		String qtlId = fields[0];
		String trait = fields[1];
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
            String trait = parts[0];
            String ontologyId = parts[1];
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
	    qtlMarker.setReference("QTL", qtl);
	    qtlMarker.setReference("marker", marker);
            qtlMarker.setReference("dataSet", dataSet);
	    if (distinction!=null) qtlMarker.setAttribute("distinction", distinction);
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
