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
import java.util.ArrayList;
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
 * #Phenotype	Identifier
 * Hastate leaf shape	Hls
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

    // things to store
    Map<String,Item> qtlMap = new HashMap<>();
    Map<String,Item> experimentMap = new HashMap<>();
    Map<String,Item> phenotypeMap = new HashMap<>();
    Map<String,Item> markerMap = new HashMap<>();
    Map<String,Item> publicationMap = new HashMap<>();
    Set<Item> qtlMarkers = new HashSet<>();
    
    /**
     * Create a new QTLFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public QTLFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
	dataSource = getDataSource();
    }

    /**
     * {@inheritDoc}
     * Process files in a QTL directory
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().endsWith("expt.tsv")) {
	    processExperiment(reader);
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
	Item publication = null;
	experiment.addToCollection("dataSets", dataSet);
	experiment.setReference("organism", organism);
	boolean experimentHasIdentifier = false;
	boolean experimentHasQTLs = false;
        BufferedReader bufferedReader = new BufferedReader(reader);
	String line;
        while ((line=bufferedReader.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) {
                continue;
            }
            String[] parts = line.split("\t");
	    if (parts[0].equals("Identifier")) {
		experimentHasIdentifier = true;
		experiment.setAttribute("primaryIdentifier", parts[1]);
		experimentMap.put(parts[1], experiment);
	    } else if (parts[0].equals("TaxonID")) {
		// do nothing, get organism from gensp in filename
	    } else if (parts[0].equals("Name")) {
		experiment.setAttribute("name", parts[1]);
	    } else if (parts[0].equals("Description")) {
		experiment.setAttribute("description", parts[1]);
	    } else if (parts[0].equals("MappingParent")) {
		Item strain = getStrain(parts[1], organism);
		experiment.addToCollection("mappingParents", strain);
	    } else if (parts[0].equals("MappingDescription")) {
		experiment.setAttribute("mappingDescription", parts[1]);
	    } else if (parts[0].equals("GenotypingPlatform")) {
		experiment.setAttribute("genotypingPlatform", parts[1]);
	    } else if (parts[0].equals("GenotypingMethod")) {
		experiment.setAttribute("genotypingMethod", parts[1]);
	    } else if (parts[0].equals("PMID")) {
		if (publication==null) {
		    publication = publicationMap.get(parts[1]);
		    if (publication==null) {
			publication = createItem("Publication");
			publication.setAttribute("pubMedId", parts[1]);
			publicationMap.put(parts[1], publication);
		    }
		} else {
		    publication.setAttribute("pubMedId", parts[1]);
		}
	    } else if (parts[0].equals("DOI")) {
		if (publication==null) {
		    publication = publicationMap.get(parts[1]);
		    if (publication==null) {
			publication = createItem("Publication");
			publication.setAttribute("doi", parts[1]);
			publicationMap.put(parts[1], publication);
		    }
		} else {
		    publication.setAttribute("doi", parts[1]);
		}
	    } else {
		// validation
		if (!experimentHasIdentifier) throw new RuntimeException("Experiment file "+getCurrentFile().getName()+" lacks the required Identifier record.");
		if (publication==null) throw new RuntimeException("Experiment file "+getCurrentFile().getName()+" lacks a publication PMID or DOI record.");
		// phenotype-QTL record
		String phenotypeId = parts[0];
		String qtlId = parts[1];
		Item phenotype = phenotypeMap.get(phenotypeId);
		if (phenotype==null) {
		    phenotype = createItem("Phenotype");
		    phenotype.setAttribute("primaryIdentifier", phenotypeId);
		    phenotypeMap.put(phenotypeId, phenotype);
		}
		phenotype.addToCollection("dataSets", dataSet);
		Item qtl = qtlMap.get(qtlId);
		if (qtl==null) {
		    experimentHasQTLs = true;
		    qtl = createItem("QTL");
		    qtl.setReference("organism", organism);
		    qtl.setAttribute("primaryIdentifier", qtlId);
		    qtl.setReference("experiment", experiment);
		    if (publication!=null) qtl.addToCollection("publications", publication);
		    qtlMap.put(qtlId, qtl);
		}
		qtl.addToCollection("dataSets", dataSet);
		qtl.addToCollection("phenotypes", phenotype);
	    }
        }
	// validation
	if (!experimentHasQTLs) throw new RuntimeException("Experiment file "+getCurrentFile().getName()+" lacks QTL records.");
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
            if (line.startsWith("#") || line.trim().length()==0) {
                continue;
            }
            String[] parts = line.split("\t");
	    String qtlId = parts[0];
	    String markerId = parts[1];
	    String distinction = null;
	    if (parts.length>2) distinction = parts[2];
	    Item qtl = qtlMap.get(qtlId);
	    if (qtl==null) throw new RuntimeException("QTL "+qtlId+" is present in "+getCurrentFile().getName()+" but is missing from corresponding expt file.");
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
	store(experimentMap.values());
	store(qtlMap.values());
	store(phenotypeMap.values());
	store(markerMap.values());
	store(publicationMap.values());
	store(qtlMarkers);
    }
}
