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
 * Store QTLs and QTL-marker relationships for a GeneticMap.
 *
 * @author Sam Hokin
 */
public class QTLFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(QTLFileConverter.class);
    
    // local things to store in close()
    List<Item> ontologyAnnotations = new LinkedList<>();
    List<Item> qtlMarkers = new LinkedList<>();
    Map<String,Item> traits = new HashMap<>();
    Map<String,Item> qtls = new HashMap<>();
    Map<String,Item> linkageGroups = new HashMap<>();
    Map<String,Item> markers = new HashMap<>();
    Map<String,Item> ontologyTerms = new HashMap<>();

    // local singletons
    Item geneticMap = createItem("GeneticMap");
    Item publication = createItem("Publication");

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
        if (getCurrentFile().getName().startsWith("README")) {
	    processReadme(reader);
        } else if (getCurrentFile().getName().endsWith("obo.tsv")) {
            processOboFile(reader);
	} else if (getCurrentFile().getName().endsWith("mrk.tsv")) {
	    processMrkFile(reader);
	} else if (getCurrentFile().getName().endsWith("qtl.tsv")) {
            processQTLFile(reader);
	} else if (getCurrentFile().getName().endsWith("trait.tsv")) {
            processTraitFile(reader);
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
        store(geneticMap);
        store(publication);
        store(linkageGroups.values());
        store(traits.values());
	store(qtls.values());
	store(markers.values());
        store(qtlMarkers);
        store(ontologyAnnotations);
        store(ontologyTerms.values());
    }

    /**
     * Process the README, which contains the GeneticMap metadata.
     *
     * README.TT_Tifrunner_x_GT-C20_c.yml
     * ----------------------------------
     * identifier: TT_Tifrunner_x_GT-C20_c
     * synopsis: "Genetic map of Tifrunner x GT-C20 for the study of early and late leaf spots (ELS and LLS) and tomato spotted wilt virus (TSWV)."
     * taxid: 3818
     * genotype:
     * - Tifrunner
     * - GT-C20
     * description: "Leaf spots, including early leaf spot (ELS) and late leaf spot (LLS), and Tomato spotted wilt virus (TSWV) are devastating diseases..."
     * publication_doi: 10.1111/pbi.12930
     * publication_title: "High-Density Genetic Map Using Whole-Genome..."
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
            readme.publication_title==null) {
            throw new RuntimeException("ERROR: a required field is missing from "+getCurrentFile().getName()+": "+
                                       "Required fields are: identifier, taxid, synopsis, description, genotype, publication_doi, publication_title");
        }
        // Organism from README taxid rather than filename
        Item organism = getOrganism(Integer.parseInt(readme.taxid));
        // GeneticMap
        geneticMap.setReference("organism", organism);
        geneticMap.setAttribute("primaryIdentifier", readme.identifier);
        geneticMap.setAttribute("synopsis", readme.synopsis);
        geneticMap.setAttribute("experimentDescription", readme.description);
        // Strain = mapping parents
        for (String mappingParent : readme.genotype) {
            Item strain = getStrain(mappingParent, organism);
            geneticMap.addToCollection("mappingParents", strain);
        }
        // Publication
        publication.setAttribute("doi", readme.publication_doi);
        publication.setAttribute("title", readme.publication_title);
        geneticMap.addToCollection("publications", publication);
        // override DataSet.description from README
        Item dataSet = getDataSet();
        dataSet.setAttribute("description", readme.description);
        dataSet.setReference("publication", publication);
        geneticMap.addToCollection("dataSets", dataSet);
    }

    /**
     * Process a trait.tsv file, creating Trait records.
     * 0                1                           2
     * #trait_id        trait_name                  description/method/etc.
     * Early leaf spot  Early leaf spot resistance  Traits bearing on the reaction of the plant...
     */
    void processTraitFile(Reader reader) throws IOException {
        Item dataSet = getDataSet();
        geneticMap.addToCollection("dataSets", dataSet);
        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
            String identifier = fields[0];
            String name = fields[1];
            String description = fields[2];
            // Trait
            Item trait = traits.get(identifier);
            if (trait==null) {
                trait = createItem("Trait");
                trait.setAttribute("primaryIdentifier", identifier);
                traits.put(identifier, trait);
            }
            trait.setAttribute("name", name);
            trait.setAttribute("description", description);
            trait.setReference("geneticMap", geneticMap);
            trait.addToCollection("dataSets", dataSet);
        }
        br.close();
    }

    /**
     * Process a qtl.tsv file. The first five columns are required.
     * 0                    1                 2                            3         4          5         6                        7    8                 9          10        11
     * #qtl_name            trait_id          lg                           left_end  right_end  qtl_peak  favorable_allele_source  lod  likelihood_ratio  marker_r2  total_r2  additivity
     * Early leaf spot 1-1  Early leaf spot   TT_Tifrunner_x_GT-C20_c-A08  100.7     102.9      102                                3.02 12.42             0.56                             
     */
    void processQTLFile(Reader reader) throws IOException {
        Item dataSet = getDataSet();
        Item organism = getOrganism();
        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank line
            String[] fields = line.split("\t");
            // required columns
            String qtlId = fields[0];
            String traitId = fields[1];
            String lgId = fields[2];
            double leftEnd = doubleOrZero(fields[3]);
            double rightEnd = doubleOrZero(fields[4]);
            // optional columns
            double peak = Double.MAX_VALUE;
            String favorableAlleleSource = null;
            double lod = Double.MAX_VALUE;
            double likelihoodRatio = Double.MAX_VALUE;
            double markerR2 = Double.MAX_VALUE;
            double totalR2 = Double.MAX_VALUE;
            double additivity = Double.MAX_VALUE;
            if (fields.length>5) peak = doubleOrZero(fields[5]);
            if (fields.length>6) favorableAlleleSource = fields[6];
            if (fields.length>7) lod = doubleOrZero(fields[7]);
            if (fields.length>8) likelihoodRatio = doubleOrZero(fields[8]);
            if (fields.length>9) markerR2 = doubleOrZero(fields[9]);
            if (fields.length>10) totalR2 = doubleOrZero(fields[10]);
            if (fields.length>11) additivity = doubleOrZero(fields[11]);
            // LinkageGroup
            Item linkageGroup = linkageGroups.get(lgId);
            if (linkageGroup==null) {
                linkageGroup = createItem("LinkageGroup");
                linkageGroup.setAttribute("identifier", lgId);
                linkageGroups.put(lgId, linkageGroup);
            }
            linkageGroup.setReference("organism", organism);
            linkageGroup.addToCollection("dataSets", dataSet);
            linkageGroup.setReference("geneticMap", geneticMap);
            // Trait
            Item trait = traits.get(traitId);
            if (trait==null) {
                trait = createItem("Trait");
                trait.setAttribute("primaryIdentifier", traitId);
                traits.put(traitId, trait);
            }
            trait.setReference("geneticMap", geneticMap);
            trait.addToCollection("dataSets", dataSet);
            // QTL
            Item qtl = qtls.get(qtlId);
            if (qtl==null) {
                qtl = createItem("QTL");
                qtl.setAttribute("identifier", qtlId);
                qtls.put(qtlId, qtl);
            }
            qtl.setReference("organism", organism);
            qtl.setReference("geneticMap", geneticMap);
            qtl.addToCollection("dataSets", dataSet);
            // required parameters
            qtl.setReference("trait", trait);
            qtl.setReference("linkageGroup", linkageGroup);
            qtl.setAttribute("start", String.valueOf(leftEnd));
            qtl.setAttribute("end", String.valueOf(rightEnd));
            // optional parameters
            if (peak<Double.MAX_VALUE) qtl.setAttribute("peak", String.valueOf(peak));
            if (lod<Double.MAX_VALUE) qtl.setAttribute("lod", String.valueOf(lod));
            if (likelihoodRatio<Double.MAX_VALUE) qtl.setAttribute("likelihoodRatio", String.valueOf(likelihoodRatio));
            if (markerR2<Double.MAX_VALUE) qtl.setAttribute("markerR2", String.valueOf(markerR2));
            // favorableAlleleSource 
            // totalR2
            // additivity = Double.MAX_VALUE;
        }
        br.close();
    }

    /**
     * Process an obo.tsv file.
     * 0            1
     * #trait_id    obo_term
     * Seed weight  TO:0000181
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
     * Process a qtlmrk.tsv file.
     * #qtl_name            marker        distinction
     * Early leaf spot 1-1  A08_35596996  flanking
     * Early leaf spot 1-1  A08_35776787  nearest
     *
     */
    void processMrkFile(Reader reader) throws IOException {
	Item dataSet = getDataSet();
	Item organism = getOrganism();
        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
	    String qtlId = fields[0];
	    String markerId = fields[1];
	    String distinction = null;
	    if (fields.length>2) distinction = fields[2];
            // QTL
	    Item qtl = qtls.get(qtlId);
	    if (qtl==null) {
                qtl = createItem("QTL");
                qtl.setReference("organism", organism);
                qtl.setAttribute("identifier", qtlId);
                qtls.put(qtlId, qtl);
            }
	    qtl.addToCollection("dataSets", dataSet);
            // GeneticMarker
	    Item marker = markers.get(markerId);
	    if (marker==null) {
		marker = createItem("GeneticMarker");
		marker.setReference("organism", organism);
		marker.setAttribute("secondaryIdentifier", markerId);
		markers.put(markerId, marker);
	    }
	    marker.addToCollection("dataSets", dataSet);
            // QTLMarker
	    Item qtlMarker = createItem("QTLMarker");
	    if (distinction!=null) qtlMarker.setAttribute("distinction", distinction);
	    qtlMarker.setReference("qtl", qtl);
	    qtlMarker.setReference("marker", marker);
            qtlMarkers.add(qtlMarker);
	}
        br.close();
    }

    static double doubleOrZero(String field) {
        if (field.length()>0) {
            return Double.parseDouble(field);
        } else {
            return 0.00;
        }
    }
}
