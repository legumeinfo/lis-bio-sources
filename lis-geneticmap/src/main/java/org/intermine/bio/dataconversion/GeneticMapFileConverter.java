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

/**
 * Store GeneticMap/marker/linkage group/QTL data from tab-delimited files.
 * 
 * @author Sam Hokin
 */
public class GeneticMapFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GeneticMapFileConverter.class);

    // local things to store
    List<Item> linkageGroupPositions = new LinkedList<>();
    Map<String,Item> geneticMapMap = new HashMap<>();
    Map<String,Item> publicationMap = new HashMap<>();
    Map<String,Item> linkageGroupMap = new HashMap<>();
    Map<String,Item> markerMap = new HashMap<>();
    Map<String,Item> phenotypeMap = new HashMap<>();
    Map<String,Item> qtlMap = new HashMap<>();

    /**
     * Create a new GeneticMapFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GeneticMapFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        // 0     1     2   3    4          5    6
        // vigun.mixed.map.SM2W.MAGIC-2017.expt.tsv
        // vigun.mixed.map.SM2W.MAGIC-2017.mrk.tsv
        // vigun.mixed.map.SM2W.MAGIC-2017.qtl.tsv
        if (!getCurrentFile().getName().endsWith(".tsv")) return;
        String filetype = getFileType(getCurrentFile().getName());
        if (filetype==null) return;
        if (filetype.equals("expt")) {
            processExptFile(reader);
        } else if (filetype.equals("mrk")) {
            processMrkFile(reader);
        } else if (filetype.equals("qtl")) {
            processQtlFile(reader);
        }
    }

    /**
     * Process a expt.tsv file
     */
    void processExptFile(Reader reader) throws IOException {
        Item dataSet = getDataSet();
        Item organism = getOrganism();
        Item geneticMap = null;
        Item publication = null;
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
            if (fields[0].equals("MapName")) {
                String mapname = fields[1].trim();
                geneticMap = geneticMapMap.get(mapname);
                if (geneticMap==null) {
                    geneticMap = createItem("GeneticMap");
                    geneticMap.setAttribute("primaryIdentifier", mapname);
                    geneticMap.setReference("organism", organism);
                    geneticMapMap.put(mapname, geneticMap);
                }
                geneticMap.addToCollection("dataSets", dataSet);
                publication = createItem("Publication");
                publicationMap.put(mapname, publication);
            } else if (fields[0].equals("Description")) {
                geneticMap.setAttribute("description", fields[1]);
            } else if (fields[0].equals("MappingParent")) {
                String[] parts = fields[1].split("\\.");
                if (parts.length!=2) {
                    throw new RuntimeException("MappingParent must have form gensp.StrainName. Aborting.");
                }
                String gensp = parts[0].trim();
                String strainName = parts[1].trim();
                geneticMap.addToCollection("mappingParents", getStrain(strainName, getOrganism(gensp)));
            } else if (fields[0].equals("PMID")) {
                publication.setAttribute("pubMedId", fields[1]);
            } else if (fields[0].equals("DOI")) {
                publication.setAttribute("doi", fields[1]);
            } else {
                // data line
                // 0                      1      2
                // LinkageGroupIdentifier Number Length
                String identifier = fields[0].trim();
                int number = Integer.parseInt(fields[1]);
                double length = Double.parseDouble(fields[2]);
                Item linkageGroup = linkageGroupMap.get(identifier);
                if (linkageGroup==null) {
                    linkageGroup = createItem("LinkageGroup");
                    linkageGroup.setAttribute("identifier", identifier);
                    linkageGroup.setReference("organism", organism);
                    linkageGroupMap.put(identifier, linkageGroup);
                }
                linkageGroup.addToCollection("dataSets", dataSet);
                linkageGroup.setReference("geneticMap", geneticMap);
                linkageGroup.setAttribute("number", String.valueOf(number));
                linkageGroup.setAttribute("length", String.valueOf(length));
            }
        }
        br.close();
    }

    /**
     * Process a mrk.tsv file
     */
    void processMrkFile(Reader reader) throws IOException {
        BufferedReader br = new BufferedReader(reader);
        Item dataSet = getDataSet();
        Item organism = getOrganism();
        Item geneticMap = null;
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
            if (fields[0].equals("MapName")) {
                String mapname = fields[1].trim();
                geneticMap = geneticMapMap.get(mapname);
                if (geneticMap==null) {
                    geneticMap = createItem("GeneticMap");
                    geneticMap.setAttribute("primaryIdentifier", mapname);
                    geneticMap.setReference("organism", organism);
                    geneticMapMap.put(mapname, geneticMap);
                }
            } else {
                // data line
                // 0                1                      2
                // MarkerIdentifier LinkageGroupIdentifier Position
                if (fields.length<3) {
                    throw new RuntimeException("File "+getCurrentFile().getName()+" data line does not contain three required fields marker, linkagegroup, position:"+line);
                }
                String markerId = fields[0].trim();
                String lgId = fields[1].trim();
                Double position = Double.parseDouble(fields[2]);
                Item linkageGroup = linkageGroupMap.get(lgId);
                if (linkageGroup==null) {
                    linkageGroup = createItem("LinkageGroup");
                    linkageGroup.setAttribute("identifier", lgId);
                    linkageGroup.setReference("organism", organism);
                    linkageGroup.setReference("geneticMap", geneticMap);
                    linkageGroupMap.put(lgId, linkageGroup);
                }
                Item marker = markerMap.get(markerId);
                if (marker==null) {
                    marker = createItem("GeneticMarker");
                    marker.setReference("organism", organism);
                }
                marker.setAttribute("secondaryIdentifier", markerId);
                marker.addToCollection("dataSets", dataSet);
                Item lgPosition = createItem("LinkageGroupPosition");
                lgPosition.setReference("marker", marker);
                lgPosition.setReference("linkageGroup", linkageGroup);
                lgPosition.setAttribute("position", String.valueOf(position));
                marker.addToCollection("linkageGroupPositions", lgPosition);
                markerMap.put(markerId, marker);
                linkageGroupPositions.add(lgPosition);
            }
        }
        br.close();
    }

    /**
     * Process a qtl.tsv file
     */
    void processQtlFile(Reader reader) throws IOException {
        BufferedReader br = new BufferedReader(reader);
        Item dataSet = getDataSet();
        Item organism = getOrganism();
        Item geneticMap = null;
        String intervalDescription = null;
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            String[] fields = line.split("\t");
            if (fields[0].equals("MapName")) {
                String mapname = fields[1].trim();
                geneticMap = geneticMapMap.get(mapname);
                if (geneticMap==null) {
                    geneticMap = createItem("GeneticMap");
                    geneticMap.setAttribute("primaryIdentifier", mapname);
                    geneticMap.setReference("organism", organism);
                    geneticMapMap.put(mapname, geneticMap);
                }
            } else if (fields[0].equals("IntervalDescription")) {
                intervalDescription = fields[1].trim();
            } else {
                // data line
                // 0          1     2            3     4   5   6
                // Identifier Trait LinkageGroup Start End LOD LikelihoodRatio
                if (fields.length<5) {
                    throw new RuntimeException("File "+getCurrentFile().getName()+" data line does not contain five required fields: identifier, trait, linkagegroup, start, end:"+line);
                }
                String identifier = fields[0].trim();
                String trait = fields[1].trim();
                String lgId = fields[2].trim();
                double start = Double.parseDouble(fields[3]);
                double end = Double.parseDouble(fields[4]);
                double lod = 0.0;
                if (fields.length>5 && fields[5]!=null && fields[5].trim().length()>0) lod = Double.parseDouble(fields[5]);
                double likelihoodRatio = 0.0;
                if (fields.length>6 && fields[6]!=null && fields[6].trim().length()>0) likelihoodRatio = Double.parseDouble(fields[6]);
                // QTL etc.
                Item qtl = createItem("QTL");
                qtl.setReference("organism", organism);
                qtl.addToCollection("dataSets", dataSet);
                qtl.setReference("geneticMap", geneticMap);
                if (intervalDescription!=null) qtl.setAttribute("intervalDescription", intervalDescription);
                qtl.setAttribute("identifier", identifier);
                Item phenotype = phenotypeMap.get(trait);
                if (phenotype==null) {
                    phenotype = createItem("Phenotype");
                    phenotype.setAttribute("primaryIdentifier", trait);
                    phenotype.setReference("organism", organism);
                    phenotype.addToCollection("dataSets", dataSet);
                    phenotypeMap.put(trait, phenotype);
                }
                qtl.setReference("phenotype", phenotype);
                Item linkageGroup = linkageGroupMap.get(lgId);
                if (linkageGroup==null) {
                    linkageGroup = createItem("LinkageGroup");
                    linkageGroup.setAttribute("identifier", lgId);
                    linkageGroup.setReference("organism", organism);
                    linkageGroupMap.put(lgId, linkageGroup);
                }
                qtl.setReference("linkageGroup", linkageGroup);
                qtl.setAttribute("start", String.valueOf(start));
                qtl.setAttribute("end", String.valueOf(end));
                if (lod>0.0) qtl.setAttribute("lod", String.valueOf(lod));
                if (likelihoodRatio>0.0) qtl.setAttribute("likelihoodRatio", String.valueOf(likelihoodRatio));
                qtlMap.put(identifier, qtl);
            }
        }
        br.close();
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
        store(strains.values());
        // local
        store(linkageGroupPositions);
        store(publicationMap.values());
        store(geneticMapMap.values());
        store(qtlMap.values());
        store(linkageGroupMap.values());
        store(markerMap.values());
        store(phenotypeMap.values());
    }

    /**
     * Extract the file type string from the filename.
     */
    String getFileType(String filename) {
        String[] fields = filename.split("\\.");
        if (fields[5].equals("expt")) {
            return "expt";
        } else if (fields[5].equals("mrk")) {
            return "mrk";
        } else if (fields[5].equals("qtl")) {
            return "qtl";
        } else {
            return null;
        }
    }
}
