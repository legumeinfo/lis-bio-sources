package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;

import java.util.Arrays;
import java.util.EmptyStackException;
import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.ncgr.jphyloio.Edge;
import org.ncgr.jphyloio.Node;

import info.bioinfweb.jphyloio.ReadWriteParameterMap;
import info.bioinfweb.jphyloio.ReadWriteParameterNames;
import info.bioinfweb.jphyloio.events.EdgeEvent;
import info.bioinfweb.jphyloio.events.JPhyloIOEvent;
import info.bioinfweb.jphyloio.events.LinkedLabeledIDEvent;
import info.bioinfweb.jphyloio.events.NodeEvent;
import info.bioinfweb.jphyloio.events.type.EventTopologyType;
import info.bioinfweb.jphyloio.factory.JPhyloIOReaderWriterFactory;
import info.bioinfweb.jphyloio.formats.newick.NewickEventReader;

import org.apache.log4j.Logger;

/**
 * Load phylogenetic tree data from LIS datastore Newick-format files.
 *
 * @author Sam Hokin
 */
public class PhylotreeFileConverter extends DatastoreFileConverter {

    private static final Logger LOG = Logger.getLogger(PhylotreeFileConverter.class);

    // local Items to store
    Map<String,Item> proteins = new HashMap<>();       // keyed by primaryIdentifier
    Map<String,Item> phylonodes = new HashMap<>();     // key=(phylotree.identifier).(Node.id)
    Map<String,Integer> childCounts = new HashMap<>(); // key=(phylotree.identifier).(Node.id)
    List<Item> geneFamilies = new LinkedList<>();      // unique per file
    List<Item> phylotrees = new LinkedList<>();        // unique per file
    List<Item> newicks = new LinkedList<>();           // unique per file

    /**
     * Create a new PhylotreeFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public PhylotreeFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) {
        // there is no README (yet) so hand-roll DataSet
        if (dataSet==null) {
            if (dataSetName==null || dataSetUrl==null || dataSetDescription==null) {
                throw new RuntimeException("ERROR: dataSetName, dataSetUrl, and dataSetDescription must be set in project.xml.");
            }
            dataSet = createItem("DataSet");
            dataSet.setAttribute("name", dataSetName);
            dataSet.setAttribute("url", dataSetUrl);
            dataSet.setAttribute("description", dataSetDescription);
            dataSet.setAttribute("synopsis", dataSetDescription);
            if (dataSetLicence!=null) {
                dataSet.setAttribute("licence", dataSetLicence);
            } else {
                dataSet.setAttribute("licence", DatastoreFileConverter.DEFAULT_DATASET_LICENCE);
            }
            dataSet.setReference("dataSource", dataSource);
        }
        try {
            processTreeFile();
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws Exception {
        store(dataSource);
        store(dataSet);
	store(geneFamilies);
        store(phylotrees);
        store(phylonodes.values());
	store(proteins.values());
        store(newicks);
    }

    /**
     * Process a file in the phylotree subdirectory, in Newick format.
     * legume.genefam.fam1.M65K.trees_ML_rooted/legfed_v1_0.L_6W30ML
     */
    void processTreeFile() throws IOException {
        JPhyloIOReaderWriterFactory factory = new JPhyloIOReaderWriterFactory();
        try {
            String formatID = factory.guessFormat(getCurrentFile());
            if (formatID==null) {
                throw new RuntimeException("File "+getCurrentFile().getName()+" format could not be determined.");
            } else if (!formatID.equals("info.bioinfweb.jphyloio.newick")) {
                throw new RuntimeException("File "+getCurrentFile().getName()+" format is not a Newick file: "+formatID);
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
        // create this Phylotree and GeneFamily
        // 0           1
        // legfed_v1_0.L_6W30ML
        String identifier = getCurrentFile().getName();
        Item phylotree = createItem("Phylotree");
        phylotrees.add(phylotree);
        phylotree.setAttribute("primaryIdentifier", identifier);
        Item geneFamily = createItem("GeneFamily");
        geneFamilies.add(geneFamily);
        geneFamily.setAttribute("primaryIdentifier", identifier);
        geneFamily.setReference("phylotree", phylotree);
        phylotree.setReference("geneFamily", geneFamily);
        // keep track of the leaf nodes
        int numLeaves = 0;
        // TODO:
        // phylotree.setAttribute("numLeaves", String.valueOf(tree.getLeafCount()));
        // phylotree.setAttribute("height", String.valueOf(tree.getHeight()));
        // set jphyloio parameters for reader
        ReadWriteParameterMap parameters = new ReadWriteParameterMap();
        // Use OTU labels as node labels if no node label is present.
        parameters.put(ReadWriteParameterNames.KEY_USE_OTU_LABEL, true);  
        // This parameter defines if cross links between nodes (defined by the clade_relation tag of PhyloXML) should be
        // modeled as metadata attached to a node or if the whole phylogeny shall be interpreted as a phylogenetic network.
        // Since the network interpretation is the default, we need to set this parameter in order to receive tree events
        // and not network events.
        parameters.put(ReadWriteParameterNames.KEY_PHYLOXML_CONSIDER_PHYLOGENY_AS_TREE, true);
        // create the reader
        NewickEventReader reader = new NewickEventReader(getCurrentFile(), parameters);
        // This loop will run until all events of the JPhyloIO reader are consumed (and the end of the document is reached). 
        while (reader.hasNextEvent()) {
            JPhyloIOEvent event = reader.next();
            // we only care about STARTs
            if (event.getType().getTopologyType().equals(EventTopologyType.START)) {
                switch (event.getType().getContentType()) {
                case EDGE:
                    // Indicates an edge in a phylogenetic tree or network.
                    Edge e = new Edge(event.asEdgeEvent());
                    processEdge(e, identifier);
                    break;
                case NODE:
                    // Indicates a node in a phylogenetic tree or network.
                    Node n = new Node(event.asNodeEvent());
                    Item phylonode = getPhylonode(n, identifier);
                    if (phylonode==null) return;
                    phylonode.setReference("tree", phylotree);
                    phylotree.addToCollection("nodes", phylonode);
                    if (n.isFeature()) numLeaves++;
                    break;
                case TREE_NETWORK_GROUP:
                    // Indicates the start or the end of a sequence of phylogenetic trees and network.
                    break;
                case TREE_NETWORK_SET:
                    // Indicates the start or end of a sequence of SET_ELEMENT events that define a set of trees and networks.
                    break;
                case UNKNOWN_COMMAND:
                    // Events of this type are used by some readers to provide the application with contents of unknown commands in a format.
                    break;
	        default:
                    // everything else
                    break;
                }
            }
        }
        reader.close();
        // set numLeaves for this tree
        phylotree.setAttribute("numLeaves", String.valueOf(numLeaves));
        // now store the Newick file contents
        Item newick = createItem("Newick");
        newicks.add(newick);
        newick.setAttribute("identifier", identifier);
        newick.setReference("phylotree", phylotree);
        newick.setReference("geneFamily", geneFamily);
        String contents = "";
        String line = null;
        BufferedReader br = new BufferedReader(new FileReader(getCurrentFile()));
        while ((line=br.readLine())!=null) {
            contents += line;
        }
        br.close();
        newick.setAttribute("contents", contents);
    }

    /**
     * Return true if the given identifier is in LIS full yuck format.
     */
    boolean isFullYuck(String primaryIdentifier) {
        return primaryIdentifier.split("\\.").length>=5;
    }
        
    /**
     * Process the source and target nodes of a given Edge.
     */
    void processEdge(Edge e, String identifier) {
        String sourceKey = identifier+"."+e.sourceId;
        String targetKey = identifier+"."+e.targetId;
        Item sourcenode;
        if (phylonodes.containsKey(sourceKey)) {
            sourcenode = phylonodes.get(sourceKey);
        } else {
            sourcenode = createItem("Phylonode");
            sourcenode.setAttribute("identifier", sourceKey);
            phylonodes.put(sourceKey, sourcenode);
        }
        Item targetnode;
        if (phylonodes.containsKey(targetKey)) {
            targetnode = phylonodes.get(targetKey);
        } else {
            targetnode = createItem("Phylonode");
            targetnode.setAttribute("identifier", targetKey);
            phylonodes.put(targetKey, targetnode);
        }
        targetnode.setAttribute("isRoot", "false"); // we have a parent
        if (e.hasLength()) {
            targetnode.setAttribute("length", String.valueOf(e.length));
        }
        targetnode.setReference("parent", sourcenode);
        sourcenode.addToCollection("children", targetnode);
        int count = 1;
        if (childCounts.containsKey(sourceKey)) {
            count = childCounts.get(sourceKey) + 1;
        }
        childCounts.put(sourceKey, count);
        sourcenode.setAttribute("numChildren", String.valueOf(count));
    }
    
    /**
     * Get/add a Protein Item, keyed by primaryIdentifier
     */
    Item getProtein(String primaryIdentifier) {
        if (proteins.containsKey(primaryIdentifier)) {
            return proteins.get(primaryIdentifier);
        } else {
            Item protein = createItem("Protein");
            proteins.put(primaryIdentifier, protein);
            protein.setAttribute("primaryIdentifier", primaryIdentifier);
            return protein;
        }
    }

    /**
     * Get/add/update a Phylonode Item, keyed by identifier.id
     *
     * @param node the Node
     * @param identifier the identifier of the PhyloTree
     * @return the Phylonode Item
     */
    Item getPhylonode(Node node, String identifier) {
        if (node==null) {
            throw new RuntimeException("null Node passed into getPhylonode.");
        }
        String key = identifier+"."+node.id;
        Item phylonode;
        if (phylonodes.containsKey(key)) {
            phylonode = phylonodes.get(key);
        } else {
            phylonode = createItem("Phylonode");
            phylonode.setAttribute("identifier", key);
            phylonode.setAttribute("isRoot", "true");  // default
            phylonode.setAttribute("isLeaf", "false"); // default
            phylonodes.put(key, phylonode);
        }
        if (node.isFeature()) {
            // convert blanks to underscores
            String label = node.label.replace(" ","_");
            // DEBUG
            System.out.println(node.label+"\t"+label);
            //
            phylonode.setAttribute("isLeaf", "true");
            phylonode.setAttribute("identifier", label);
            // some nodes hold proteins
            if (isFullYuck(label)) {
                Item protein = getProtein(label);
                phylonode.setReference("protein", protein);
            }
        }
        return phylonode;
    }

}
