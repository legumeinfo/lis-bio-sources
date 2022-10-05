package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.io.IOException;
import java.io.InputStream;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.Properties;

import org.apache.log4j.Logger;

import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Loads data from an LIS InterPro Scan GFF file.
 *
 * @author Sam Hokin
 */
public class IPRGFFConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(IPRGFFConverter.class);

    // local things to store
    Map<String,Item> proteins = new HashMap<>();
    List<Item> features = new ArrayList<>();

    /**
     * Create a new IPRGFFConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public IPRGFFConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
            setStrain();
        } else if (getCurrentFile().getName().endsWith("iprscan.gff3")) {
            System.out.println("## Processing "+getCurrentFile().getName());
            processGFF(reader);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        // references and collections
        for (Item protein : proteins.values()) {
            protein.setReference("organism", organism);
            protein.setReference("strain", strain);
        }
        for (Item feature : features) {
            feature.setReference("organism", organism);
            feature.setReference("strain", strain);
        }
        // collection items
        storeCollectionItems();
        // local stuff
        store(proteins.values());
        store(features);
    }

    /**
     * Process an IPR GFF file which has Proteins in the sequence column.
     * 0     1            2    3    4    5       6
     * medsa.XinJiangDaYe.gnm1.ann1.RKB9.iprscan.gff3
     */
    void processGFF(Reader reader) throws IOException {
        // spin through the file
        String line = null;
        BufferedReader br = new BufferedReader(reader);
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue;
            GFF3Record rec = new GFF3Record(line);
            // medsa.XinJiangDaYe.gnm1.ann1.MS_gene000000.t1 PANTHER protein_match 1 463 . + . Name=PTHR10178:SF14;status=T;ID=match$1047424_1_463;date=04-02-2021
            // medsa.XinJiangDaYe.gnm1.ann1.MS_gene000008.t1 ProSiteProfiles protein_match  10 264   . + . Name=PS50294;status=T;ID=match$461640_10_264;date=03-02-2021;
            //    signature_desc=Trp-Asp (WD) repeats circular profile.
            // medsa.XinJiangDaYe.gnm1.ann1.MS_gene000000.t1 Pfam    protein_hmm_match 3 88 4.8E-11 + . Name=PF03732;status=T;ID=match$1047423_3_88;date=04-02-2021;
            //    signature_desc=Retrotransposon gag protein;Target=PF03732 6 96;
            Item protein = getProtein(rec.getSequenceID());
            Item feature;
            if (rec.getType().equals("protein_match")) {
                feature = createItem("ProteinMatch");
                protein.addToCollection("proteinMatches", feature);
            } else if (rec.getType().equals("protein_hmm_match")) {
                feature = createItem("ProteinHmmMatch");
                protein.addToCollection("proteinHmmMatches", feature);
            } else {
                throw new RuntimeException("GFF record type "+rec.getType()+" is not supported by this loader.");
            }
            features.add(feature);
            // ID
            feature.setAttribute("primaryIdentifier", rec.getId());
            // location
            Item location = createItem("Location");
            location.setAttribute("start", String.valueOf(rec.getStart()));
            location.setAttribute("end", String.valueOf(rec.getEnd()));
            location.setReference("feature", feature);
            location.setReference("locatedOn", protein);
            // source
            feature.setAttribute("source", rec.getSource());
            // accession=Name
            feature.setAttribute("accession", rec.getNames().get(0));
            // status
            if (rec.getAttributes().containsKey("status")) {
                feature.setAttribute("status", rec.getAttributes().get("status").get(0));
            }
            // date
            if (rec.getAttributes().containsKey("date")) {
                feature.setAttribute("date", rec.getAttributes().get("date").get(0));
            }
            // target
            if (rec.getAttributes().containsKey("Target")) {
                feature.setAttribute("target", rec.getAttributes().get("Target").get(0));
            }
            // signatureDesc
            if (rec.getAttributes().containsKey("signature_desc")) {
                feature.setAttribute("signatureDesc", rec.getAttributes().get("signature_desc").get(0));
            }
            // store this line
            try {
                store(location);
            } catch (ObjectStoreException ex) {
                throw new RuntimeException(ex);
            }
        }
        br.close();
    }

    /**
     * Get/add a Protein Item, keyed by primaryIdentifier
     */
    public Item getProtein(String primaryIdentifier) {
        if (proteins.containsKey(primaryIdentifier)) {
            return proteins.get(primaryIdentifier);
        } else {
            Item protein = createItem("Protein");
            protein.setAttribute("primaryIdentifier", primaryIdentifier);
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, true);
	    if (secondaryIdentifier!=null) protein.setAttribute("secondaryIdentifier", secondaryIdentifier);
            proteins.put(primaryIdentifier, protein);
            return protein;
        }
    }
}
