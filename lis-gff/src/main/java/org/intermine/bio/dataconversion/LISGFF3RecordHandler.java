package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2019 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.util.List;
import java.util.Map;
import java.util.HashMap;

import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreWriter;
import org.intermine.objectstore.ObjectStoreWriterFactory;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.ItemFactory;
import org.intermine.xml.full.ItemHelper;

/**
 * Handle special cases when converting LIS data store GFF3 files.
 *
 * gene            ID=phavu.G19833.gnm1.ann1.Phvul.001G000100;
 * exon
 * mRNA            ID=phavu.G19833.gnm1.ann1.Phvul.001G000100.1;Parent=phavu.G19833.gnm1.ann1.Phvul.001G000100;
 * CDS             ID=phavu.G19833.gnm1.ann1.Phvul.001G000100.1.v1.0.CDS.1;Parent=phavu.G19833.gnm1.ann1.Phvul.001G000100.1;
 * five_prime_UTR  ID=phavu.G19833.gnm1.ann1.Phvul.001G000100.1.v1.0.five_prime_UTR.1;Parent=phavu.G19833.gnm1.ann1.Phvul.001G000100.1;
 * three_prime_UTR ID=phavu.G19833.gnm1.ann1.Phvul.001G000100.1.v1.0.three_prime_UTR.1;Parent=phavu.G19833.gnm1.ann1.Phvul.001G000100.1;
 * rRNA_primary_transcript ID=glyma.Zh13.gnm1.ann1.SoyZH13_CG009100.rRNA1;Name=SoyZH13_CG009100.rRNA1;Parent=glyma.Zh13.gnm1.ann1.SoyZH13_CG009100
 * tRNA_primary_transcript ID=glyma.Zh13.gnm1.ann1.SoyZH13_CG003000.tRNA1;Name=SoyZH13_CG003000.tRNA1;Parent=glyma.Zh13.gnm1.ann1.SoyZH13_CG003000
 * genetic_marker  ID=phavu.G19833.gnm1.mrk.3_12345;Name=3_12345;Alleles=A/T
 *
 * Parents:
 *
 * Exon.transcripts
 * CDS.transcript
 * FivePrimeUTR.transcripts
 * ThreePrimeUTR.transcripts
 * MRNA.gene
 * RRNAPrimaryTranscript.gene
 * TRNAPrimaryTranscript.gene
 *
 * @author Richard Smith
 * @author Sam Hokin
 */
public class LISGFF3RecordHandler extends GFF3RecordHandler {

    ItemFactory itemFactory;
    Map<String, String> aliases = new HashMap<>();
    Map<String, Integer> ids = new HashMap<>();
    int nextClsId = 0;

    // ProteinDomain Items
    Map<String,Item> proteinDomainMap = new HashMap<>();

    // OntologyTerm Items
    Map<String,Item> ontologyTermMap = new HashMap<>();

    /**
     * Create a new LISGFF3RecordHandler object.
     * @param tgtModel the target Model
     */
    public LISGFF3RecordHandler(Model tgtModel) {
        super(tgtModel);
        // refsAndCollections controls references and collections that are set from the Parent= attributes in the GFF3 file.
	// these should be consistent with sequence ontology!
        refsAndCollections.put("Exon", "transcripts");            // SO: exon->transcript_region->transcript
        refsAndCollections.put("CDS", "transcript");              // SO: CDS->mRNA_region->transcript
        refsAndCollections.put("FivePrimeUTR", "transcripts");    // SO: five_prime_UTR->UTR->mRNA_region
        refsAndCollections.put("ThreePrimeUTR", "transcripts");   // SO: three_prime_UTR->UTR->mRNA_region
        refsAndCollections.put("MRNA", "gene");                   // SO: mRNA->mature_transcript->transcript->gene_member_region
	refsAndCollections.put("RRNAPrimaryTranscript", "gene");  // SO: rRNA->ncRNA->mature_transcript->transcript->gene_member_region
	refsAndCollections.put("TRNAPrimaryTranscript", "gene");  // SO: tRNA->ncRNA->mature_transcript->transcript->gene_member_region
    }

    /**
     * {@inheritDoc}
     * type             id                                                         
     * ---------------- -----------------------------------------------------------
     * gene             phavu.G19833.gnm2.ann1.Phvul.003G111200                    
     * mRNA             phavu.G19833.gnm2.ann1.Phvul.003G111200.1
     * CDS              phavu.G19833.gnm2.ann1.Phvul.003G111100.1.CDS.1            
     * five_prime_UTR   phavu.G19833.gnm2.ann1.Phvul.003G111100.1.five_prime_UTR.3 
     * three_prime_UTR  phavu.G19833.gnm2.ann1.Phvul.003G111100.1.three_prime_UTR.1
     */
    @Override
    public void process(GFF3Record record) {
        String id = record.getId();
        String type = record.getType();
        Item feature = getFeature();
        String className = feature.getClassName();
        Map<String,List<String>> attributesMap = record.getAttributes();
        // only update feature if ID is present
        if (id!=null) {
            // 0     1      2    3    4     5          6 7   8
            // phavu.G19833.gnm2.ann1.Phvul.003G111100.1.CDS.1
            // 0     1            2    3      4                         5
            // medtr.jemalong_A17.gnm5.ann1_6.exon:MtrunA17Chr1g0187771.1
            // 0     1            2    3      4
            // medtr.jemalong_A17.gnm5.ann1_6.gene:MtrunA17CPg0492171
	    // 0     1            2    3
	    // vigun.IT97K-499-35.gnm1.1_0052;
            String[] parts = id.split("\\.");
            if (parts.length<4) {
                throw new RuntimeException("ID has too few dot-separated parts:"+id);
            }
            String gensp = parts[0];
            String strainId = parts[1];
            String assemblyVersion = parts[2];
	    String annotationVersion = null;
	    if (parts.length>3) annotationVersion = parts[3];

            // set other standard attributes
            feature.setAttribute("assemblyVersion", assemblyVersion);
            if (annotationVersion!=null) feature.setAttribute("annotationVersion", annotationVersion);

            // add marker type = SNP if it is a marker with length 1
            if (type.equals("genetic_marker") && (record.getStart()-record.getEnd())==0) {
                feature.setAttribute("type", "SNP");
            }

            // more specific attributes
            for (String key : attributesMap.keySet()) {
                List<String> attributes = attributesMap.get(key);
		if (key.equals("Name")) {
		    feature.setAttribute("secondaryIdentifier", attributes.get(0));
		} else if (key.equals("Note")) {
                    // Note=ATP binding protein... IPR002624 (...)%2C IPR027417 (...)%3B GO:0005524 (...)%2C GO:0006139 (...);
                    feature.setAttribute("description", attributes.get(0));
                } else if (type.equals("gene") && key.equals("Dbxref")) {
                    // Dbxref=Gene3D:G3DSA:3.40.50.300,InterPro:IPR002624,InterPro:IPR027417,PANTHER:PTHR10513,PANTHER:PTHR10513:SF6,Pfam:PF01712,Superfamily:SSF52540;
                    for (String term : attributes) {
                        String[] pieces = term.split(":");
                        if (pieces[0].equals("Gene3D")) {
                            // need Gene3D ProteinDomain Item
                        } else if (pieces[0].equals("InterPro")) {
                            String identifier = pieces[1];
			    Item proteinDomain = proteinDomainMap.get(identifier);
			    if (proteinDomain==null) {
				proteinDomain = converter.createItem("ProteinDomain");
				proteinDomain.setAttribute("primaryIdentifier", identifier);
				addItem(proteinDomain);
				proteinDomainMap.put(identifier, proteinDomain);
			    }
			    feature.addToCollection("proteinDomains", proteinDomain);
                        } else if (pieces[0].equals("PANTHER")) {
                            String identifier = pieces[1];
                            // need PANTHER ProteinDomain Item
                        } else if (pieces[0].equals("Pfam")) {
                            String identifier = pieces[1];
                            // need Pfam ProteinDomain Item
                        } else if (pieces[0].equals("Superfamily")) {
                            // need Superfamily ProteinDomain Item
                        }
                    }
                } else if (type.equals("gene") && key.equals("Ontology_term")) {
                    // Ontology_term=GO:0005524,GO:0006139,GO:0016773;
                    for (String term : attributes) {
                        if (term.startsWith("GO:")) {
			    String identifier = term;
			    Item goTerm = ontologyTermMap.get(identifier);
			    if (goTerm==null) {
				goTerm = converter.createItem("OntologyTerm");
				goTerm.setAttribute("identifier", identifier);
				addItem(goTerm);
				ontologyTermMap.put(identifier, goTerm);
			    }
			    Item goAnnotation = converter.createItem("OntologyAnnotation");
			    goAnnotation.setReference("subject", feature);
			    goAnnotation.setReference("ontologyTerm", goTerm);
			    addItem(goAnnotation);
                        }
                    }
		} else if (type.equals("genetic_marker") && key.equals("Alleles")) {
		    feature.setAttribute("alleles", attributes.get(0));
                } else if (key.equals("evid_id")) {
                    // [GAR_10012494]
                    // do nothing
                } else if (key.equals("geneFamily")) {
                    // geneFamily=phytozome_10_2.59141255
                    // gene family association
                }
            }
        }
    }
}
