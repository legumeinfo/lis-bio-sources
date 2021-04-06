package org.intermine.bio.dataconversion;

import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.xml.full.Item;

/**
 * This implementation creates Chromosome items by default and Supercontig items if the identifier contains appropriate text.
 *
 * @author Sam Hokin
 */
public class LISGFF3SeqHandler extends GFF3SeqHandler {
    /**
     * For the given GFF3Converter and sequence identifier, make a new sequence object. 
     * This implementation creates Chromosome items by default and Supercontig items if identified as such by DatastoreUtils.isSupercontig().
     * @param converter the current GFF3Converter
     * @param identifier the identifier of the sequence from the GFF file
     * @return a new sequence Item
     * @override
     */
    public Item makeSequenceItem(GFF3Converter converter, String identifier) {
        DatastoreUtils dsu = new DatastoreUtils();
        String taxonId = converter.getTaxonId();
        String strainIdentifier = converter.getStrainIdentifier();
        String assemblyVersion = DatastoreUtils.extractAssemblyVersion(identifier);
        String gensp = dsu.getGensp(taxonId);
        if (dsu.isSupercontig(gensp, strainIdentifier, identifier)) {
            if (assemblyVersion==null) throw new RuntimeException("assemblyVersion=null for identifier="+identifier);
            Item seq = converter.createItem("Supercontig");
            seq.setAttribute("primaryIdentifier", identifier);
            seq.setAttribute("assemblyVersion", assemblyVersion);
            return seq;
        } else {
            if (assemblyVersion==null) throw new RuntimeException("assemblyVersion=null for identifier="+identifier);
            Item seq = converter.createItem("Chromosome");
            seq.setAttribute("primaryIdentifier", identifier);
            seq.setAttribute("assemblyVersion", assemblyVersion);
            return seq;
        }
    }

    /**
     * For the given GFF3Converter and sequence identifier, make a new sequence object.
     * This implementation creates Chromosome items by default and Supercontig items if the identifier contains "scaffold", etc.
     * @param converter the current GFF3Converter
     * @param record record from the GFF file
     * @param identifier the identifier of the sequence from the GFF file
     * @return a new sequence Item
     * @override
     */
    public Item makeSequenceItem(GFF3Converter converter, String identifier, GFF3Record record) {
        DatastoreUtils dsu = new DatastoreUtils();
        String taxonId = converter.getTaxonId();
        String strainIdentifier = converter.getStrainIdentifier();
        String assemblyVersion = DatastoreUtils.extractAssemblyVersion(identifier);
        String gensp = dsu.getGensp(taxonId);
        if (dsu.isSupercontig(gensp, strainIdentifier, identifier)) {
            if (assemblyVersion==null) throw new RuntimeException("assemblyVersion=null for identifier="+identifier);
            Item seq = converter.createItem("Supercontig");
            seq.setAttribute("primaryIdentifier", identifier);
            seq.setAttribute("assemblyVersion", assemblyVersion);
            return seq;
        } else {
            if (assemblyVersion==null) throw new RuntimeException("assemblyVersion=null for identifier="+identifier);
            Item seq = converter.createItem("Chromosome");
            seq.setAttribute("primaryIdentifier", identifier);
            seq.setAttribute("assemblyVersion", assemblyVersion);
            return seq;
        }
    }
}
