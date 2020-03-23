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
        if (dsu.isSupercontig(converter.getTaxonId(), converter.getStrainIdentifier(), identifier)) {
            Item seq = converter.createItem("Supercontig");
            seq.setAttribute("primaryIdentifier", identifier);
            return seq;
        } else {
            Item seq = converter.createItem("Chromosome");
            seq.setAttribute("primaryIdentifier", identifier);
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
        if (dsu.isSupercontig(converter.getTaxonId(), converter.getStrainIdentifier(), identifier)) {
            Item seq = converter.createItem("Supercontig");
            seq.setAttribute("primaryIdentifier", identifier);
            return seq;
        } else {
            Item seq = converter.createItem("Chromosome");
            seq.setAttribute("primaryIdentifier", identifier);
            return seq;
        }
    }

}
