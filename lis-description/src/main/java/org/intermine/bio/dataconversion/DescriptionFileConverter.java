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
import java.io.File;
import java.io.Reader;
import java.io.IOException;

import java.util.List;
import java.util.ArrayList;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.ncgr.datastore.DescriptionSpecies;
import org.ncgr.datastore.DescriptionStrain;

/**
 * Loads data from an LIS datastore description_Genus_species.yml file.
 *
 * This creates a single Organism (species) with multiple associated Strains.
 * Resources are ignored.
 *
 * @author Sam Hokin
 */
public class DescriptionFileConverter extends FileConverter {

    Item organism;
    List<Item> strains = new ArrayList<>();
	
    private static final Logger LOG = Logger.getLogger(DescriptionFileConverter.class);

    /**
     * Create a new DescriptionFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public DescriptionFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (!getCurrentFile().getName().startsWith("description") && !getCurrentFile().getName().endsWith("yml")) return;
        DescriptionSpecies species = DescriptionSpecies.parse(reader);
        organism = createItem("Organism");
        organism.setAttribute("taxonId", String.valueOf(species.taxid));
        organism.setAttribute("genus", species.genus);
        organism.setAttribute("species", species.species);
        organism.setAttribute("abbreviation", species.abbrev);
        organism.setAttribute("description", species.description);
        if (species.commonName!=null) organism.setAttribute("commonName", species.commonName);
        // DEBUG
        System.out.println("taxonId="+species.taxid);
        System.out.println("genus="+species.genus);
        System.out.println("abbreviation="+species.abbrev);
        System.out.println("description="+species.description);
        System.out.println("commonName="+species.commonName);
        //
        // extras
        organism.setAttribute("name", species.genus+" "+species.species);
        organism.setAttribute("shortName", species.genus.charAt(0)+". "+species.species);
        for (DescriptionStrain strain : species.strains) {
            Item strainItem = createItem("Strain");
            strainItem.setAttribute("identifier", strain.identifier);
            strainItem.setAttribute("accession", strain.accession);
            strainItem.setAttribute("name", strain.name);
            strainItem.setAttribute("origin", strain.origin);
            strainItem.setAttribute("description", strain.description);
            strainItem.setReference("organism", organism);
            // DEBUG
            System.out.println("identifier="+strain.identifier);
            System.out.println("accession="+strain.accession);
            System.out.println("name="+strain.name);
            System.out.println("origin="+strain.origin);
            System.out.println("description="+strain.description);
            //
            strains.add(strainItem);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
	store(organism);
        store(strains);
    }
}
