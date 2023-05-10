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

import org.ncgr.datastore.DescriptionGenus;

/**
 * Loads data from an LIS datastore description_Genus.yml file.
 *
 * This creates a single genus Organism.
 * Species are ignored since their taxon ID is not provided.
 * Resources are ignored as well.
 *
 * @author Sam Hokin
 */
public class GenusDescriptionFileConverter extends FileConverter {

    Item dataSet;
    Item dataSource;
    Item organism;

    String dataSetName, dataSetUrl, dataSetDescription;

    /**
     * dataSetUrl is set in project.xml
     */
    public void setDataSetName(String name) {
        this.dataSetName = name;
    }
    /**
     * dataSetUrl is set in project.xml
     */
    public void setDataSetUrl(String url) {
        this.dataSetUrl = url;
    }
    /**
     * dataSetDescription is set in project.xml
     */
    public void setDataSetDescription(String description) {
        this.dataSetDescription = description;
    }

    /**
     * Create a new GenusDescriptionFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GenusDescriptionFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (!getCurrentFile().getName().startsWith("description") && !getCurrentFile().getName().endsWith("yml")) return;
        DescriptionGenus genus = DescriptionGenus.parse(reader);
        // genus Organism
        organism = createItem("Organism");
        organism.setAttribute("taxonId", String.valueOf(genus.taxid));
        organism.setAttribute("genus", genus.genus);
        organism.setAttribute("description", genus.description);
        if (genus.commonName!=null) organism.setAttribute("commonName", genus.commonName);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        // stock DataSource
        Item dataSource = createItem("DataSource");
        dataSource.setAttribute("name", DatastoreFileConverter.DEFAULT_DATASOURCE_NAME);
        dataSource.setAttribute("url", DatastoreFileConverter.DEFAULT_DATASOURCE_URL);
        dataSource.setAttribute("description", DatastoreFileConverter.DEFAULT_DATASOURCE_DESCRIPTION);
        // DataSet from project.xml params
        Item dataSet = createItem("DataSet");
        dataSet = createItem("DataSet");
        dataSet.setReference("dataSource", dataSource);
        dataSet.setAttribute("name", dataSetName);
        dataSet.setAttribute("description", dataSetDescription);
        dataSet.setAttribute("synopsis", dataSetDescription);
        dataSet.setAttribute("url", dataSetUrl);
        // organism and strains
        organism.addToCollection("dataSets", dataSet);
        // store
        store(dataSource);
        store(dataSet);
	store(organism);
    }
}
