package org.intermine.bio.dataconversion;

import java.io.File;
import java.io.IOException;
import java.io.Reader;

import java.util.Map;
import java.util.HashMap;

import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

/**
 * Class providing standard methods for Datastore file converters. Extend for your specific converter.
 *
 * NOTE: Items are NOT stored in these methods, just created. Your extending class should include the
 * following in the close() method:
 *
 * store(dataSource);
 * store(dataSets.values());
 * store(organisms.values());
 * store(strains.values());
 *
 * @author Sam Hokin
 */
public class DatastoreFileConverter extends FileConverter {

    // defaults for LIS datasource
    public static final String DEFAULT_DATASOURCE_NAME = "LIS Datastore";
    public static final String DEFAULT_DATASOURCE_URL = "https://legumeinfo.org/data/v2/";
    public static final String DEFAULT_DATASOURCE_DESCRIPTION =
        "A collaborative, community resource to facilitate crop improvement by integrating genetic, genomic, and trait data across legume species.";
    public static final String DEFAULT_DATASET_LICENCE = "ODC Public Domain Dedication and Licence (PDDL)";
	
    // Items to be stored by extending classes
    Item dataSource;
    Map<String,Item> organisms = new HashMap<>();  // keyed by TaxonID
    Map<String,Item> strains = new HashMap<>();    // keyed by identifier
    Map<String,Item> dataSets = new HashMap<>();   // keyed by name = filename

    // there is only one dataSource, set in project.xml
    String dataSourceName;         // optional
    String dataSourceUrl;          // required if Name given
    String dataSourceDescription;  // required if Name given

    // Set some DataSet attributes in project.xml
    String dataSetName;            // optional
    String dataSetLicence;         // optional
    String dataSetUrl;             // required
    String dataSetDescription;     // required

    DatastoreUtils dsu;

    /**
     * Create a new DatastoreFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public DatastoreFileConverter(ItemWriter writer, Model model) {
	super(writer, model);
        dsu = new DatastoreUtils();
        dataSource = getDataSource();
    }

    // Set DataSource fields in project.xml
    public void setDataSourceName(String name) {
        this.dataSourceName = name;
    }
    public void setDataSourceUrl(String url) {
        this.dataSourceUrl = url;
    }
    public void setDataSourceDescription(String description) {
        this.dataSourceDescription = description;
    }

    // Set DataSet fields in project.xml
    public void setDataSetName(String name) {
        this.dataSetName = name;
    }
    public void setDataSetUrl(String url) {
        this.dataSetUrl = url;
    }
    public void setDataSetDescription(String description) {
        this.dataSetDescription = description;
    }
    public void setDataSetLicence(String licence) {
	this.dataSetLicence = licence;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws IOException {
        // do nothing here
    }

    /**
     * Print out info about the current file or directory being processed.
     */
    static void printInfoBlurb(String blurb) {
	String outBlurb = "Processing "+blurb;
	String hashes = "";
	for (int i=0; i<outBlurb.length(); i++) hashes += "#";
        System.out.println(hashes);
        System.out.println(outBlurb);
        System.out.println(hashes);
    }

    /**
     * Get/add the DataSource Item.
     */
    public Item getDataSource() {
        if (dataSource==null) {
            // set defaults for LIS if not given
            if (dataSourceName==null) {
		dataSourceName = DEFAULT_DATASOURCE_NAME;
                dataSourceUrl = DEFAULT_DATASOURCE_URL;
                dataSourceDescription = DEFAULT_DATASOURCE_DESCRIPTION;
            }
            // create the DataSource once
            dataSource = createItem("DataSource");
            dataSource.setAttribute("name", dataSourceName);
            if (dataSourceUrl!=null) dataSource.setAttribute("url", dataSourceUrl);
            if (dataSourceDescription!=null) dataSource.setAttribute("description", dataSourceDescription);
        }
	return dataSource;
    }

    /**
     * Get/add a DataSet Item from the current file's containing directory, if not set.
     */
    public Item getDataSet() {
        if (dataSetName!=null) {
            return getDataSet(dataSetName);
        } else {
            return getDataSet(getCurrentFile().getParentFile().getName());
        }
    }

    /**
     * Get/add a DataSet Item from the given filename.
     */
    public Item getDataSet(String name) {
	if (dataSets.containsKey(name)) return dataSets.get(name);
        if (dataSource==null) {
	    throw new RuntimeException("DataSource is not initialized.");
	}
	if (dataSetUrl==null || dataSetDescription==null) {
	    throw new RuntimeException("You must set dataSetUrl and dataSetDescription in project.xml.");
	}
	// extract attributes
	String assemblyVersion = DatastoreUtils.extractAssemblyVersion(name);
	String annotationVersion = DatastoreUtils.extractAnnotationVersion(name);
	String dataSetVersion = null;
	if (assemblyVersion!=null && annotationVersion!=null) {
	    dataSetVersion = assemblyVersion+"."+annotationVersion;
	} else if (assemblyVersion!=null) {
	    dataSetVersion = assemblyVersion;
	}
	if (dataSetLicence==null) dataSetLicence = DEFAULT_DATASET_LICENCE;
	// create and return a new DataSet
	Item dataSet = createItem("DataSet");
	dataSet.setReference("dataSource", dataSource);
	dataSet.setAttribute("name", name);
        dataSet.setAttribute("url", dataSetUrl);
	dataSet.setAttribute("description", dataSetDescription);
	dataSet.setAttribute("licence", dataSetLicence);
	if (dataSetVersion!=null) dataSet.setAttribute("version", dataSetVersion);
	dataSets.put(name, dataSet);
	return dataSet;
    }

    /**
     * Get/add the organism Item by extracting gensp value from the current filename.
     */
    Item getOrganism() {
        return getOrganism(getGensp());
    }

    /**
     * Get the gensp string for the current filename.
     */
    String getGensp() {
        return DatastoreUtils.extractGensp(getCurrentFile().getName());
    }

    /**
     * Get/add the organism Item by referencing its taxonId value, as an integer to 
     * distinguish this method from the gensp method. (TaxonIDs are always integers.)
     */
    Item getOrganism(int intTaxonId) {
	String taxonId = String.valueOf(intTaxonId);
	if (organisms.containsKey(taxonId)) {
	    return organisms.get(taxonId);
	} else {
	    Item organism = createItem("Organism");
	    organism.setAttribute("taxonId", taxonId);
	    organisms.put(taxonId, organism);
	    return organism;
	}
    }
    
    /**
     * Get/add the organism Item associated with the given gensp value (e.g. "phavu").
     * Returns null if the gensp isn't resolvable to a taxonId.
     */
    Item getOrganism(String gensp) {
	String taxonId = dsu.getTaxonId(gensp);
	if (taxonId==null) {
	    return null;
	} else {
	    Item organism = getOrganism(Integer.parseInt(taxonId));
	    String genus = dsu.getGenus(gensp);
	    String species = dsu.getSpecies(gensp);
            organism.setAttribute("abbreviation", gensp);
            organism.setAttribute("genus", genus);
            organism.setAttribute("species", species);
            organism.setAttribute("name", genus+" "+species);
	    organism.setAttribute("shortName", genus.substring(0,1)+". "+species);
	    return organism;
        }
    }

    /**
     * Get the organism Item for a genus and species.
     */
    Item getOrganism(String genus, String species) {
        String gensp = genus.toLowerCase().substring(0,3) + species.toLowerCase().substring(0,2);
        return getOrganism(gensp);
    }

    /**
     * Get/add the strainItem associated with the current filename and organism.
     */
    Item getStrain(Item organism) {
	String identifier = DatastoreUtils.extractStrainIdentifier(getCurrentFile().getName());
	return getStrain(identifier, organism);
    }
    
    /**
     * Get/add the strain Item associated with the given strain name.
     * Sets the organism reference if created.
     */
    Item getStrain(String identifier, Item organism) {
        if (strains.containsKey(identifier)) {
            return strains.get(identifier);
        } else {
            Item strain = createItem("Strain");
            strain.setAttribute("identifier", identifier);
            strain.setReference("organism", organism);
            strains.put(identifier, strain);
	    return strain;
        }
    }
}
