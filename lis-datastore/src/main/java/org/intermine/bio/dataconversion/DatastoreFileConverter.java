package org.intermine.bio.dataconversion;

import java.io.File;
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
 * NOTE: Items are NOT stored in these methods, just created.
 *
 * @author Sam Hokin
 */
public class DatastoreFileConverter extends FileConverter {

    // defaults for LIS datasource
    public static final String DEFAULT_DATASOURCE_NAME = "LIS Datastore";
    public static final String DEFAULT_DATASOURCE_URL = "https://legumeinfo.org/data/public/";
    public static final String DEFAULT_DATASOURCE_DESCRIPTION =
        "A collaborative, community resource to facilitate crop improvement by integrating genetic, genomic, and trait data across legume species.";
    public static final String DEFAULT_DATASET_LICENCE = "ODC Public Domain Dedication and Licence (PDDL)";
	
    // store Items in maps if they may be read more than once
    Map<String,Item> organisms = new HashMap<>();
    Map<String,Item> strains = new HashMap<>();
    Map<String,Item> dataSets = new HashMap<>();

    // there is only one dataSource, set in project.xml
    Item dataSource;
    String dataSourceName;         // optional
    String dataSourceUrl;          // required if Name given
    String dataSourceDescription;  // required if Name given

    // Set some DataSet attributes in project.xml
    String dataSetName;            // optional
    String dataSetLicence;         // optional
    String dataSetUrl;             // required
    String dataSetDescription;     // required

    /**
     * Create a new DatastoreFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public DatastoreFileConverter(ItemWriter writer, Model model) {
	super(writer, model);
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
    public void process(Reader reader) throws Exception {
	// do nothing here
    }

    /**
     * Print out info about the current file or directory being processed.
     */
    static void printInfoBlurb(String blurb) {
        System.out.println("####################################################################################################################");
        System.out.println("Processing file/dir:"+blurb);
        System.out.println("####################################################################################################################");
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
     * Get/add the DataSet Item.
     * dataSetUrl and dataSetDescription MUST be set in project.xml.
     */
    public Item getDataSet() {
	if (dataSource==null) {
	    throw new RuntimeException("DataSource is not initialized.");
	}
	if (dataSetUrl==null || dataSetDescription==null) {
	    throw new RuntimeException("You must set dataset.url and dataset.description in project.xml.");
	}
	// support optional data.set.name in project.xml
	String assemblyVersion = null;
	String annotationVersion = null;
	String name = dataSetName;
	if (name==null) {
	    name = getCurrentFile().getName();
	    assemblyVersion = extractAssemblyVersion(name);
	    annotationVersion = extractAnnotationVersion(name);
	}
	if (dataSetLicence==null) dataSetLicence = DEFAULT_DATASET_LICENCE;
	// return an existing DataSet
        if (dataSets.containsKey(name)) return dataSets.get(name);
	// create and return a new DataSet
	Item dataSet = createItem("DataSet");
	dataSet.setAttribute("name", name);
        dataSet.setAttribute("url", dataSetUrl);
	dataSet.setAttribute("description", dataSetDescription);
	dataSet.setAttribute("licence", dataSetLicence);
	if (assemblyVersion!=null && annotationVersion!=null) {
	    dataSet.setAttribute("version", assemblyVersion+"."+annotationVersion);
	} else if (assemblyVersion!=null) {
	    dataSet.setAttribute("version", assemblyVersion);
	}
	dataSet.setReference("dataSource", dataSource);
	dataSets.put(name, dataSet);
	return dataSet;
    }

    /**
     * Get/add the organism Item by extracting gensp value from the current filename.
     */
    Item getOrganism() {
	return getOrganism(extractGensp(getCurrentFile().getName()));
    }
    
    /**
     * Get/add the organism Item associated with the given gensp value (e.g. "phavu").
     * Returns null if the gensp isn't resolvable.
     */
    Item getOrganism(String gensp) {
        Item organism = null;
        if (organisms.containsKey(gensp)) {
            organism = organisms.get(gensp);
        } else {
	    DatastoreUtils dsu = new DatastoreUtils();
            String taxonId = dsu.getTaxonId(gensp);
            String genus = dsu.getGenus(gensp);
            String species = dsu.getSpecies(gensp);
            organism = createItem("Organism");
            organism.setAttribute("abbreviation", gensp);
            organism.setAttribute("taxonId", taxonId);
            organism.setAttribute("genus", genus);
            organism.setAttribute("species", species);
            organism.setAttribute("name", genus+" "+species);
            organisms.put(gensp, organism);
        }
        return organism;
    }

    /**
     * Get the organism Item for a string array like ["something","Genus","species"]
     */
    Item getOrganism(String[] threeparts) {
        String genus = threeparts[1];
        String species = threeparts[2];
        String gensp = genus.toLowerCase().substring(0,3) + species.toLowerCase().substring(0,2);
        return getOrganism(gensp);
    }

    /**
     * Get/add the strain Item associated with the given strain name.
     * Sets the organism reference if created.
     */
    Item getStrain(String strainId, Item organism) {
        Item strain;
        if (strains.containsKey(strainId)) {
            strain = strains.get(strainId);
        } else {
            strain = createItem("Strain");
            strain.setAttribute("identifier", strainId);
            strain.setReference("organism", organism);
            strains.put(strainId, strain);
        }
        return strain;
    }

    /**
     * Extract the gensp string from the given filename.
     * gensp.strain.assy.anno.key.content.ext
     */
    public static String extractGensp(String filename) {
	String[] fields = filename.split("\\.");
	if (fields.length>1) {
	    return fields[0];
	} else {
	    return null;
	}
    }

    /**
     * Extract the Strain identifier from the given filename.
     * gensp.strain.assy.anno.key.content.ext
     */
    public static String extractStrainIdentifier(String filename) {
	String[] fields = filename.split("\\.");
	if (fields.length>1) {
	    return fields[1];
	} else {
	    return null;
	}
    }

    /**
     * Extract the assembly version from the given filename.
     * gensp.strain.assy.anno.key.content.ext
     */
    public static String extractAssemblyVersion(String filename) {
        String[] fields = filename.split("\\.");
        if (fields.length>2 && fields[2].startsWith("gnm")) {
            return fields[2];
        } else {
            return null;
        }
    }

    /**
     * Extract the annotation version from the given filename.
     * gensp.strain.assy.anno.key.content.ext
     */
    public static String extractAnnotationVersion(String filename) {
        String[] fields = filename.split("\\.");
        if (fields.length>3 && fields[3].startsWith("ann")) {
            return fields[3];
        } else {
            return null;
        }
    }

    /**
     * Extract the secondaryIdentifier from a full-yuck LIS identifier.
     * Set isAnnotationFeature=true if this is an annotation feature (at least five dot-separated parts) as opposed to an assembly feature (at least four dot-separated parts).
     *
     * isAnnotationFeature==true:
     * 0      1      2   3 
     * genesp.strain.gnm.secondaryIdentifier
     *
     * isAnnotationFeature==false:
     * 0      1      2   3   4
     * genesp.strain.gnm.ann.secondaryIdentifier
     *
     * @param   lisIdentifier the LIS full-yuck identifier
     * @param   isAnnotationFeature true if this is an annotation feature with five or more dot-separated parts
     * @returns the secondaryIdentifier
     */
    public static String extractSecondaryIdentifier(String lisIdentifier, boolean isAnnotationFeature) {
	String[] fields = lisIdentifier.split("\\.");
	if (isAnnotationFeature && fields.length>=5) {
	    String secondaryIdentifier = fields[4];
	    for (int i=5; i<fields.length; i++) {
		secondaryIdentifier += "."+fields[i];
	    }
	    return secondaryIdentifier;
	} else if (!isAnnotationFeature && fields.length>=4) {
	    String secondaryIdentifier = fields[3];
	    for (int i=4; i<fields.length; i++) {
		secondaryIdentifier += "."+fields[i];
	    }
	    return secondaryIdentifier;
	} else {
	    return null;
	}
    }

    /**
     * Form a primaryIdentifier from the secondaryIdentifier, gensp, strain, assembly and annotation
     */
    public static String formPrimaryIdentifier(String gensp, String strain, String assembly, String annotation, String secondaryIdentifier) {
        return gensp+"."+strain+"."+assembly+"."+annotation+"."+secondaryIdentifier;
    }
}
