package org.intermine.bio.dataconversion;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;
import java.io.UnsupportedEncodingException;

import java.net.MalformedURLException;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.Properties;

import javax.xml.parsers.ParserConfigurationException;

import org.intermine.bio.util.BioConverterUtil;
import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.JSONValue;
import org.json.simple.parser.ParseException;

import org.xml.sax.SAXException;

import org.ncgr.datastore.Readme;
import org.ncgr.crossref.WorksQuery;

/**
 * Class providing standard objects and methods for Datastore file converters. Extend for your specific converter.
 *
 * A DatastoreFileConverter operates a single Datastore collection, which is a DataSet.
 * There is one organism, as indicated in the README taxid.
 * There is one strain, indicated (if so) by the first piece of the collection identifier.
 *
 * NOTE 1: collection Items created here are stored by calling storeCollectionItems() in the extending class.
 *
 * NOTE 2: Extending classes need not set a dataSet reference or add to the dataSets collection.
 *         That is automatically done during store() by the BioStoreHook when the README is parsed.
 *
 * @author Sam Hokin
 */
public abstract class DatastoreFileConverter extends FileConverter {

    // defaults for LIS datasource
    public static final String DEFAULT_DATASOURCE_NAME = "LIS Datastore";
    public static final String DEFAULT_DATASOURCE_URL = "https://legumeinfo.org/data/v2/";
    public static final String DEFAULT_DATASOURCE_DESCRIPTION =
        "A collaborative, community resource to facilitate crop improvement by integrating genetic, genomic, and trait data across legume species.";
    public static final String DEFAULT_DATASET_LICENCE = "ODC Public Domain Dedication and Licence (PDDL)";
	
    static final String ORGANISM_PROP_FILE = "organism_config.properties";

    // the README file content
    Readme readme;
    
    // Items common to collections; these are stored with storeCollectionItems()
    Item dataSource;
    Item dataSet;
    Item organism;
    Item strain;
    Item publication;
    List<Item> authors = new ArrayList<>();

    // There is only one DataSource, set in project.xml, or default
    String dataSourceName;         // required
    String dataSourceUrl;          // required
    String dataSourceDescription;  // required

    // These can be set in project.xml
    String dataSetName;            // optional
    String dataSetDescription;     // optional
    String dataSetUrl;             // required in project.xml
    String dataSetLicence;         // optional

    // other attributes we may need for Items
    int taxonId;
    String gensp;
    String strainIdentifier;
    String assemblyVersion;
    String annotationVersion;

    Map<String,String> taxonIdGenus = new HashMap<>();
    Map<String,String> taxonIdSpecies = new HashMap<>();
    Map<String,String> taxonIdGensp = new HashMap<>();
    Map<String,String> genspTaxonId = new HashMap<>();
    Map<String,String> genusSpeciesTaxonId = new HashMap<>();

    List<String> chromosomePrefixes = new ArrayList<>();              // from README
    List<String> supercontigPrefixes = new ArrayList<>();             // from README

    /**
     * Create a new DatastoreFileConverter
     * Call with null,null arguments for use in classes that do not extend DatastoreFileConverter
     *
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public DatastoreFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
        setDataSource();
        loadOrganismProperties();
    }

    /**
     * dataSourceName can be set in project.xml but is usually the default
     */
    public void setDataSourceName(String name) {
        this.dataSourceName = name;
    }
    /**
     * dataSourceUrl can be set in project.xml but is usually the default
     */
    public void setDataSourceUrl(String url) {
        this.dataSourceUrl = url;
    }
    /**
     * dataSourceDescription can be set in project.xml but is usually the default
     */
    public void setDataSourceDescription(String description) {
        this.dataSourceDescription = description;
    }

    /**
     * dataSetName is usually README.identifier
     */
    public void setDataSetName(String name) {
        this.dataSetName = name;
    }

    /**
     * dataSetUrl must be set in project.xml
     */
    public void setDataSetUrl(String url) {
        this.dataSetUrl = url;
    }

    /**
     * dataSetDescription is usually README.description
     */
    public void setDataSetDescription(String description) {
        this.dataSetDescription = description;
    }
    
    /**
     * dataSetLicence can be set in project.xml but is usually the default
     */
    public void setDataSetLicence(String licence) {
	this.dataSetLicence = licence;
    }

    /**
     * Process the README file for the common collection Items stored in this class.
     * This also sets the BioStoreHook which automatically associates dataSet with stored Items.
     * Note: we get strain from the collection identifier in setStrain(), not from the README.
     */
    void processReadme(Reader reader) throws IOException {
        System.out.println("## Processing "+getCurrentFile().getName());
        readme = Readme.parse(reader);
        // check required stuff
        if (readme.identifier==null ||
            readme.taxid==0 ||
            readme.scientific_name_abbrev==null ||
            readme.synopsis==null ||
            readme.description==null
            ) {
            throw new RuntimeException("ERROR in README: a required field is missing. "+
                                       "Required fields are: identifier, taxid, scientific_name_abbrev, synopsis, description");
        }
        if (dataSetUrl==null) {
            throw new RuntimeException("ERROR: dataSetUrl must be set in project.xml.");
        }
        // local vars
        taxonId = readme.taxid;
        gensp = readme.scientific_name_abbrev;
        dataSetName = readme.identifier;
        dataSetDescription = readme.description;
        strainIdentifier = DatastoreUtils.extractStrainIdentifierFromCollection(readme.identifier);
        assemblyVersion = DatastoreUtils.extractAssemblyVersionFromCollection(readme.identifier);
        annotationVersion = DatastoreUtils.extractAnnotationVersionFromCollection(readme.identifier);
        if (readme.chromosome_prefix != null) {
            for (String prefix : readme.chromosome_prefix.split(",")) {
                chromosomePrefixes.add(prefix.trim());
            }
        }
        if (readme.supercontig_prefix != null) {
            for (String prefix : readme.supercontig_prefix.split(",")) {
                supercontigPrefixes.add(prefix.trim());
            }
        }
        // Organism - if not Fabacea, which isn't an organism
        if (taxonId != 3803) {
            organism = createItem("Organism");
            organism.setAttribute("taxonId", String.valueOf(taxonId));
        }
        // Publication - optional!
        if (readme.publication_doi == null) {
            System.err.println("### WARNING: readme.publication_doi IS NULL!");
        } else {
            try {
                createPublication();
            } catch (Exception ex) {
                throw new RuntimeException(ex);
            }
        }
        // DataSet
        dataSet = createItem("DataSet");
        dataSet.setReference("dataSource", dataSource);
        dataSet.setAttribute("name", dataSetName);
        dataSet.setAttribute("description", dataSetDescription);
        dataSet.setAttribute("synopsis", readme.synopsis);
        if (publication!=null) dataSet.setReference("publication", publication);
        if (assemblyVersion!=null && annotationVersion!=null) {
            dataSet.setAttribute("version", assemblyVersion+"."+annotationVersion);
        } else if (assemblyVersion!=null) {
            dataSet.setAttribute("version", assemblyVersion);
        }
        if (dataSetLicence!=null) {
            dataSet.setAttribute("licence", dataSetLicence);
        } else {
            dataSet.setAttribute("licence", DEFAULT_DATASET_LICENCE);
        }
        dataSet.setAttribute("url", dataSetUrl);
        // set BioStoreHook
        setStoreHook(new BioStoreHook(getModel(), dataSet.getIdentifier(), dataSource.getIdentifier(), BioConverterUtil.getOntology(this)));
    }

    /**
     * Set the instance DataSource from project.xml properties or defaults.
     */
    void setDataSource() {
        if (dataSource!=null) return;
        // set defaults for LIS if not given in project.xml
        if (dataSourceName==null) dataSourceName = DEFAULT_DATASOURCE_NAME;
        if (dataSourceUrl==null) dataSourceUrl = DEFAULT_DATASOURCE_URL;
        if (dataSourceDescription==null) dataSourceDescription = DEFAULT_DATASOURCE_DESCRIPTION;
        // create the DataSource once
        dataSource = createItem("DataSource");
        dataSource.setAttribute("name", dataSourceName);
        dataSource.setAttribute("url", dataSourceUrl);
        dataSource.setAttribute("description", dataSourceDescription);
    }

    /**
     * Set the instance Strain Item from readme.identifier.
     * NOTE: processReadme() must already have been run.
     * Since collections aren't always associated with a single strain, (e.g. "mixed" or "strain1_x_strain2") we don't set strain in processReadme().
     */
    void setStrain() {
        if (strain!=null) return;
        if (readme==null) {
            throw new RuntimeException("ERROR: attempted to setStrain() before README has been read.");
        }
        if (strainIdentifier==null) {
            throw new RuntimeException("ERROR: could not extract strain identifier from "+readme.identifier+".");
        }
        strain = createItem("Strain");
        strain.setAttribute("identifier", strainIdentifier);
        strain.setReference("organism", organism);
    }

    /**
     * Store core collection items:
     *
     *   dataSource
     *   dataSet
     *   organism
     *   strain (if non-null)
     *   publication (if non-null)
     *
     * @throws RuntimeException if README not read.
     */
    void storeCollectionItems() throws ObjectStoreException {
        // bail if README not read
        if (readme==null) {
            throw new RuntimeException("README not read. Aborting.");
        }
        // store stuff
        store(dataSource);
        store(dataSet);
        if (organism != null) store(organism);
        if (strain != null) store(strain);
        if (publication != null) store(publication);
        if (authors.size() > 0) store(authors);
    }

    /**
     * Form an LIS "full-yuck" primary identifier from a supplied secondaryIdentifier.
     * primaryIdentifier = gensp.strain.assy.annot.secondaryIdentifier
     */
    String formPrimaryIdentifier(String secondaryIdentifier) {
        if (readme==null) {
            throw new RuntimeException("ERROR: cannot formPrimaryIdentifier() because README has not yet been read.");
        }
        if (assemblyVersion==null || annotationVersion==null) {
            throw new RuntimeException("ERROR: cannot formPrimaryIdentifier() because assemblyVersion="+assemblyVersion+" and annotationVersion="+annotationVersion);
        }
        setStrain();
        return gensp+"."+strainIdentifier+"."+assemblyVersion+"."+annotationVersion+"."+secondaryIdentifier;
    }

    /**
     * Return true if the provided identifier matches the current collection.
     * 0     1      2   3   4
     * gensp.strain.gnm.ann.Name matches strain.gnm.ann.KEY4
     * 0     1      2   3       
     * gensp.strain.gnm.Name matches strain.gnm.KEY4
     */
    boolean matchesCollection(String identifier) {
        if (readme==null) {
            throw new RuntimeException("ERROR: DatastoreFileConverter.matchesCollection(identifier) cannot run - README has not yet been read.");
        }
        setStrain();
        String[] fields = identifier.split("\\.");
        if (fields.length>4) {
            return identifier.startsWith(gensp+"."+strainIdentifier+"."+assemblyVersion+"."+annotationVersion);
        } else if (fields.length==4) {
            return identifier.startsWith(gensp+"."+strainIdentifier+"."+assemblyVersion);
        } else {
            return false;
        }
    }

    /**
     * Return true if the strain and assembly version match between identifier and collection.
     * 0     1      2    3                                       1      2    !   !
     * phavu.G19833.gnm1.TOG905303_749 does not match collection G19833.gnm1.mrk.PvCookUCDavis2009
     */
    boolean matchesStrainAndAssembly(String identifier) {
        if (readme==null) {
            throw new RuntimeException("ERROR: DatastoreFileConverter.matchesStrainAndAssembly(identifier) cannot run - README has not yet been read.");
        }
        setStrain();
        String[] fields = identifier.split("\\.");
        try {
            boolean matches = fields[1].equals(strainIdentifier) && fields[2].equals(assemblyVersion);
            return matches;
        } catch (Exception ex) {
            System.err.println("ERROR in DatastoreFileConverter.matchesStrainAndAssembly: ID="+identifier);
            throw new RuntimeException(ex);
        }
    }

    /**
     * Return true if the given primaryIdentifier is for a Chromosome based on chromosome_prefix in genome README.
     *
     * README:
     * chromosome_prefix: chr
     */
    public boolean isChromosome(String primaryIdentifier) {
        for (String prefix : chromosomePrefixes) {
            if (primaryIdentifier.startsWith(gensp+"."+strainIdentifier+"."+assemblyVersion+"."+prefix)) return true;
        }
        return false;
    }

    /**
     * Return true if the given primaryIdentifier is for a Supercontig based on supercontig_prefix in genome README.
     *
     * README:
     * supercontig_prefix: chr
     */
    public boolean isSupercontig(String primaryIdentifier) {
        for (String prefix : supercontigPrefixes) {
            if (primaryIdentifier.startsWith(gensp+"."+strainIdentifier+"."+assemblyVersion+"."+prefix)) return true;
        }
        return false;
    }
    
    /**
     * Get the taxon ID for a gensp string like "phavu".
     */
    public String getTaxonId(String gensp) {
        String taxonId = genspTaxonId.get(gensp);
        if (taxonId==null) {
            throw new RuntimeException("Taxon ID not available for gensp="+gensp);
        }
        return taxonId;
    }

    /**
     * Load the properties from the mine's organism properties file.
     */
    void loadOrganismProperties() {
        Properties orgProps = new Properties();
        try {
            // load the organism properties into maps
            InputStream orgPropsResource = getClass().getClassLoader().getResourceAsStream(ORGANISM_PROP_FILE);
            if (orgPropsResource == null) {
                throw new RuntimeException("Did not find organism properties file:"+ORGANISM_PROP_FILE);
            }
            orgProps.load(orgPropsResource);
            for (Object obj : orgProps.keySet()) {
                String key = (String) obj;
                String value = orgProps.getProperty(key);
                String[] dotparts = key.split("\\.");
                if (dotparts[0].equals("taxon")) {
                    String taxonId = dotparts[1];
                    if (dotparts[2].equals("genus")) {
                        taxonIdGenus.put(taxonId, value);
                    } else if (dotparts[2].equals("species")) {
                        taxonIdSpecies.put(taxonId, value);
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("Problem loading properties from:"+ORGANISM_PROP_FILE, e);
        }
        // map gensp, Genus_species to taxonId
        for (String taxonId : taxonIdGenus.keySet()) {
            String genus = taxonIdGenus.get(taxonId);
            String species = taxonIdSpecies.get(taxonId);
            String gensp = genus.substring(0,3).toLowerCase()+species.substring(0,2).toLowerCase();
            String genusSpecies = genus+"_"+species;
            genspTaxonId.put(gensp, taxonId);
            taxonIdGensp.put(taxonId, gensp);
            genusSpeciesTaxonId.put(genusSpecies, taxonId);
        }
        
    }
    
    /**
     * Create the instance publication using CrossRef data. If CrossRef doesn't have the publication, use the README.publication_title.
     */
    void createPublication() throws UnsupportedEncodingException, MalformedURLException, ParseException, IOException, ParserConfigurationException, SAXException {
        publication = createItem("Publication");
        publication.setAttribute("doi", readme.publication_doi);
        // query CrossRef entry from DOI
        WorksQuery wq = new WorksQuery(readme.publication_doi);
        if (wq.getStatus() != null && wq.getStatus().equals("ok")) {
            String title = wq.getTitle();
            int year = 0;
            try { year = wq.getJournalIssueYear(); } catch (Exception ex) { }
            if (year==0) {
                try { year = wq.getIssuedYear(); } catch (Exception ex) { }
            }
            String month = null;
            try { month = String.valueOf(wq.getJournalIssueMonth()); } catch (Exception ex) { }
            if (month==null) {
                try { month = String.valueOf(wq.getIssuedMonth()); } catch (Exception ex) { }
            }
            String journal = null;
            if (wq.getShortContainerTitle()!=null) {
                journal = wq.getShortContainerTitle();
            } else if (wq.getContainerTitle()!=null) {
                journal = wq.getContainerTitle();
            }
            String volume = wq.getVolume();
            String issue = wq.getIssue();
            String pages = wq.getPage();
            String doi = wq.getDOI();
            if (!doi.equals(readme.publication_doi)) {
                // use CrossRef DOI
                System.err.println("### WARNING: CrossRef DOI " + doi + " does not match README.publication_doi " + readme.publication_doi);
                publication.setAttribute("doi", doi);
            }
            JSONArray authorsJSON = wq.getAuthors();
            String firstAuthor = null;
            if (authorsJSON!=null && authorsJSON.size()>0) {
                JSONObject firstAuthorObject = (JSONObject) authorsJSON.get(0);
                firstAuthor = (String) firstAuthorObject.get("family");
                if (firstAuthorObject.get("given")!=null) firstAuthor += ", " + (String) firstAuthorObject.get("given");
            }
            // get PubMed ID from PubMed API
            int pubMedId = DatastoreUtils.getPubMedId(doi);
            // core IM model does not contain lastAuthor
            // if (authorsJSON.size()>1) {
            //     JSONObject lastAuthorObject = (JSONObject) authorsJSON.get(authorsJSON.size()-1);
            //     lastAuthor = lastAuthorObject.get("family")+", "+lastAuthorObject.get("given");
            // }
            // update publication object
            if (title!=null) publication.setAttribute("title", title);
            if (firstAuthor!=null) publication.setAttribute("firstAuthor", firstAuthor);
            if (month!=null && !month.equals("0")) publication.setAttribute("month", month);
            if (journal!=null) publication.setAttribute("journal", journal);
            if (volume!=null) publication.setAttribute("volume", volume);
            if (issue!=null) publication.setAttribute("issue", issue);
            if (pages!=null) publication.setAttribute("pages", pages);
            if (year>0) publication.setAttribute("year", String.valueOf(year));
            if (pubMedId>0) publication.setAttribute("pubMedId", String.valueOf(pubMedId));
            // core IM model does not contain lastAuthor
            // if (lastAuthor!=null) publication.setAttribute("lastAuthor", lastAuthor);
                
            // populate publication.authors
            if (authorsJSON!=null) {
                for (Object authorObject : authorsJSON)  {
                    JSONObject authorJSON = (JSONObject) authorObject;
                    // IM Author attributes from CrossRef fields
                    String firstName = null; // there are rare occasions when firstName is missing, so we'll fill that in with a placeholder "X"
                    if (authorJSON.get("given")==null) {
                        firstName = "X";
                    } else {
                        firstName = (String) authorJSON.get("given");
                    }
                    // we require lastName, so if it's missing then bail on this author
                    if (authorJSON.get("family")==null) continue;
                    String lastName = (String) authorJSON.get("family");
                    // split out initials if present
                    // R. K. => R K
                    // R.K.  => R K
                    // Douglas R => Douglas R
                    // Douglas R. => Douglas R
                    String initials = null;
                    // deal with space
                    String[] parts = firstName.split(" ");
                    if (parts.length==2) {
                        if (parts[1].length()==1) {
                            firstName = parts[0];
                            initials = parts[1];
                        } else if (parts[1].length()==2 && parts[1].endsWith(".")) {
                            firstName = parts[0];
                            initials = parts[1].substring(0,1);
                        }
                    }
                    // pull initial out if it's an R.K. style first name (but not M.V.K.)
                    if (initials==null && firstName.length()==4 && firstName.charAt(1)=='.' && firstName.charAt(3)=='.') {
                        initials = String.valueOf(firstName.charAt(2));
                        firstName = String.valueOf(firstName.charAt(0));
                    }
                    // remove trailing period from a remaining R. style first name
                    if (firstName.length()==2 && firstName.charAt(1)=='.') {
                        firstName = String.valueOf(firstName.charAt(0));
                    }
                    String name = firstName+" "+lastName;
                    Item author = createItem("Author");
                    author.setAttribute("firstName", firstName);
                    author.setAttribute("lastName", lastName);
                    author.setAttribute("name", name);
                    if (initials!=null) author.setAttribute("initials", initials);
                    author.addToCollection("publications", publication);
                    authors.add(author);
                }
            }
        } else {
            // bail on WorksQuery and populate publication from README
            System.err.println("### DOI: " + readme.publication_doi + "  WorksQuery status: " + wq.getStatus());
            if (readme.publication_title != null) publication.setAttribute("title", readme.publication_title);
        }
    }

    /**
     * Grab chromosome/supercontig prefixes from corresponding genome collection Strain.gnm.KEY4.
     *
     * @param currentFile the current README file being processed (e.g. in an annotation collection)
     */
    void processGenomeReadme(File currentFile) throws IOException {
        File strainDir = currentFile.getParentFile().getParentFile().getParentFile();
        File genomesDir = new File(strainDir, "genomes");
        String[] genomes = genomesDir.list();
        String genomePrefix = strainIdentifier + "." + assemblyVersion;
        for (String genome : genomes) {
            if (genome.startsWith(genomePrefix)) {
                File genomeDir = new File(genomesDir, genome);
                String genomeReadmeFilename = "README." + genome + ".yml";
                File genomeReadmeFile = new File(genomeDir, genomeReadmeFilename);
                System.out.println("## Processing " + genomeReadmeFilename);
                Readme genomeReadme = Readme.parse(genomeReadmeFile);
                if (genomeReadme.chromosome_prefix == null && genomeReadme.supercontig_prefix == null) {
                    throw new RuntimeException("Genome README "+genomeReadmeFilename+" has neither chromosome_prefix or supercontig_prefix.");
                }
                if (genomeReadme.chromosome_prefix == null) {
                    System.out.println("## Genome README.chromosome_prefix is missing.");
                } else {
                    for (String prefix : genomeReadme.chromosome_prefix.split(",")) {
                        chromosomePrefixes.add(prefix.trim());
                    }
                }
                if (genomeReadme.supercontig_prefix == null) {
                    System.out.println("## Genome README.supercontig_prefix is missing.");
                } else {
                    for (String prefix : genomeReadme.supercontig_prefix.split(",")) {
                        supercontigPrefixes.add(prefix.trim());
                    }
                }
            }
        }
    }
    
}
