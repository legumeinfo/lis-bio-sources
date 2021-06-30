package org.intermine.bio.dataconversion;

import java.io.File;
import java.io.IOException;
import java.io.Reader;

import com.fasterxml.jackson.dataformat.yaml.YAMLMapper;

/**
 * Encapsulates an LIS README file to be read in using the Jackson YAML reader.
 *
 * @author Sam Hokin
 */
public class Readme {
    // general attributes
    public String identifier; 
    public String provenance; 
    public String source; 
    public String synopsis; 
    public String related_to; 
    public String scientific_name; 
    public String taxid; 
    public String scientific_name_abbrev; 
    public String[] genotype; 
    public String description; 
    public String bioproject;
    public String dataset_doi; 
    public String genbank_accession; 
    public String original_file_creation_date; 
    public String local_file_creation_date;
    public String publication_doi;
    public String dataset_release_date;
    public String publication_title;
    public String contributors;
    public String data_curators;
    public String public_access_level;
    public String license;
    public String keywords; 
    public String citations;
    // source-specific attributes
    public String genotyping_platform;
    public String genotyping_method;
    public String expression_unit;
    public String sraproject;
    public String geoseries;

    /**
     * Create a Readme from a file.
     */
    public static Readme getReadme(File file) throws IOException {
        YAMLMapper mapper = new YAMLMapper();
        return mapper.readValue(file, Readme.class);
    }

    /**
     * Create a Readme from a Reader.
     */
    public static Readme getReadme(Reader reader) throws IOException {
        YAMLMapper mapper = new YAMLMapper();
        return mapper.readValue(reader, Readme.class);
    }
}

