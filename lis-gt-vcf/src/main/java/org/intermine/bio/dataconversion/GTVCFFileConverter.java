package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2021 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 */

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.LinkedList;
import java.util.Set;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

import org.apache.log4j.Logger;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
/**
 * Load genotyping study data from a VCF file.
 *
 * @author Sam Hokin
 */
public class GTVCFFileConverter extends DatastoreFileConverter {
    private static final Logger LOG = Logger.getLogger(GTVCFFileConverter.class);

    // full classes
    Map<String,Item> genotypingStudies = new HashMap<>();
    Map<String,Item> chromosomes = new HashMap<>();
    Map<String,Item> samples = new HashMap<>();
    List<Item> publications = new LinkedList<>();

    /**
     * Constructor
     *
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public GTVCFFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     */
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().endsWith("vcf.gz")) {
            try {
                processVCF(reader);
            } catch (ObjectStoreException ex) {
                throw new RuntimeException(ex);
            }
        } else if (getCurrentFile().getName().startsWith("README")) {
            processREADME(reader);
        }
    }

    /**
     * Process the README, which contains the GenotypingStudy attributes.
     *
     * README.0SZD.yml
     * ---------------
     * identifier: 0SZD
     * subject: VCF file containing genotype information for 31 wild and cultivated soybeans received from Meng Ni, Tin Hang, and Hon-Ming Lam.
     * description: VCF file from resequencing 31 wild and cultivated chinese soybean accessions ....
     * genbank_accession: SRA020131
     * publication_doi: 10.1038/ng.715
     */
    void processREADME(Reader reader) throws IOException {
        // <attribute name="identifier" type="java.lang.String"/>
        // <attribute name="subject" type="java.lang.String"/>
        // <attribute name="description" type="java.lang.String"/>
        // <attribute name="genbank" type=="java.lang.String"/>
        // <reference name="publication" referenced-type="Publication"/>
        Readme readme = Readme.getReadme(getCurrentFile());
        Item study = genotypingStudies.get(readme.identifier);
        if (study==null) {
            study = createItem("GenotypingStudy");
            study.setAttribute("identifier", readme.identifier);
            genotypingStudies.put(readme.identifier, study);
        }
        study.setAttribute("subject", readme.subject);
        study.setAttribute("description", readme.description);
        study.setAttribute("genbank", readme.genbank_accession);
        Item publication = createItem("Publication");
        publication.setAttribute("doi", readme.publication_doi);
        publications.add(publication);
        study.setReference("publication", publication);
    }

    /**
     * Process the VCF file.
     *
     * 0     1    2    3   KEY4=identifier
     * glyma.Wm82.gnm2.div.0SZD.SNPData.vcf.gz
     */
    void processVCF(Reader reader) throws IOException, ObjectStoreException {
        // Items from the file name
        Item dataSet = getDataSet();
        Item organism = getOrganism();
        Item strain = getStrain(organism);
        String key4 = DatastoreUtils.extractKEY4(getCurrentFile().getName());
        Item study = genotypingStudies.get(key4);
        if (study==null) {
            study = createItem("GenotypingStudy");
            study.setAttribute("identifier", key4);
            genotypingStudies.put(key4, study);
        }
        study.setReference("organism", organism);
        study.setReference("strain", strain);
        study.setReference("dataSet", dataSet);
        // spin through the VCF file
        // NOTE: index must be present!
        VCFFileReader vcfReader = new VCFFileReader(getCurrentFile());
        VCFHeader header = vcfReader.getFileHeader();
        // load samples
        List<String> sampleNames = header.getGenotypeSamples();
        for (String sampleName : sampleNames) {
            Item sample = createItem("GenotypingSample");
            samples.put(sampleName, sample);
            study.addToCollection("samples", sample);
            sample.setAttribute("identifier", sampleName);
            sample.setReference("organism", organism);
            sample.setReference("strain", strain);
        }
        LOG.info("Loaded "+samples.size()+" samples from VCF header.");
        // VCFRecords and VCFSampleRecords
        // #CHROM  POS ID REF ALT QUAL FILTER INFO                       FORMAT  C01         C02        C08          C12          C14         ...
        // Chr01   5   .  T   G   97   .      DP=327;VDB=1.18178e-13;... GT:PL   0/0:0,12,85 0/0:0,9,61 0/0:0,25,126 0/0:0,60,170 0/0:0,3,124 ...
        for (VariantContext vc : vcfReader.iterator().toList()) {
            String identifier = vc.getID();
            String contig = vc.getContig();
            Allele ref = vc.getReference();
            // ALT - if there is more than one, store the first
            List<Allele> alts = vc.getAlternateAlleles();
            Allele alt = alts.get(0);
            // QUAL
            double qual = vc.getPhredScaledQual();
            // FILTER
            String filters = "";
            Set<String> filterSet = vc.getFilters();
            for (String filter : filterSet) {
                filters += filter;
            }
            // INFO data
            String info = "";
            Map<String,Object> attributes = vc.getCommonInfo().getAttributes();
            for (String key : attributes.keySet()) {
                info += key+"="+attributes.get(key).toString()+";";
            }
            // chromosome and location
            Item chromosome = chromosomes.get(vc.getContig());
            if (chromosome==null) {
                chromosome = createItem("Chromosome");
                chromosomes.put(vc.getContig(), chromosome);
                chromosome.setAttribute("secondaryIdentifier", vc.getContig());
                chromosome.setReference("organism", organism);
                chromosome.setReference("strain", strain);
                chromosome.addToCollection("dataSets", dataSet);
            }
            Item location = createItem("Location");
            location.setAttribute("start", String.valueOf(vc.getStart()));
            location.setAttribute("end", String.valueOf(vc.getEnd()));
            location.setReference("locatedOn", chromosome);
            // location.setReference("feature", bioEntity);
            location.addToCollection("dataSets", dataSet);
            store(location);
            // VCFRecord
            Item vcfRecord = createItem("VCFRecord");
            if (identifier!=null && !identifier.equals(".")) vcfRecord.setAttribute("identifier", identifier);
            vcfRecord.setAttribute("ref", ref.getBaseString());
            vcfRecord.setAttribute("alt", alts.get(0).getBaseString());
            vcfRecord.setAttribute("qual", String.valueOf(qual));
            if (filters.length()>0) vcfRecord.setAttribute("filter", filters);
            vcfRecord.setAttribute("info", info);
            vcfRecord.setReference("genotypingStudy", study);
            vcfRecord.setReference("chromosome", chromosome);
            vcfRecord.setReference("location", location);
            store(vcfRecord);
            // VCFSampleRecord
            for (String sampleName : samples.keySet()) {
                Item sample = samples.get(sampleName);
                String genotype = vc.getGenotype(sampleName).getGenotypeString();
                String likelihoods = vc.getGenotype(sampleName).getLikelihoods().getAsString();
                Item vcfSampleRecord = createItem("VCFSampleRecord");
                vcfSampleRecord.setReference("sample", sample);
                vcfSampleRecord.setAttribute("genotype", genotype);
                vcfSampleRecord.setAttribute("likelihoods", likelihoods);
                vcfSampleRecord.setReference("vcfRecord", vcfRecord);
                store(vcfSampleRecord);
            }
        }
        vcfReader.close();
    }

    /**
     * {@inheritDoc}
     */
    public void close() throws ObjectStoreException {
        // DatastoreFileConverter
        store(dataSource);
        store(dataSets.values());
        store(organisms.values());
        store(strains.values());
        // local not already stored
        store(genotypingStudies.values());
        store(publications);
        store(samples.values());
        store(chromosomes.values());
    }
}
