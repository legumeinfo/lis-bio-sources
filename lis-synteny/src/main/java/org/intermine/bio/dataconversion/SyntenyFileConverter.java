package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Read synteny blocks for two organisms from a DAGchainer synteny GFF file and store them as SyntenyBlock items,
 * each related to two SyntenicRegions. This is designed to use the GFF file produced by DAGchainer.
 *
 * ##gff-version 3
 * ##other comments
 * glyma.Wm82.gnm2.Gm01 DAGchainer syntenic_region 1122016 1208621 226.0 - . Name=phavu.G19833.gnm2.Chr02;matches=phavu.G19833.gnm2.Chr02:27981238..28077179;median_Ks=0.3600
 *
 * Source and target strains are given by the file name, for example: glyma.Wm82.gnm2.x.aradu.V14167.gnm1.gff
 *
 * @author Sam Hokin
 */
public class SyntenyFileConverter extends DatastoreFileConverter {
	
    private static final Logger LOG = Logger.getLogger(SyntenyFileConverter.class);

    static final String RANGE_SEPARATOR = "-";
    static final String REGION_SEPARATOR = "|";

    // stuff to store (that isn't declared in DatastoreFileConverter)
    List<Item> syntenyBlocks = new ArrayList<>();
    List<Item> sourceRegions = new ArrayList<>();
    List<Item> sourceChromosomeLocations = new ArrayList<>();
    List<Item> targetRegions = new ArrayList<>();
    List<Item> targetChromosomeLocations = new ArrayList<>();
    Map<String,Item> organisms = new HashMap<>();     // keyed by taxonId
    Map<String,Item> strains = new HashMap<>();       // keyed by strain identifier
    Map<String,Item> chromosomeMap = new HashMap<>(); // keyed by primaryIdentifier

    // keep to only one source/target pair
    Map<String,String> syntenyBlockIds = new HashMap<>();

    // for getting taxonId from gensp, etc.
    DatastoreUtils dsu;

    /**
     * Create a new SyntenyFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public SyntenyFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Process each GFF file by creating SyntenyBlock and SyntenicRegion items.
     * 0     1           2    3    4 5     6    7    8    (9)  (10)
     * cicar.CDCFrontier.gnm1.ann1.x.lotja.MG20.gnm3.ann1.7Bqh.gff3
     */
    @Override
    public void process(Reader reader) throws IOException {
        if (getCurrentFile().getName().startsWith("README")) {
            processReadme(reader);
            setStrain();
        } else if (getCurrentFile().getName().endsWith("gff") || getCurrentFile().getName().endsWith("gff3")) {
            processGFF(reader);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws ObjectStoreException {
        storeCollectionItems();
        store(organisms.values());
        store(strains.values());
	store(chromosomeMap.values());
	store(syntenyBlocks);
	store(sourceRegions);
	store(sourceChromosomeLocations);
	store(targetRegions);
	store(targetChromosomeLocations);
    }

    /**
     * Process a synteny GFF file.
     */
    public void processGFF(Reader reader) throws IOException {
        System.out.println("Processing "+getCurrentFile().getName());
        dsu = new DatastoreUtils();
        String[] fileNameParts = getCurrentFile().getName().split("\\.");
        if (fileNameParts.length<9) {
            System.err.println(getCurrentFile().getName()+" does not have at least 9 dot-separated parts.");
            System.err.println("Example: cicar.CDCFrontier.gnm1.ann1.x.lotja.MG20.gnm3.ann1.7Bqh.gff3");
            return;
        }
        // get the identifiers from the file name
        String sourceGensp = fileNameParts[0];
        String sourceStrainId = fileNameParts[1];
        String sourceAssy = fileNameParts[2];
        String sourceAnnot = fileNameParts[3];
        String targetGensp = fileNameParts[5];
        String targetStrainId = fileNameParts[6];
        String targetAssy = fileNameParts[7];
        String targetAnnot = fileNameParts[8];
        // get the organisms and strains
        Item sourceOrganism = getOrganism(sourceGensp);
        Item sourceStrain = getStrain(sourceStrainId, sourceOrganism);
        Item targetOrganism = getOrganism(targetGensp);
        Item targetStrain = getStrain(targetStrainId, targetOrganism);
        // -------------------------------------------------------------------------------------------------------
        // Load the GFF data into a map. Adds new chromosomes to chromosomeMap, keyed by primaryIdentifier.
        // -------------------------------------------------------------------------------------------------------
        BufferedReader gffReader = new BufferedReader(reader);
        String line = null;
        while ((line=gffReader.readLine())!=null) {
            // comment
            if (line.startsWith("#") || line.trim().length()==0) continue;
            // load the GFF line
            SyntenyGFF3Record gff = new SyntenyGFF3Record(line);
            if (gff.getType().equals("syntenic_region")) {
                String sourceChrName = getSourceChromosomeName(gff);
                String targetChrName = getTargetChromosomeName(gff);
                if (targetChrName==null) {
                    throw new RuntimeException("GFF syntenic_region record is missing target attribute:"+line);
                }
                // ignore this record if source or target are on a supercontig/scaffold
                if (dsu.isSupercontig(sourceChrName)) continue;
                if (dsu.isSupercontig(targetChrName)) continue;
                // get the source and target chromosomes
                Item sourceChromosome = getChromosome(sourceChrName, sourceOrganism, sourceStrain);
                Item targetChromosome = getChromosome(targetChrName, targetOrganism, targetStrain);
                // populate the source region and its location
                Item sourceRegion = createItem("SyntenicRegion");
                sourceRegion.setAttribute("assemblyVersion", sourceAssy);
                if (sourceAnnot!=null) sourceRegion.setAttribute("annotationVersion", sourceAnnot);
                Item sourceChromosomeLocation = createItem("Location");
                populateSourceRegion(sourceRegion, gff, sourceOrganism, sourceStrain, sourceChromosome, sourceChromosomeLocation);
                String sourceIdentifier = getSourceRegionName(gff);
                // populate the target region and its location
                Item targetRegion = createItem("SyntenicRegion");
                targetRegion.setAttribute("assemblyVersion", targetAssy);
                if (targetAnnot!=null) targetRegion.setAttribute("annotationVersion", targetAnnot);
                Item targetChromosomeLocation = createItem("Location");
                populateTargetRegion(targetRegion, gff, targetOrganism, targetStrain, targetChromosome, targetChromosomeLocation);
                String targetIdentifier = getTargetRegionName(gff);
                // only continue if we haven't stored a synteny block with these regions
                if (syntenyBlockIds.containsKey(sourceIdentifier) && syntenyBlockIds.get(sourceIdentifier).equals(targetIdentifier) ||
                    syntenyBlockIds.containsKey(targetIdentifier) && syntenyBlockIds.get(targetIdentifier).equals(sourceIdentifier)) {
                    // do nothing
                } else {
                    // store the regions in the map for future non-duplication
                    syntenyBlockIds.put(sourceIdentifier, targetIdentifier);
                    // get the medianKs value for this block
                    Map<String, List<String>> attributes = gff.getAttributes();
                    String medianKs = attributes.get("median_Ks").get(0);
                    // associate the two regions with this synteny block
                    Item syntenyBlock = createItem("SyntenyBlock");
                    syntenyBlock.setAttribute("medianKs", medianKs);
                    syntenyBlock.addToCollection("syntenicRegions", sourceRegion);
                    syntenyBlock.addToCollection("syntenicRegions", targetRegion);
                    syntenyBlocks.add(syntenyBlock);
                    // associate the block with the regions and store them
                    sourceRegion.setReference("syntenyBlock", syntenyBlock);
                    sourceRegions.add(sourceRegion);
                    sourceChromosomeLocations.add(sourceChromosomeLocation);
                    targetRegion.setReference("syntenyBlock", syntenyBlock);
                    targetRegions.add(targetRegion);
                    targetChromosomeLocations.add(targetChromosomeLocation);
                }
            }
        }
        gffReader.close();
    }
    
    /**
     * Populate the attributes of a source SyntenicRegion with a SyntenyGFF3Record's data.
     *
     * @param syntenicRegion the SyntenicRegion Item
     * @param gff the SyntenyGFF3Record holding the data
     * @param chromosome the source Chromosome Item
     * @param chromosomeLocation the source Location Item to be filled in
     */
    void populateSourceRegion(Item syntenicRegion, SyntenyGFF3Record gff, Item organism, Item strain, Item chromosome, Item chromosomeLocation) {
	syntenicRegion.setAttribute("primaryIdentifier", getSourceRegionName(gff));
	syntenicRegion.setAttribute("secondaryIdentifier", DatastoreUtils.extractSecondaryIdentifier(getSourceRegionName(gff), false));
	syntenicRegion.setAttribute("length", String.valueOf(getSourceEnd(gff)-getSourceStart(gff)+1));
	syntenicRegion.setAttribute("score", String.valueOf(gff.getScore()));
	syntenicRegion.setReference("organism", organism);
	syntenicRegion.setReference("strain", strain);
	syntenicRegion.setReference("chromosome", chromosome);
	syntenicRegion.setReference("chromosomeLocation", chromosomeLocation);
	chromosomeLocation.setAttribute("start", String.valueOf(getSourceStart(gff)));
	chromosomeLocation.setAttribute("end", String.valueOf(getSourceEnd(gff)));
	chromosomeLocation.setAttribute("strand", String.valueOf(gff.getStrand()));
	chromosomeLocation.setReference("feature", syntenicRegion);
	chromosomeLocation.setReference("locatedOn", chromosome);
    }

    /**
     * Populate the attributes of a target SyntenicRegion with a SyntenyGFF3Record's DAGchainer attributes data; Organism, Chromosome and ChromosomeLocation Items must be passed in as well.
     *
     * @param syntenicRegion the SyntenicRegion Item
     * @param gff the SyntenyGFF3Record holding the data
     * @param chromosome the target Chromosome Item
     * @param chromosomeLocation the target Location Item to be filled in
     */
    void populateTargetRegion(Item syntenicRegion, SyntenyGFF3Record gff, Item organism, Item strain, Item chromosome, Item chromosomeLocation) {
	syntenicRegion.setAttribute("primaryIdentifier", getTargetRegionName(gff));
	syntenicRegion.setAttribute("secondaryIdentifier", DatastoreUtils.extractSecondaryIdentifier(getTargetRegionName(gff), false));
	syntenicRegion.setAttribute("length", String.valueOf(getTargetEnd(gff)-getTargetStart(gff)+1));
	syntenicRegion.setAttribute("score", String.valueOf(gff.getScore()));
	syntenicRegion.setReference("organism", organism);
	syntenicRegion.setReference("strain", strain);
	syntenicRegion.setReference("chromosome", chromosome);
	syntenicRegion.setReference("chromosomeLocation", chromosomeLocation);
	chromosomeLocation.setAttribute("start", String.valueOf(getTargetStart(gff)));
	chromosomeLocation.setAttribute("end", String.valueOf(getTargetEnd(gff)));
	if (getTargetStrand(gff)!='\0') chromosomeLocation.setAttribute("strand", String.valueOf(getTargetStrand(gff)));
	chromosomeLocation.setReference("feature", syntenicRegion);
	chromosomeLocation.setReference("locatedOn", chromosome);
    }

    /**
     * Return the DAGchainer target strand from a DAGchainer Name attribute, checking for ' ' instead of '+' since SyntenyGFF3Record converts a plus to space.
     */
    char getTargetStrand(SyntenyGFF3Record gff) {
	String name = gff.getNames().get(0);
	char endChar = name.charAt(name.length()-1);
	if (endChar==' ' || endChar=='+') {
	    return '+';
	} else if (endChar=='-') {
	    return '-';
	} else {
	    return '\0';
	}
    }

    /**
     * Return the source chromosome nanme - opportunity for tweaks here
     */
    String getSourceChromosomeName(SyntenyGFF3Record gff) {
	return gff.getSequenceID();
    }

    /**
     * Return the DAGchainer target chromosome from a DAGchainer SyntenyGFF3Record
     * Target=Araip.B01:17125379..17229197
     */
    String getTargetChromosomeName(SyntenyGFF3Record gff) {
        String target = gff.getTarget();
        if (target==null) {
            System.err.println("GFF file has null target in the following record, so aborting:");
            System.err.println(gff.toString());
            System.exit(1);
        }
        String[] chunks = gff.getTarget().split(":");
        return chunks[0];
    }
    
    /**
     * Return the source sequence start from a GFF3 record - just echoes SyntenyGFF3Record.getStart() for clarity
     */
    int getSourceStart(SyntenyGFF3Record gff) {
	return gff.getStart();
    }

    /**
     * Return the source sequence end from a GFF3 record - just echoes SyntenyGFF3Record.getEnd() for clarity
     */
    int getSourceEnd(SyntenyGFF3Record gff) {
	return gff.getEnd();
    }

    
    /**
     * Return the target sequence start from a DAGchainer SyntenyGFF3Record
     * Target=Araip.B01:17125379..17229197
     */
    int getTargetStart(SyntenyGFF3Record gff) {
	String[] chunks = gff.getTarget().split(":");
	String range = chunks[1];
	String[] pieces = range.split("\\.\\.");
	return Integer.parseInt(pieces[0]);
    }

    /**
     * Return the target sequence end from a DAGchainer SyntenyGFF3Record
     */
    int getTargetEnd(SyntenyGFF3Record gff) {
	String[] chunks = gff.getTarget().split(":");
	String range = chunks[1];
	String[] pieces = range.split("\\.\\.");
	return Integer.parseInt(pieces[1]);
    }

    /**
     * Return a source syntenic region primary identifier from a GFF record; same as GBrowse standard
     */
    String getSourceRegionName(SyntenyGFF3Record gff) {
	String sourceChrName = getSourceChromosomeName(gff);
	int sourceStart = getSourceStart(gff);
	int sourceEnd = getSourceEnd(gff);
	return sourceChrName+":"+sourceStart+RANGE_SEPARATOR+sourceEnd;
    }

    /**
     * Return a target syntenic region primary identifier from a GFF record; same as GBrowse standard
     */
    String getTargetRegionName(SyntenyGFF3Record gff) {
	String targetChrName = getTargetChromosomeName(gff);
	int targetStart = getTargetStart(gff);
	int targetEnd = getTargetEnd(gff);
	return targetChrName+":"+targetStart+RANGE_SEPARATOR+targetEnd;
    }

    /**
     * Return a synteny block primary identifier formed from the source and target names with a separator
     */
    String getSyntenyBlockName(SyntenyGFF3Record gff) {
	return getSourceRegionName(gff)+REGION_SEPARATOR+getTargetRegionName(gff);
    }

    /**
     * Create/get a chromosome, extracting the secondaryIdentifier from the primaryIdentifier.
     */
    Item getChromosome(String primaryIdentifier, Item organism, Item strain) {
	if (chromosomeMap.containsKey(primaryIdentifier)) {
	    return chromosomeMap.get(primaryIdentifier);
	} else {
	    String secondaryIdentifier = DatastoreUtils.extractSecondaryIdentifier(primaryIdentifier, false);
	    if (secondaryIdentifier==null) {
		throw new RuntimeException("Error in chromosome primaryIdentifier:"+primaryIdentifier);
	    }
            String assemblyVersion = DatastoreUtils.extractAssemblyVersionFromFeature(primaryIdentifier);
            if (assemblyVersion==null) {
                throw new RuntimeException("Error in chromosome assemblyVersion:"+primaryIdentifier);
            }
	    Item chromosome = createItem("Chromosome");
	    chromosome.setAttribute("primaryIdentifier", primaryIdentifier);
	    chromosome.setAttribute("secondaryIdentifier", secondaryIdentifier);
            chromosome.setAttribute("assemblyVersion", assemblyVersion);
	    chromosome.setReference("organism", organism);
	    chromosome.setReference("strain", strain);
	    chromosomeMap.put(primaryIdentifier, chromosome);
	    return chromosome;
	}
    }

    /**
     * Return an organism Item associated with the given gensp.
     * There will be an organism associated with the collection, so don't duplicate that one.
     */
    Item getOrganism(String gs) {
        if (gs.equals(gensp)) {
            // return collection organism
            return organism;
        } else {
            // return a non-collection organism keyed by taxonId
            String taxonId = dsu.getTaxonId(gs);
            if (taxonId==null) {
                throw new RuntimeException("ERROR: could not get taxon ID from gensp="+gs);
            }
            if (organisms.containsKey(taxonId)) {
                return organisms.get(taxonId);
            } else {
                Item organism = createItem("Organism");
                organism.setAttribute("taxonId", taxonId);
                organisms.put(taxonId, organism);
                return organism;
            }
        }
    }

    /**
     * Return a Strain Item associated with the given identifier.
     * It may be the collection strain if it matches strainIdentifier, otherwise we create it here.
     */
    Item getStrain(String strainId, Item organism) {
        if (strainId.equals(strainIdentifier)) {
            return strain;
        } else if (strains.containsKey(strainId)) {
            return strains.get(strainId);
        } else {
            Item strain = createItem("Strain");
            strain.setAttribute("identifier", strainId);
            strain.setReference("organism", organism);
            strains.put(strainId, strain);
            return strain;
        }
    }
}
