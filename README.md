# lis-bio-sources
Legume Information System data sources for building LIS mines from the LIS datastore.

#### lis-datastore
Abstract class extended by most LIS loaders. Contains method to read README files and includes BioStoreHook to associate DataSet.

#### lis-expression
Loads gene expression from LIS datastore /expression/ collections.

#### lis-fasta
Loads DNA or amino acid data from LIS FASTA files: chromosomes, supercontigs, proteins, mRNA, CDSes, etc.

#### lis-genefamily
Loads LIS gene family identifiers and descriptions (but not the members).

#### lis-genetic
Loads QTL/LG/marker data from LIS /genetic/ collections.

#### lis-gfa
Loads LIS gene family associations from a gfa.tsv file in an /annotations/ collection.

#### lis-gff
Loads GFF data from a GFF3 file in an LIS /annotations/ collection.

#### lis-gt-vcf
Loads VCF genotyping data from an LIS /diversity/ collection.

#### lis-gwas
Loads GWAS experiments from LIS /genetic/ collections.

### lis-info-annot
Loads info-annot.tsv files from an LIS /annotations/ collection.

#### lis-synteny
Loads synteny data from GFF files in an LIS /synteny/ collection.

#### lis-hsh
Loads data from an LIS datastore pan-gene set hsh.tsv file, e.g. glysp.mixed.pan2.TV81.hsh.tsv.

#### lis-ipr-gff
Loads an InterPro GFF file from an LIS /annotations/ collection.

#### lis-mstmap
Loads genotyping study data from an MSTmap file (from UC-Riverside, mostly).

#### lis-pathway
Loads gene pathway membership from an LIS /annotations/ collection.

#### lis-phylotree
Loads the LIS gene family phylotrees.

