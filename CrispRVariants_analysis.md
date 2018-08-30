# CrispRVariants

Commands used to investigate crispr edits in lines of strawberry using the CrispRVariants package.


# Data transfer

Documented in the investigate-edits.md file

```bash
  ProjectDir=/home/groups/harrisonlab/project_files/strawberry_crispr
  cd $ProjectDir
```


#### QC of MiSeq data

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
  for RawData in $(ls raw_dna/paired/*/*/*/*/*/*.fastq.gz); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    echo $RawData;
    Jobs=$(qstat | grep 'run_fastq' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]; do
    sleep 1m
    printf "."
    Jobs=$(qstat | grep 'run_fastq' | grep 'qw' | wc -l)
    done		
    printf "\n"
    qsub $ProgDir/run_fastqc.sh $RawData
  done
```

# Align reads to the reference

```bash
for StrainPath in $(ls -d raw_dna/paired/*/*/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
done		
printf "\n"
ReadsF=$(ls $StrainPath/F/*.fastq*)
ReadsR=$(ls $StrainPath/R/*.fastq*)
echo $ReadsF
echo $ReadsR
OutDir=$(echo $StrainPath | sed 's/raw_dna/qc_dna/g')
qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA $OutDir
done
```

The number of reads for each treatment was determined:

```bash
  for F_read in $(ls qc_dna/paired/*/*/*/*/F/*fq.gz); do
    Locus=$(echo $F_read | cut -f3 -d '/')
    Pool=$(echo $F_read | cut -f4 -d '/')
    Dilution=$(echo $F_read | cut -f5 -d '/')
    TechRep=$(echo $F_read | cut -f6 -d '/')
    ReadCount=$(cat $F_read | gunzip -cf | awk '{s++}END{print s/4}')
    printf "$Locus\t$Pool\t$Dilution\t$TechRep\t$ReadCount\n"
  done > qc_dna/reads_per_sample.tsv
```


Move the reference file into the project directory
```bash
mkdir -p alignment/reference
cat ../fragaria_vesca/v4.0.a1/Fragaria_vesca_v4.0.a1.fasta | sed 's/^[^>]\s*$//g' | grep -v -E "^$" > alignment/reference/Fragaria_vesca_v4.0.a1.fasta
```


Run BWA-mem

```bash
CurDir=$PWD
Reference=$(ls alignment/reference/Fragaria_vesca_v4.0.a1.fasta )
for StrainPath in $(ls -d qc_dna/paired/*/*/*/* | tail -n+2); do
  Jobs=$(qstat | grep 'sub_bwa' | grep 'qw' | wc -l)
    while [ $Jobs -gt 2 ]; do
    sleep 20s
    printf "."
    Jobs=$(qstat | grep 'sub_bwa' | grep 'qw' | wc -l)
  done		
  printf "\n"
  ReadsF=$(ls $StrainPath/F/*fq.gz)
  ReadsR=$(ls $StrainPath/R/*fq.gz)
  Prefix=$(echo $ReadsF | rev | cut -f1 -d '/' | rev | cut -f1 -d '_')
  OutDir=alignment/bwa/vs_vesca/$(echo $StrainPath | cut -f3,4,5,6 -d '/')
  ProgDir=/home/armita/git_repos/emr_repos/scripts/strawberry_crispr/scripts
  qsub $ProgDir/sub_bwa_ampseq.sh $Prefix $CurDir/$Reference $ReadsF $ReadsR $OutDir
done
```

Download F.vesca gene models:

```bash
CurDir=$PWD
OutDir=alignment/reference
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.bioinfo.wsu.edu/www.rosaceae.org/Fragaria_vesca/Fvesca-genome.v4.0.a1/genes/Fragaria_vesca_v4.0.a1.transcripts.gff3.gz
gunzip *.gz
# Extract just the genes on the target contig:
cat Fragaria_vesca_v4.0.a1.transcripts.gff3 | grep -e 'Fvb4' -e 'gff-version 3' > Fragaria_vesca_v4.0.a1.transcripts_parsed.gff3
cd $CurDir
```


# Perform analysis using CrispRVariants

This is based upon commands found at:
https://bioconductor.org/packages/release/bioc/vignettes/CrispRVariants/inst/doc/user_guide.pdf

# Create a metadata file:
Create a file detailing experimental conditions for samples:

```bash
mkdir -p analysis/variants
printf "FileName\tLocation\tLocus\tCultivar\tLine\tSampleID\tPrefix\tGroup\n" > analysis/variants/metadata.tsv
# for Bam in $(ls alignment/bwa/vs_vesca/*/*/*/*/*_sorted.bam); do
# for Bam in $(ls alignment/bwa/vs_vesca/PDS74/Hawaii_4/*/*/*_sorted.bam | grep -w -e  'WT' -e 'Hawaii_4/10'); do
# for Bam in $(ls alignment/bwa/vs_vesca/PDS74/Calypso/*/*/*_sorted.bam); do
for Bam in $(ls alignment/bwa/vs_vesca/PDS74/Hawaii_4/*/*/*_sorted.bam); do
  Directory=$(dirname $Bam)
  FileName=$(basename $Bam)
  Locus=$(echo $Bam | cut -f4 -d '/')
  Cultivar=$(echo $Bam | cut -f5 -d '/')
  Line=$(echo $Bam | cut -f6 -d '/')
  SampleID=$(echo $Bam | cut -f7 -d '/')
  Prefix=$(echo $Bam | cut -f8 -d '/' | cut -f1 -d '_')
  Group=$(printf "${Cultivar} ${Line} ${SampleID}" | sed 's/Hawaii_4/Hawaii 4/g' | sed 's/Albino//g' | sed 's/Hawaii 4 WT 0/Hawaii 4 WT/g' | sed 's/FV//g' | sed 's/Hawaii 4 1 0/Hawaii 4 1 A/g' | sed 's/Calypso WT 0/Calypso WT/g')
  printf "$FileName\t$Bam\t$Locus\t$Cultivar\t$Line\t$SampleID\t$Prefix\t$Group\n"
done >> analysis/variants/metadata.tsv
```

Open R and start the analysis pipeline:
<!--
```R
library(CrispRVariants)
library("gdata")
setwd('/Users/armita/cluster_mount/groups/harrisonlab/project_files/strawberry_crispr')

#---
# Step 1 Load the metadata file
#---

# Read in the metadata file
md <- read.table("analysis/variants/metadata.tsv", sep = '\t',header = TRUE)
# Get the bam filenames from the metadata table
bam_fnames <- file.path(md$Location)
# check that all files exist
all( file.exists(bam_fnames))

#---
# Step 2 Create the target location and reference sequence
#---
# Note this bit was performed on the cluster where Samtools is installed
# cd /home/groups/harrisonlab/project_files/strawberry_crispr
# R
# setwd('/home/groups/harrisonlab/project_files/strawberry_crispr')
library(GenomicRanges)

df1 <- data.frame(chr="Fvb4", start=16333700, end=16333722,
                 strand=c("-"), score=0)
gd <- makeGRangesFromDataFrame(df1)
# Expande to detect variants within 10bp of the target site:
gdl <- GenomicRanges::resize(gd, width(gd) + 10, fix = "center")
# Retrieve the target sequence:
system("samtools faidx alignment/reference/Fragaria_vesca_v4.0.a1.fasta")
reference=system(sprintf("samtools faidx alignment/reference/Fragaria_vesca_v4.0.a1.fasta %s:%s-%s",
                         seqnames(gdl)[1], start(gdl)[1], end(gdl)[1]),
                 intern = TRUE)[[2]]
# The guide is on the negative strand, so the reference needs to be reverse comp
reference=Biostrings::reverseComplement(Biostrings::DNAString(reference))
save(reference, file = "alignment/reference/Fv4_PDS74.rda")


#---
# Step 2b Loading a previously saved reference sequence object
#---
# We’ll load the previously saved reference sequence.
# Note the NGG sequence (here, TGG) is present with the 5 extra bases on the end.
#ref_fname <- system.file(package="CrispRVariants", "/Users/armita/cluster_mount/groups/harrisonlab/project_files/strawberry_crispr/alignment/reference/Fv4_PDS74.rda")
#load(ref_fname)
load("alignment/reference/Fv4_PDS74.rda")
reference
##   33-letter "DNAString" instance
## seq: GGGATACCTGATCGAGTAACTACTGAGGTGTTT


#---
# Step 3 Creating a CriprSet
#---
# This step loads data into a crispr_set object for analysis

# in our case a local version of R studio is being used fro the analysis
# As the reference region was created on the cluster we need to reset some
# of the variabes used above:
df1 <- data.frame(chr="Fvb4", start=16333700, end=16333722, strand=c("-"), score=0)
gd <- makeGRangesFromDataFrame(df1)
gdl <- GenomicRanges::resize(gd, width(gd) + 10, fix = "center")

# Data can now be loaded into the object
# "Note that the zero point (target.loc parameter) is 22"
# ^ I dont know what this refers to!
crispr_set <- readsToTarget(bam_fnames, target = gdl, reference = reference,
                            names = md$Group, target.loc = 22)

# show the object:
crispr_set

# The counts table can be accessed with the "variantCounts" function
vc <- variantCounts(crispr_set)
print(class(vc))
nrow(vc)
# 1416

# You can see that in the table of variant counts, variants are summarised by the location of their insertions and deletions with respect to the target site. Non-variant sequences and sequences with a single nucleotide variant (SNV) but no insertion or deletion (indel) are displayed first, followed by the indel variants from most to least frequent For example, the most frequent non-wild-type variant, “-1:4D” is a 4 base pair deletion starting 1 base upstream of the zero point.

# Check the consensus sequence labelled as no variant matches the reference sequence:
sqs <- consensusSeqs(crispr_set)
Biostrings::reverseComplement(sqs[["no variant"]]) == reference

#---
# Step 4 Adding gene annotation
#---
library(GenomicFeatures)
gtf_fname <- "alignment/reference/Fragaria_vesca_v4.0.a1.transcripts_parsed.gff3"
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_fname, format = "gff")
saveDb(txdb, file= "alignment/reference/Fv4_PDS74_transcripts_txdb.sqlite")

#---
# Step 4 Creating summary plots of variants
#---
# We now load the the previously saved database
# plotVariants() is a wrapper function that groups together a plot of the transcripts of the gene/s overlapping the guide (optional), CrispRVariants::plotAlignments(), which dis- plays the alignments of the consensus variant sequences to the reference, and CrispRVari- ants::plotFreqHeatmap(), which produces a table of the variant counts per sample, coloured by either their counts or percentage contribution to the variants observed for a given sample. If a transcript database is supplied, the transcript plot is annotated with the guide location. Arguments for plotAlignments() and plotFreqHeatmap() can be passed to plotVariants() as lists named plotAlignments.args and plotFreqHeatmap.args, respectively.

# The gridExtra package is required to specify the legend.key.height
# as a "unit" object. It is not needed to call plotVariants() with defaults
library(gridExtra)

# Match the clutch id to the column names of the variants
group <- md$Group

# Set variable to colour all x labels black
ncolumns <- ncol(vc)
grp <- rep("1", each = ncolumns)


p <- plotVariants(crispr_set, txdb = txdb, gene.text.size = 8,
    row.ht.ratio = c(1,4), col.wdth.ratio = c(4,2),
    plotAlignments.args = list(line.weight = 0.5, ins.size = 2,
                               legend.symbol.size = 4
                               , min.freq = 10),
    plotFreqHeatmap.args = list(plot.text.size = 3, x.size = 8, group = group,
                                group.colours = grp, legend.text.size = 8,
                                legend.key.height = grid::unit(0.5, "lines")
                                , min.freq = 10, type = "counts", header = "counts"),
    left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines")
                              )
# plot_name <- "/Users/armita/Downloads/strawberry_crispr_plot1.pdf"
# plot_name <- "/Users/armita/Downloads/strawberry_crispr_Hawaii_4.pdf"
# plot_name <- "/Users/armita/Downloads/strawberry_crispr_Hawaii_4_WT.pdf"
plot_name <- "/Users/armita/Downloads/strawberry_crispr_Hawaii_4_line10.pdf"
ggsave(plot_name, plot = p, width =30, height = 15, units = "cm", limitsize = FALSE)
```
 -->



Run the R analysis above iterating through cultivar and lines

```R
library(CrispRVariants)
library("gdata")
library(gridExtra)
setwd('/Users/armita/cluster_mount/groups/harrisonlab/project_files/strawberry_crispr')

#---
# Step 1 Load the metadata file
#---

md_all <- read.table("analysis/variants/metadata.tsv", sep = '\t',header = TRUE)
md_all$analysis <- as.factor(paste(md_all$Cultivar, md_all$Line, sep = '_'))

for (t in levels(md_all$analysis)){
  md <- subset(md_all, md_all$analysis == t | md_all$analysis == 'Hawaii_4_WT' | md_all$analysis == 'Calypso_WT')
  crispr_function(md, t)
}



crispr_function <- function(md, prefix) {
  bam_fnames <- file.path(md$Location)
  all( file.exists(bam_fnames))
  #---
  # Step 2 Create the target location and reference sequence
  #---

  #---
  # Step 2b Loading a previously saved reference sequence object
  #---
  # We’ll load the previously saved reference sequence.
  # Note the NGG sequence (here, TGG) is present with the 5 extra bases on the end.
  load("alignment/reference/Fv4_PDS74.rda")
  reference

  #---
  # Step 3 Creating a CriprSet
  #---
  # This step loads data into a crispr_set object for analysis
  df1 <- data.frame(chr="Fvb4", start=16333700, end=16333722, strand=c("-"), score=0)
  gd <- makeGRangesFromDataFrame(df1)
  gdl <- GenomicRanges::resize(gd, width(gd) + 10, fix = "center")

  # Data can now be loaded into the object
  # "Note that the zero point (target.loc parameter) is 22"
  # ^ I dont know what this refers to!
  crispr_set <- readsToTarget(bam_fnames, target = gdl, reference = reference,
                              names = md$Group, target.loc = 22)


  # The counts table can be accessed with the "variantCounts" function
  vc <- variantCounts(crispr_set)
  print(class(vc))
  nrow(vc)

  #---
  # Step 4a Adding gene annotation
  #---
  library(GenomicFeatures)
  gtf_fname <- "alignment/reference/Fragaria_vesca_v4.0.a1.transcripts_parsed.gff3"
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_fname, format = "gff")
  saveDb(txdb, file= "alignment/reference/Fv4_PDS74_transcripts_txdb.sqlite")

  #---
  # Step 4b Creating summary plots of variants
  #---
  # We now load the the previously saved database

  # The gridExtra package is required to specify the legend.key.height
  # as a "unit" object. It is not needed to call plotVariants() with defaults

  # Match the clutch id to the column names of the variants
  group <- as.factor(md$Group)
  group <- droplevels(group)

  # Set variable to colour all x labels black
  ncolumns <- ncol(vc)
  grp <- rep("black", each = ncolumns)


  p <- plotVariants(crispr_set, txdb = txdb, gene.text.size = 6,
      row.ht.ratio = c(2,5), col.wdth.ratio = c(4,2),
      plotAlignments.args = list(line.weight = 1.0, ins.size = 2,
                                 legend.symbol.size = 4
                                 , min.freq = 10, plot.text.size = 4, axis.text.size = 12),
      plotFreqHeatmap.args = list(plot.text.size = 3, x.size = 12, x.angle = 45, group = group,
                                  group.colours = grp, legend.text.size = 10,
                                  legend.key.height = grid::unit(1.0, "lines")
                                  , min.freq = 10, type = "counts", header = "counts"),
      left.plot.margin = ggplot2::unit(c(0.1,0.1,5,0.5), "lines")
                                )

  w = 20 + (1.5 * ncolumns)
  h = 5 + (1 * ncolumns)

  plot_name <- paste("/Users/armita/Downloads/", t ,".pdf", sep = '')
  ggsave(plot_name, plot = p, width =w, height = h, units = "cm", limitsize = FALSE)
  plot_name <- paste("/Users/armita/Downloads/", t ,".eps", sep = '')
  ggsave(plot_name, plot = p, width =w, height = h, units = "cm", limitsize = FALSE)
  plot_name <- paste("/Users/armita/Downloads/", t ,".jpg", sep = '')
  ggsave(plot_name, plot = p, width =w, height = h, units = "cm", limitsize = FALSE)
}

```


## Make specific plots for figures in the paper:

```bash
mkdir -p analysis/variants
printf "FileName\tLocation\tLocus\tCultivar\tLine\tSampleID\tPrefix\tGroup\n" > analysis/variants/metadata_figs.tsv
for Bam in $(ls alignment/bwa/vs_vesca/PDS74/Calypso/*/*/*_sorted.bam | grep -e 'Calypso/1/' -e 'Calypso/5/' -e 'Calypso/94/' -e 'Calypso/WT/'); do
# for Bam in $(ls alignment/bwa/vs_vesca/PDS74/Hawaii_4/*/*/*_sorted.bam | grep -w -e '89' -e '39' -e 'WT'); do
  Directory=$(dirname $Bam)
  FileName=$(basename $Bam)
  Locus=$(echo $Bam | cut -f4 -d '/')
  Cultivar=$(echo $Bam | cut -f5 -d '/')
  Line=$(echo $Bam | cut -f6 -d '/')
  SampleID=$(echo $Bam | cut -f7 -d '/')
  Prefix=$(echo $Bam | cut -f8 -d '/' | cut -f1 -d '_')
  Group=$(printf "${Cultivar} ${Line} ${SampleID}" | sed 's/Hawaii_4/Hawaii 4/g' | sed 's/Albino//g' | sed 's/Hawaii 4 WT 0/Hawaii 4 WT/g' | sed 's/FV//g' | sed 's/Hawaii 4 1 0/Hawaii 4 1 A/g' | sed 's/Calypso WT 0/Calypso WT/g')
  printf "$FileName\t$Bam\t$Locus\t$Cultivar\t$Line\t$SampleID\t$Prefix\t$Group\n"
done >> analysis/variants/metadata_figs.tsv
```

Run the R analysis above iterating through cultivar and lines

```R
library(CrispRVariants)
library("gdata")
library(gridExtra)
setwd('/Users/armita/cluster_mount/groups/harrisonlab/project_files/strawberry_crispr')

#---
# Step 1 Load the metadata file
#---

md_all <- read.table("analysis/variants/metadata_figs.tsv", sep = '\t',header = TRUE)
md_all$analysis <- as.factor(paste(md_all$Cultivar, md_all$Line, sep = '_'))


md <- md_all
t <- "Calypso_FigX"
# t <- "Hawaii_4_FigX"

crispr_function(md, t)


crispr_function <- function(md, prefix) {
  bam_fnames <- file.path(md$Location)
  all( file.exists(bam_fnames))
  #---
  # Step 2 Create the target location and reference sequence
  #---

  #---
  # Step 2b Loading a previously saved reference sequence object
  #---
  # We’ll load the previously saved reference sequence.
  # Note the NGG sequence (here, TGG) is present with the 5 extra bases on the end.
  load("alignment/reference/Fv4_PDS74.rda")
  reference

  #---
  # Step 3 Creating a CriprSet
  #---
  # This step loads data into a crispr_set object for analysis
  df1 <- data.frame(chr="Fvb4", start=16333700, end=16333722, strand=c("-"), score=0)
  gd <- makeGRangesFromDataFrame(df1)
  gdl <- GenomicRanges::resize(gd, width(gd) + 10, fix = "center")

  # Data can now be loaded into the object
  # "Note that the zero point (target.loc parameter) is 22"
  # ^ I dont know what this refers to!
  crispr_set <- readsToTarget(bam_fnames, target = gdl, reference = reference,
                              names = md$Group, target.loc = 22)


  # The counts table can be accessed with the "variantCounts" function
  vc <- variantCounts(crispr_set)
  print(class(vc))
  nrow(vc)

  #---
  # Step 4a Adding gene annotation
  #---
  library(GenomicFeatures)
  gtf_fname <- "alignment/reference/Fragaria_vesca_v4.0.a1.transcripts_parsed.gff3"
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_fname, format = "gff")
  saveDb(txdb, file= "alignment/reference/Fv4_PDS74_transcripts_txdb.sqlite")

  #---
  # Step 4b Creating summary plots of variants
  #---
  # We now load the the previously saved database

  # The gridExtra package is required to specify the legend.key.height
  # as a "unit" object. It is not needed to call plotVariants() with defaults

  # Match the clutch id to the column names of the variants
  group <- as.factor(md$Group)
  group <- droplevels(group)

  # Set variable to colour all x labels black
  ncolumns <- ncol(vc)
  grp <- rep("black", each = ncolumns)


  p <- plotVariants(crispr_set, txdb = txdb, gene.text.size = 6,
      row.ht.ratio = c(0,5), col.wdth.ratio = c(4,2),
      plotAlignments.args = list(line.weight = 1.0, ins.size = 2,
                                 legend.symbol.size = 4
                                 , min.freq = 10, plot.text.size = 4, axis.text.size = 12),
      plotFreqHeatmap.args = list(plot.text.size = 3, x.size = 12, x.angle = 45, group = group,
                                  group.colours = grp, legend.text.size = 10,
                                  legend.key.height = grid::unit(1.0, "lines")
                                  , min.freq = 10, type = "counts", header = "counts"),
      left.plot.margin = ggplot2::unit(c(0.1,0.1,5,0.5), "lines")
                                )

  # w = 20 + (1.5 * ncolumns)
  # h = 5 + (1 * ncolumns)
  w = 20 + (1.7 * ncolumns)
  h = 5 + (0.8 * ncolumns)

  plot_name <- paste("/Users/armita/Downloads/", t ,".pdf", sep = '')
  ggsave(plot_name, plot = p, width =w, height = h, units = "cm", limitsize = FALSE)
  plot_name <- paste("/Users/armita/Downloads/", t ,".eps", sep = '')
  ggsave(plot_name, plot = p, width =w, height = h, units = "cm", limitsize = FALSE)
  plot_name <- paste("/Users/armita/Downloads/", t ,".jpg", sep = '')
  ggsave(plot_name, plot = p, width =w, height = h, units = "cm", limitsize = FALSE)
}
```
