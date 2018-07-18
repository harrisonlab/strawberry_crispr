# Invesigate-edits
Commands used to investigate crispr edits in lines of strawberry.

# Data transfer

Samples sequenced were:

```bash
  RawDatDir=/data/seq_data/miseq/2018/RAW/180629_M04465_0081_000000000-D4LBH/Data/Intensities/BaseCalls/
  ls $RawDatDir
```

Sample names referred to the following:

Loci:
PDS74


Raw sequencing data was symbolicly linked to the working directory:

```bash
  RawDatDir=/data/seq_data/miseq/2018/RAW/180629_M04465_0081_000000000-D4LBH/Data/Intensities/BaseCalls
  ProjectDir=/home/groups/harrisonlab/project_files/strawberry_crispr
  PlateDesign=$(ls raw_data/plate_design.txt)
  for ReadF in $(ls $RawDatDir/*_L001_R1_001.fastq.gz | grep -v 'Undetermined'); do
    ReadR=$(echo $ReadF | sed 's/_L001_R1_001.fastq.gz/_L001_R2_001.fastq.gz/g')
    Name=$(basename $ReadF _L001_R1_001.fastq.gz)
    Well=$(echo $Name | cut -f1 -d '_')
    Row=$(echo $Well | head -c 1)
    Column=$(echo $Well | tail -c +2)
    # Set Locus
    Locus=PDS74
    # Set cultivar
    Cultivar=$(cat $PlateDesign | grep -P "^${Row}\t${Column}" | cut -f3 | sed 's/ /_/g')
    # Set Line
    Line=$(cat $PlateDesign | grep -P "^${Row}\t${Column}" | cut -f4)
    # Set SampleID
    SampleID=$(cat $PlateDesign | grep -P "^${Row}\t${Column}" | cut -f5)
    printf "$Row\t$Column\t$Cultivar\t$Line\t$SampleID\n"
    OutDir=$ProjectDir/raw_dna/paired/$Locus/$Cultivar/$Line/$SampleID
    mkdir -p $OutDir/F
    cp -s $ReadF $OutDir/F/.
    mkdir -p $OutDir/R
    cp -s $ReadR $OutDir/R/.
  done
```

# Metabarcoding analysis

These commands are based upon the workflow described at:
https://github.com/eastmallingresearch/Metabarcoding_pipeline

### Set pipeline variables

```bash
# set MBPL variable to pipeline folder
MBPL=/home/deakig/metabarcoding_pipeline
```

# Analysis of individual loci

## Demultiplexing

This script demultiplexs mixed (e.g. ITS and 16S) libraries based on the primer sequence. Any sequence which has mismatches is written to ambiguous.fq (f & r seperately). Primer sequences
are detailed in the submission wrapper.
*Note* Regex are used to describe degenerate bases in the primer.

Run below to demultiplex:

```bash
  for DataDir in $(ls -d raw_dna/paired/*/*/*/*); do
    Jobs=$(qstat | grep 'submit_dem' | grep 'qw'| wc -l)
    while [ $Jobs -gt 5 ]; do
    sleep 1m
    printf "."
    Jobs=$(qstat | grep 'submit_dem' | grep 'qw'| wc -l)
    done
    printf "\n"
    WorkDir=/data/scratch/armita/fusarium_ampseq
    R1=$(ls $DataDir/F/*.fastq.gz)
    R2=$(ls $DataDir/R/*.fastq.gz)
    echo $DataDir
    echo $(basename $R1)
    echo $(basename $R2)
    ProgDir=/home/armita/git_repos/emr_repos/scripts/strawberry_crispr/scripts
    OutDir=demulti_dna/$(echo $R1 | rev | cut -f3,4,5,6 -d '/' | rev)
    qsub $ProgDir/submit_demulti.sh $R1 $R2 ${OutDir}
  done
```

 Summise reads demultiplexed:
```bash
printf "RunName\tLocus/Cultivar/Line/SampleID\tPDS74\tAmbiguous\n" > demulti_dna/demultiplex_summary.tsv
for RunDir in $(ls -d demulti_dna/*/*/*/*); do
  Locus=$(echo $RunDir | rev | cut -f4 -d '/' | rev)
  Cultivar=$(echo $RunDir | rev | cut -f3 -d '/' | rev)
  Line=$(echo $RunDir | rev | cut -f2 -d '/' | rev)
  SampleID=$(echo $RunDir | rev | cut -f1 -d '/' | rev)
  RunName=$(basename $RunDir/*_ambiguous.fq | sed 's/_ambiguous.fq//g')
  Pds74=$(ls $RunDir/*PDS74.fq)
  Pds74Lines=$(cat $Pds74 | wc -l)
  Ambiguous=$(ls $RunDir/*ambiguous.fq)
  AmbLines=$(cat $Ambiguous | wc -l)
  printf "$RunName\t$Locus\t$Pool\t$Dilution\t$Rep\t$Pds74Lines\t$AmbLines\n"
done >> demulti_dna/demultiplex_summary.tsv
```

<!-- From this data thresholding values of cross-contamination as a result of illumina
adapter read hopping were determined. For each row in the dataset, the most
abundant locus was assumed to be the target locus and reads attributed to other
loci assumed to be contaminant reads. Reads from another experiment with the same
locus could be contaminating the sample. The 2nd most abundant locus was
identified for each run (the highest contaminant locus) and the maximum value
identified accross the entire plate. A threshold for a minimum abundance of reads
attributed to an OTU was set at this value.

In the case of our plate, this was:
```bash
Threshold=206
``` -->


## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.


```bash
for DataDir in $(ls -d demulti_dna/*/*/*/*); do
  Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
  while [ $Jobs -gt 5 ]; do
  sleep 20s
  printf "."
  Jobs=$(qstat | grep 'sub_pro' | grep 'qw'| wc -l)
  done
  printf "\n"
  Locus=$(echo $DataDir | cut -f2 -d '/')
  Prefix=$(echo $DataDir | cut -f2,3,4,5 -d '/' | sed 's&/&_&g')
  WorkDir=/data/scratch/armita/fusarium_ampseq
  R1=$(ls $DataDir/*_${Locus}.fq | head -n1 | tail -n1)
  R2=$(ls $DataDir/*_${Locus}.fq | head -n2 | tail -n1)
  echo $DataDir
  echo $(basename $R1)
  echo $(basename $R2)
  OutDir="processed_dna/"$(echo $R1 | rev | cut -f2,3,4,5 -d '/' | rev)
  ProgDir=/home/armita/git_repos/emr_repos/scripts/strawberry_crispr/scripts
  qsub $ProgDir/sub_process_reads.sh $R1 $R2 $Locus $Prefix $OutDir
done
```


## OTU assignment
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. I've written my own scripts to do the dereplication and sorting

 * All read files associated with the each locus in the project are concatenated
 * Clustering run on all the data for this locus within the project
 * Quantification can then be performed for each treatment against the total set

```bash
  # Concatenate files
  for Locus in "PDS74"; do
    OutDir=clustering/$Locus
    mkdir -p $OutDir
    cat processed_dna/$Locus/*/*/*/filtered/*.filtered.fa > $OutDir/${Locus}_concatenated.temp.fa
    ls -lh $OutDir/${Locus}_concatenated.temp.fa
    # Submit clustering
    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
    qsub $ProgDir/submit_clustering.sh $OutDir/${Locus}_concatenated.temp.fa $OutDir $Locus
  done
```

<!-- # Submit Taxonomy -->



# Quantify variants

Create OTU tables

```bash
for Locus in "PDS74"; do
  for Cultivar in Hawaii_4 Calypso; do
  # for Pool in Fusarium_spp; do
    OutDir=quantified/$Locus/$Cultivar
    mkdir -p $OutDir
    cat processed_dna/$Locus/$Cultivar/*/*/merged/*.fa | cut -f1 -d '.' > $OutDir/${Locus}_reads_appended.fa
    QueryReads=$(ls $OutDir/${Locus}_reads_appended.fa)
    for OtuType in zOTUs; do
      RefDb=$(ls clustering/PDS74/PDS74_zOTUs.fa)
      Prefix=$Locus
      ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
      qsub $ProgDir/submit_quantification.sh $QueryReads $RefDb $OtuType $Prefix $OutDir
    done
  done
done
```
