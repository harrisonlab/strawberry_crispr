#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace


F=$1
R=$2
Locus=$3
Prefix=$4
OutDir=$5
# Adapters=$5
# MINL=150

CurDir=$PWD
WorkDir=${TMPDIR}/process_reads
mkdir -p $WorkDir
cd $WorkDir

# LABEL=${Prefix}

# ---
# Step1 merge reads and trim poor quality data
# ---
# The min length may need to be modified by locus.
# Maximum differences per 100bp have a large effect on how many bases are merged.
# The final number will still be limited by number of potential errors below.
MINL=150
MAXDIFF=5
# MAXDIFF=10
# LEN=$(( $MINL - $RPL ))
# Prefix='processing_test'
ProgDir=/home/deakig/usr/local/bin
# $ProgDir/usearch -fastq_mergepairs tmpF_ITS.fq -reverse tmpR_ITS.fq -fastqout ${Prefix}.t1  -fastq_pctid 0 -fastq_maxdiffs $(($MINL*${MAXDIFF}/100)) -fastq_minlen $MINL -fastq_minovlen 0 | tee 2>&1 ${Prefix}_mergelog.txt #-fastq_trunctail 25
$ProgDir/usearch -fastq_mergepairs $CurDir/$F -reverse $CurDir/$R -fastqout ${Prefix}.t1  -fastq_pctid 0 -fastq_maxdiffs $(($MINL*${MAXDIFF}/100)) -fastq_minlen $MINL -fastq_minovlen 0 | tee 2>&1 ${Prefix}_merge.log #-fastq_trunctail 25

# ---
# Step2 Identify reads containing illimuna adapters
# ---
Adapters=adapters.fa
printf \
">adapter1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>adapter2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" \
> $Adapters
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -search_oligodb ${Prefix}.t1 -db $Adapters -strand both -userout ${Prefix}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots
# ---
# Delete sequencing adapters adapters
# ---
# header lines of reads containing illumina sequencing adapters
# (identified in the previous step) are used to filter out reads.
ProgDir=/home/armita/git_repos/emr_repos/scripts/Metabarcoding_pipeline/scripts
awk -F"\t" '{print $1}' ${Prefix}.t1.txt|sort|uniq|$ProgDir/adapt_delete.pl ${Prefix}.t1 > ${Prefix}.t2

# ---
# Trim adapters off of fastq reads
# ---
# As the reads have come through the demultiplexing pipline which ensures that
# the reads start with PCR primers with a 100% match, we can simply trim the
# correspopnding number of bp off of the beginning and end of each joined sequence.
# FPL is forward primer length
# RPL is reverse primer length
Primers=primers.fa
printf \
">PDS74_F
TGAAGGGCTATTAGAAAATG
>PDS74_R
GGAGTTTTGGAATGGTTAAA
" \
> $Primers
# Locus=$(basename $F | cut -f1 -d '-')
FPL=$(cat $Primers | grep -A1 "${Locus}_F" | tail -n1 | wc -c)
RPL=$(cat $Primers | grep -A1 "${Locus}_R" | tail -n1 | wc -c)
awk  -v SL="$FPL" -v SR="$RPL" -F" " '{if(NR % 2 == 0){$1=substr($1,(SL+1),(length($1)-SL-SR))};print $1}' ${Prefix}.t2 > ${Prefix}.t3

# ---
# Identify reads which may contain errors based upon quality scores
# ---
# Number of expected errors in the fastq total read
# It does this by summing all the fastq quality scores for the read
# This step will also rename reads
QUAL=0.5
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -fastq_filter ${Prefix}.t3 -fastq_maxee $QUAL -relabel $Prefix -fastaout ${Prefix}.t3.fa

# Filter these potentially error-containing reads
cat ${Prefix}.t3.fa | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | sed -e '1d' > ${Prefix}.filtered.fa


# Prepare the unfiltered reads for later quantification steps:
cat ${Prefix}.t3 | awk -v S="$Prefix" -F" " '{if(NR % 4 == 1){print ">" S "." count+1 ";"$1;count=count+1} if(NR % 4 == 2){print $1}}' > ${Prefix}_prefilter.fa

# ---
# Cleanup
# ---
mkdir -p $CurDir/$OutDir/merged
mkdir -p $CurDir/$OutDir/filtered
mkdir -p $CurDir/$OutDir/unfiltered

mv *.log $CurDir/$OutDir/.
mv ${Prefix}.filtered.fa $CurDir/$OutDir/filtered/.
mv ${Prefix}.t3 $CurDir/$OutDir/unfiltered/${Prefix}.unfiltered.fastq
mv ${Prefix}_prefilter.fa $CurDir/$OutDir/merged/.

rm ${Prefix}.t1.txt ${Prefix}.t1 ${Prefix}.t3.fa
