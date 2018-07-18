#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

MergedSeqs=$1
OutDir=$2
Prefix=$3

CurDir=$PWD
WorkDir=${TMPDIR}/cluster_reads
mkdir -p $WorkDir
cd $WorkDir


#### Remove multiplex primers and pad reads to same length
X=`awk '{if ($1!~/>/){mylen=mylen+length($0)}else{print mylen;mylen=0};}' $CurDir/$MergedSeqs|awk '$0>x{x=$0};END{print x}'`
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -fastx_truncate $CurDir/$MergedSeqs -stripleft $SL -stripright $SR -trunclen $X -padlen $X -fastaout ${Prefix}_no_multiplex.fa
rm ${Prefix}.temp.fa

# Prefix=processing_test_ITS
# SL=0
# SR=0
# X=`awk '{if ($1!~/>/){mylen=mylen+length($0)}else{print mylen;mylen=0};}' processing_test.filtered.fa |awk '$0>x{x=$0};END{print x}'`
# ProgDir=/home/deakig/usr/local/bin
# $ProgDir/usearch -fastx_truncate processing_test.filtered.fa -stripleft $SL -stripright $SR -trunclen $X -padlen $X -fastaout ${Prefix}_no_multiplex.fa


#### Dereplication
ProgDir=/home/armita/git_repos/emr_repos/scripts/Metabarcoding_pipeline/scripts
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${Prefix}_no_multiplex.fa | $ProgDir/get_uniq.pl > ${Prefix}_dereplicated.fasta
#rm ${Prefix}.fa

#### Clustering (Cluster dereplicated seqeunces and produce OTU fasta (also filters for chimeras))

# ---
# Identify OTUs
# ---
# With 97% clustering, an OTU sequence should be at least 3% different from all
# other OTUs, and should be the most abundant sequences in its neighborhood.
# This is done by the cluster_otus command, which is an implementation of the
# UPARSE algorithm.
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -cluster_otus ${Prefix}_dereplicated.fasta -otus ${Prefix}_OTUs.fa -relabel OTU -minsize 4


# ---
# Identify zOTUs
# ---
# Denoising attempts to identify all correct biological sequences in the reads.
# This is done by the unoise3 command, which is an implementation of the UNOISE
# algorithm. A denoised sequence is called a "ZOTU" (zero-radius OTU).
# ZOTUs are valid OTUs for diversity analysis etc., though the interpretation of
# the results is a bit different from the usual 97% OTUs. For example, it is
# expected that one species may have more than one ZOTU, and with 97% OTUs it is
# expected than an OTU may have more than one species.
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -unoise3 ${Prefix}_dereplicated.fasta -zotus ${Prefix}_zOTUs_unparsed.fa #-relabel OTU #-minampsize 8
# usearch may need "zOUTs be names as OTUs in later steps."
cat ${Prefix}_zOTUs_unparsed.fa | sed -e 's/Zotu/OTU/' > $CurDir/$OutDir/${Prefix}_zOTUs.fa

#usearch -unoise ${Prefix}.sorted.fasta -tabbedout ${Prefix}.txt -fastaout ${Prefix}.otus.fa -relabel OTU #-minampsize 8

#perl -pi -e 's/uniq.*/OTU . ++$n/ge' ${Prefix}.otus.fa

#rm ${Prefix}.sorted.fasta
mv ${Prefix}_dereplicated.fasta $CurDir/$OutDir/.
mv ${Prefix}_OTUs.fa $CurDir/$OutDir/.
mv ${Prefix}_zOTUs_unparsed.fa $CurDir/$OutDir/.
