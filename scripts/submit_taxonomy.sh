#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

OTUs=$1
RefDB=$2
Prefix=$3
OutDir=$4

CurDir=$PWD
WorkDir=${TMPDIR}/taxonomy
mkdir -p $WorkDir
cd $WorkDir

#### Assign Taxonomy
# ProgDir=/home/deakig/usr/local/bin
# $ProgDir/usearch -sintax clustering/ITS/ITS_OTUs.fa -db databases/ITS/utax_fungi_ITS_sintax.udp -strand both -tabbedout tmp/test.sintax
# ProgDir=/home/armita/git_repos/emr_repos/scripts/Metabarcoding_pipeline/scripts
# cat tmp/test.sintax | $ProgDir/mod_taxa_sintax.pl > tmp/test.taxa

ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -sintax $CurDir/$OTUs -db $CurDir/$RefDB -strand both -tabbedout ${Prefix}.sintax
ProgDir=/home/armita/git_repos/emr_repos/scripts/Metabarcoding_pipeline/scripts
cat ${Prefix}.sintax | $ProgDir/mod_taxa_sintax.pl > ${Prefix}.taxa

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
$ProgDir/add_taxa.py --fasta $CurDir/$OTUs --taxa ${Prefix}.taxa > ${Prefix}_taxa.fa

mkdir -p $CurDir/$OutDir
mv ${Prefix}.sintax $CurDir/$OutDir/.
mv ${Prefix}.taxa $CurDir/$OutDir/.
mv ${Prefix}_taxa.fa $CurDir/$OutDir/.
