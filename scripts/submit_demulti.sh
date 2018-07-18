#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

R1=$1
R2=$2
OutDir=$3
# shift
# shift
# PrimerSets=$@

CurDir=$PWD
WorkDir=${TMPDIR}/demulti
mkdir -p $WorkDir
cd $WorkDir

F=$(basename $R1 | cut -f1 -d '.')
R=$(basename $R2 | cut -f1 -d '.')
# F=$(echo $R1|awk -F"/" '{print $NF}')
# R=$(echo $R2|awk -F"/" '{print $NF}')

# OUTDIR=$(echo $R1 |sed "s/$F//" | sed 's/raw_dna/demulti_dna/g')
echo "Output directory:"
echo $CurDir/${OutDir}
mkdir -p $CurDir/"$OutDir"


mkfifo $F.fa
mkfifo $R.fa

zcat -f -- $CurDir/$R1 > $F.fa &
zcat -f -- $CurDir/$R2 > $R.fa &



#PDS74
P1=PDS74
P1F=TGAAGGGCTATTAGAAAATG
P1R=GGAGTTTTGGAATGGTTAAA

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
$ProgDir/demulti.py --FastqF $F.fa --FastqR $R.fa --primer_loci $P1 --primersF $P1F --primersR $P1R 2>&1 | tee ${Prefix}_demulti.log

rm $F.fa $R.fa

cp *.fq $CurDir/$OutDir/.
cp *.fastq $CurDir/$OutDir/.
