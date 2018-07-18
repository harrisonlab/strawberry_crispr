
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G


QueryReads=$1
RefDb=$2
OtuType=$3
Prefix=$4
OutDir=$5
Threshold=$6

CurDir=$PWD
WorkDir=${TMPDIR}/quantification
mkdir -p $WorkDir
cd $WorkDir

# quantify taxa
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -otutab $CurDir/$QueryReads -db $CurDir/$RefDb -strand plus -id 0.97 -biomout ${Prefix}_${OtuType}_table.biom -otutabout ${Prefix}_${OtuType}_table.txt -notmatched ${Prefix}_${OtuType}_nomatch.fa -userout ${Prefix}_${OtuType}_hits.out -userfields query+target
# Normalise results to 10000 reads
$ProgDir/usearch -otutab_norm ${Prefix}_${OtuType}_table.txt -sample_size 10000 -output ${Prefix}_${OtuType}_table_norm.txt

# Combine reads by species and filter taxa with reads fewer than a given threshold
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium_ampseq/scripts
$ProgDir/filter_OTU_tables.py --table ${Prefix}_${OtuType}_table.txt --threshold 0 > ${Prefix}_${OtuType}_table_by_spp.txt
$ProgDir/filter_OTU_tables.py --table ${Prefix}_${OtuType}_table.txt --threshold 206 > ${Prefix}_${OtuType}_table_by_spp_thresholded.txt
# Normalise these filtered samples
cat ${Prefix}_${OtuType}_table_by_spp_thresholded.txt | sed "s/^Species/#OTU ID/g" > tmp.txt
ProgDir=/home/deakig/usr/local/bin
$ProgDir/usearch -otutab_norm tmp.txt -sample_size 10000 -output tmp2.txt
cat tmp2.txt | sed "s/^#OTU ID/Species/g"  > ${Prefix}_${OtuType}_table_by_spp_thresholded_norm.txt

mkdir -p $CurDir/$OutDir
cp ${Prefix}_${OtuType}* $CurDir/$OutDir/.
