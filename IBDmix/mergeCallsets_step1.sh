#!/bin/bash

#SBATCH --time 2:00:00       # wall time
#SBATCH --array=1-26       # array indices and concurrency
#SBATCH --mem=20G


cd /tigress/AKEY/akey_vol1/home/luchenuw/WenqingEmpiricaltest_chr22/Finalcallset/2013pubAltai/
popfile=/tigress/AKEY/akey_vol1/home/luchenuw/WenqingEmpiricaltest_chr22/code/popfile.txt
pop=`head -n $SLURM_ARRAY_TASK_ID $popfile | tail -n 1`

module load anaconda

for i in {1..22}
do
        python /tigress/AKEY/akey_vol2/IntroSeg_WQFu/FUIntroSeg/FUIntroSeg_Summary.py -f Altai2013_1KGP3_chr${i}_${pop}.txt -l 30000 -d 4 -i 1
done

cat ALL_D4.0_L30000_Altai2013_1KGP3_chr*_$pop.txt | grep -v "start" | awk -v population="$pop" '{OFS="\t"} {$8=population; print $0}' > allchr.$pop.tmp



#cat MergedAmbig_Archaic0_D4.0_L30000_ADVmulti2016pub_1KGP3_chr*_$pop.txt | grep -v "start" | awk -v population="$pop" '{OFS="\t"} {$8=population; print $0}' > allchr.AltaiAmbig.$pop.tmp

#cat MergedAmbig_Archaic1_D4.0_L30000_ADVmulti2016pub_1KGP3_chr*_$pop.txt | grep -v "start" | awk -v population="$pop" '{OFS="\t"} {$8=population; print $0}' > allchr.VindiAmbig.$pop.tmp

#cat MergedAmbig_Archaic2_D4.0_L30000_ADVmulti2016pub_1KGP3_chr*_$pop.txt | grep -v "start" | awk -v population="$pop" '{OFS="\t"} {$8=population; print $0}' > allchr.DeniAmbig.$pop.tmp


