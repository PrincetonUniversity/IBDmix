#!/bin/bash

#SBATCH --time 6:00:00       # wall time
#SBATCH --array=1-26       # array indices and concurrency
#SBATCH --mem=20G



cd /tigress/AKEY/akey_vol1/home/luchenuw/WenqingEmpiricaltest_chr22/Finalcallset/2013pubAltai
popfile=/tigress/AKEY/akey_vol1/home/luchenuw/WenqingEmpiricaltest_chr22/code/popfile.txt
pop=`head -n $SLURM_ARRAY_TASK_ID $popfile | tail -n 1`

gtfile=/tigress/AKEY/akey_vol1/home/luchenuw/data/VCFs/Merged_Archaic_1KGP3/2013pubAltai/Altai_1KGP3_chr${chrnum}.txt
bedfile=/tigress/AKEY/akey_vol1/home/luchenuw/data/FilterBed/FullVersionforAltai2013_1KGP3strict/Altai2013_1KGP3strict.mask.chr${chrnum}.bed

/tigress/AKEY/akey_vol2/IntroSeg_WQFu/FUIntroSeg/FUIntroSeg -g $gtfile -s /tigress/AKEY/akey_vol2/wqfu/nobackup/1KGP/1KGP.$pop.txt -o ./Altai2013pub_1KGP3_chr${chrnum}_$pop.txt -i 1 -n 2504 -l $(wc -l <$gtfile) -d 4 -e 0.002 -r $bedfile

