#!/bin/bash

#SBATCH --time 6:00:00       # wall time
#SBATCH --mem=20G



cd /tigress/AKEY/akey_vol1/home/luchenuw/WenqingEmpiricaltest_chr22/Finalcallset/2013pubAltai/resampling

gtfile=/tigress/AKEY/akey_vol1/home/luchenuw/data/VCFs/Merged_Archaic_1KGP3/2013pubAltai/Altai_1KGP3_chr${chrnum}.txt
bedfile=/tigress/AKEY/akey_vol1/home/luchenuw/data/FilterBed/FullVersionforAltai2013_1KGP3strict/Altai2013_1KGP3strict.mask.chr${chrnum}.bed

/tigress/AKEY/akey_vol2/IntroSeg_WQFu/FUIntroSeg/FUIntroSeg -g $gtfile -s 1KGP.CEU10.txt -o ./Altai2013pub_1KGP3_chr${chrnum}_CEU10.txt -i 1 -n 2504 -l $(wc -l <$gtfile) -d 4 -e 0.002 -r $bedfile

