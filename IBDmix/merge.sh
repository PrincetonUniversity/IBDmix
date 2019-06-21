#!/bin/bash

#SBATCH --time 24:00:00       # wall time
#SBATCH --mem=50G

#set -euox pipefail
date +%F" "%T


cd /tigress/AKEY/akey_vol1/home/luchenuw/SimonsProject/
module load samtools

vcftools --gzvcf /tigress/limingli/data/old-Altai/gz/chr${chrnum}.vcf.gz --out ./2013pubAltai/Altai_chr${chrnum} --remove-indels --min-alleles 1 --max-alleles 2 --recode

GZ_INPUT_FILE=/tigress/AKEY/akey_vol2/Reich_SGDP_2016/PS3_multisample_public/cteam_extended.v4.PS3_phase.public.chr${chrnum}.vcf.gz
VCF_INPUT_FILE=./2013pubAltai/Altai_chr${chrnum}.recode.vcf

tabix -h $GZ_INPUT_FILE ${chrnum}:1-250000000 | ./mergeVCF -a $VCF_INPUT_FILE -l $(($(wc -l <$VCF_INPUT_FILE) -135)) -d 1 -n 279 -o ./2013pubAltai/Altai2013pub_Simons_chr${chrnum}.txt

#for i in {1..22}; do echo $i; sbatch --export=chrnum=$i merge.sh; done
