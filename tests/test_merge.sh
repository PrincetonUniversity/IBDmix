#!/usr/bin/bash

MOD_INPUT_FILE=/tigress/AKEY/akey_vol2/wqfu/nobackup/1KGP/ALL.chr
MOD_SUFF=.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
original_out=/tigress/AKEY/akey_vol1/home/luchenuw/data/VCFs/Merged_Archaic_1KGP3/2013pubAltai/Altai_1KGP3_chr
ARCH_INPUT_FILE=/tigress/limingli/data/old-Altai/gz/chr
CPP_OUT2=/tigress/tcomi/ibdmix_temp/alti_1kg.gz


#/usr/bin/time -v bash -c "cat $MOD_INPUT_FILE | ./mergeVCF -a $ARCH_INPUT_FILE -l $(wc -l <$ARCH_INPUT_FILE) -d 1 -n 279 -o $CPP_OUT"
#mv gmon.out gmon_old.out
g++ -std=c++11 IBDmix/generate_gt.cpp -o tests/generate_gt

    #-a <(zcat $ARCH_INPUT_FILE | head -100000) \
for chr in {1..20}; do
    echo "starting on $chr"
    line=$(cmp <(tests/generate_gt \
        -a <(zcat ${ARCH_INPUT_FILE}${chr}.vcf.gz) \
        -m <(zcat $MOD_INPUT_FILE$chr$MOD_SUFF) \
        -o - | tee >(gzip > $CPP_OUT2) \
            | tail -n +2) \
        ${original_out}${chr}.txt | awk '{print $NF}')

    if [ ! -z $line ]; then
        echo "MISMATCH! Line: $line"
        awk -v line=$line 'NR==line{print "Expected: "$1, $2, $3, $4, $5; exit}' "${original_out}${chr}.txt"
        zcat $CPP_OUT2 | awk -v line=$line 'NR==line+1{print "Actual: "$1, $2, $3, $4, $5; exit}'
        exit
    else
        echo "Matching with chromosome ${chr}!"
    fi
done
