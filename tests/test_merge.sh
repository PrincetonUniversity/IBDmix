#!/usr/bin/bash

MOD_INPUT_FILE=/tigress/AKEY/akey_vol2/wqfu/nobackup/1KGP/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
original_out=/tigress/AKEY/akey_vol1/home/luchenuw/data/VCFs/Merged_Archaic_1KGP3/2013pubAltai/Altai_1KGP3_chr1.txt
ARCH_INPUT_FILE=/tigress/limingli/data/old-Altai/gz/chr1.vcf.gz
CPP_OUT2=/tigress/tcomi/ibdmix_temp/alti_1kg.gz


#/usr/bin/time -v bash -c "cat $MOD_INPUT_FILE | ./mergeVCF -a $ARCH_INPUT_FILE -l $(wc -l <$ARCH_INPUT_FILE) -d 1 -n 279 -o $CPP_OUT"
#mv gmon.out gmon_old.out
g++ -std=c++11 IBDmix/generate_gt.cpp -o tests/generate_gt
tests/generate_gt \
    -a <(zcat $ARCH_INPUT_FILE | head -1000) \
    -m <(zcat $MOD_INPUT_FILE | head -1000) \
    -o - | gzip > $CPP_OUT2
line=$(cmp <(zcat $CPP_OUT2 | tail -n +2) $original_out | awk '{print $NF}')
if [ ! -z $line ]; then
    awk -v line=$line 'NR==line{print "Expected: "$1, $2, $3, $4, $5; exit}' "$original_out"
    zcat $CPP_OUT2 | awk -v line=$line 'NR==line+1{print "Actual: "$1, $2, $3, $4, $5; exit}'
else
    echo "Matching!"
fi
exit

EX=$CPP_OUT
ACT=${CPP_OUT2}.tail
tail -n +2 $CPP_OUT2 > $ACT

line=$(cmp "$ACT" "$EX" | awk '{print $NF}')
if [ ! -z $line ]; then
    awk -v file="$EX" -v line=$line 'NR==line{print "In file "file": "$1, $2, $3, $4, $5; exit}' "$EX"
    awk -v file="$ACT" -v line=$line 'NR==line{print "In file "file": "$1, $2, $3, $4, $5; exit}' "$ACT"
else
    echo "Matching!"
fi
