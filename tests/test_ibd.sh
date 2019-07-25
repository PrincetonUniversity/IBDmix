#!/usr/bin/bash

export popfile=/tigress/AKEY/akey_vol1/home/luchenuw/WenqingEmpiricaltest_chr22/code/popfile.txt
export GT_HEAD=/tigress/tcomi/ibdmix_temp/alt_head.txt

g++ -std=c++11 IBDmix/IBD_Collection.cpp \
   IBDmix/Genotype_Reader.cpp  \
   IBDmix/IBD_Segment.cpp  \
   IBDmix/IBD_Stack.cpp  \
   IBDmix/IBDmix.cpp -o tests/ibd_new

g++ -std=c++11 IBDmix/FUIntroSeg_v2.1.cpp -o tests/ibd_orig
#set -euo pipefail

check_one (){
    local count=$1
    NEW_OUT=/tigress/tcomi/ibdmix_temp/new_ibd_${count}.txt
    pop=`head -n $count $popfile | tail -n 1`
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        chrom "gt lines" "in mask" "fail maf" "sel==0" \
        "lines recovered" "avg affected sites" "over estimate" > $pop.out
    INDIV=/tigress/AKEY/akey_vol2/wqfu/nobackup/1KGP/1KGP.${pop}.txt 
    echo "Starting on ${pop}"
    for chrnum in {1..22}; do
        GT_IN=/tigress/AKEY/akey_vol1/home/luchenuw/data/VCFs/Merged_Archaic_1KGP3/2013pubAltai/Altai_1KGP3_chr${chrnum}.txt
        OLD_OUT=/tigress/AKEY/akey_vol1/home/luchenuw/WenqingEmpiricaltest_chr22/Finalcallset/2013pubAltai/calling_intermediatefiles/Altai2013_1KGP3_chr${chrnum}_${pop}.txt 
        MASK=/tigress/AKEY/akey_vol1/home/luchenuw/data/FilterBed/FullVersionforAltai2013_1KGP3strict/Altai2013_1KGP3strict.mask.chr${chrnum}.bed

        printf "%d\t" $chrnum >> $pop.out
        tests/ibd_new -g <(cat $GT_HEAD $GT_IN) -o $NEW_OUT \
            -s <(awk '{if($2 == 1) print $1}' $INDIV) \
            -r $MASK -d 4 -e 0.002 >> $pop.out

        # convert new output to old, sorted
        awk -v OFS='\t' '{print $0, $5}' $NEW_OUT | sort > ${NEW_OUT}.sort
        line=$(cmp ${NEW_OUT}.sort <(sort $OLD_OUT) | awk '{print $NF}')

        if [ ! -z $line ]; then
            echo "MISMATCH! Line: $line Pop: $pop Chrom: $chrnum"
            awk -v line=$line 'NR==line{print "Expected: "$0; exit}' <(sort ${OLD_OUT})
            awk -v line=$line 'NR==line{print "Actual:   "$0; exit}' "${NEW_OUT}.sort"
            exit
        fi
    done
    echo "${pop} all matching!"
}

check_chrom (){
    local count=$1
    local chrnum=$2
    NEW_OUT=/tigress/tcomi/ibdmix_temp/new_ibd_${count}.txt
    pop=`head -n $count $popfile | tail -n 1`
    INDIV=/tigress/AKEY/akey_vol2/wqfu/nobackup/1KGP/1KGP.${pop}.txt 
    echo "Starting on ${pop} chrom ${chrnum}"
    GT_IN=/tigress/AKEY/akey_vol1/home/luchenuw/data/VCFs/Merged_Archaic_1KGP3/2013pubAltai/Altai_1KGP3_chr${chrnum}.txt
    OLD_OUT=/tigress/AKEY/akey_vol1/home/luchenuw/WenqingEmpiricaltest_chr22/Finalcallset/2013pubAltai/calling_intermediatefiles/Altai2013_1KGP3_chr${chrnum}_${pop}.txt 
    MASK=/tigress/AKEY/akey_vol1/home/luchenuw/data/FilterBed/FullVersionforAltai2013_1KGP3strict/Altai2013_1KGP3strict.mask.chr${chrnum}.bed

    time tests/ibd_new -g <(cat $GT_HEAD $GT_IN) -o $NEW_OUT \
        -s <(awk '{if($2 == 1) print $1}' $INDIV) \
        -r $MASK -d 4 -e 0.002

    exit
    # convert new output to old, sorted
    awk -v OFS='\t' '{print $0, $5}' $NEW_OUT | sort > ${NEW_OUT}.sort
    line=$(cmp ${NEW_OUT}.sort <(sort $OLD_OUT) | awk '{print $NF}')

    if [ ! -z $line ]; then
        echo "MISMATCH! Line: $line Pop: $pop Chrom: $chrnum"
        awk -v line=$line 'NR==line{print "Expected: "$0; exit}' <(sort ${OLD_OUT})
        awk -v line=$line 'NR==line{print "Actual:   "$0; exit}' "${NEW_OUT}.sort"
        exit
    fi
    echo "${pop} matching!"
}

check_old (){
    local count=$1
    local chrnum=$2
    NEW_OUT=/tigress/tcomi/ibdmix_temp/old_ibd_${count}.txt
    pop=`head -n $count $popfile | tail -n 1`
    INDIV=/tigress/AKEY/akey_vol2/wqfu/nobackup/1KGP/1KGP.${pop}.txt 
    echo "Starting on ${pop} chrom ${chrnum}"
    GT_IN=/tigress/AKEY/akey_vol1/home/luchenuw/data/VCFs/Merged_Archaic_1KGP3/2013pubAltai/Altai_1KGP3_chr${chrnum}.txt
    OLD_OUT=/tigress/AKEY/akey_vol1/home/luchenuw/WenqingEmpiricaltest_chr22/Finalcallset/2013pubAltai/calling_intermediatefiles/Altai2013_1KGP3_chr${chrnum}_${pop}.txt 
    MASK=/tigress/AKEY/akey_vol1/home/luchenuw/data/FilterBed/FullVersionforAltai2013_1KGP3strict/Altai2013_1KGP3strict.mask.chr${chrnum}.bed

    tests/ibd_orig -g $GT_IN -o $NEW_OUT \
        -s $INDIV \
        -i 1 -n $(wc -l $INDIV) -l $(wc -l $GT_IN) \
        -r $MASK -d 4 -e 0.002

    # convert new output to old, sorted
    line=$(cmp ${NEW_OUT} $OLD_OUT | awk '{print $NF}')

    if [ ! -z $line ]; then
        echo "MISMATCH! Line: $line Pop: $pop Chrom: $chrnum"
        awk -v line=$line 'NR==line{print "Expected: "$0; exit}' <(sort ${OLD_OUT})
        awk -v line=$line 'NR==line{print "Actual:   "$0; exit}' "${NEW_OUT}.sort"
        exit
    fi
    echo "${pop} matching!"
}
export -f check_one

check_one 1
# check_chrom 1 22
# seq 1 26 | xargs -P 1 -I {} bash -c 'check_one "{}"'
