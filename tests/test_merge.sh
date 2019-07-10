MOD_INPUT_FILE=/tigress/tcomi/ibdmix_temp/mod.vcf
ARCH_INPUT_FILE=/tigress/tcomi/ibdmix_temp/arch.vcf
CPP_OUT=/tigress/tcomi/ibdmix_temp/cpp.txt
CPP_OUT2=/tigress/tcomi/ibdmix_temp/cpp2.txt
PY_OUT="/tigress/tcomi/ibdmix_temp/py_*.txt"
PY_OUT2=/tigress/tcomi/ibdmix_temp/py.txt

# module load anaconda3
# conda activate IBDmix
# 
# cd ../IBDmix
# 
# time python IBDmix.py merge-genotype \
#     --modern-vcf $MOD_INPUT_FILE \
#     --archaic-vcf $ARCH_INPUT_FILE \
#     --output "$PY_OUT"
# cd -
# # combine outputs from python, stripping header, moving columns, adding extra tab
# ls /tigress/tcomi/ibdmix_temp/py_*.txt | xargs tail -q -n +2 | \
#     awk -F'\t' '{ t = $1; $1 = $2; $2 = t; print $0, "";}' OFS=$'\t' \
#     > $PY_OUT2

#/usr/bin/time -v bash -c "cat $MOD_INPUT_FILE | ./mergeVCF -a $ARCH_INPUT_FILE -l $(wc -l <$ARCH_INPUT_FILE) -d 1 -n 279 -o $CPP_OUT"
#mv gmon.out gmon_old.out
g++ -std=c++11 IBDmix/generate_gt.cpp -o tests/generate_gt
/usr/bin/time -v bash -c "tests/generate_gt -a $ARCH_INPUT_FILE -m $MOD_INPUT_FILE -o $CPP_OUT2"
#valgrind tests/generate_gt -a $ARCH_INPUT_FILE -m $MOD_INPUT_FILE -o $CPP_OUT2


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
