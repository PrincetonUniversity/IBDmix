module load boost
set -euo pipefail
g++ tests/test_mergeVCF.cpp -std=c++11 -lboost_system -lboost_thread -lboost_unit_test_framework -o tests/merge
tests/merge

g++ IBDmix/IBD_Stack.cpp tests/test_IBD_node.cpp \
    -std=c++11 \
    -o tests/node -lboost_system -lboost_thread -lboost_unit_test_framework
tests/node

g++ IBDmix/IBD_Stack.cpp IBDmix/IBD_Segment.cpp tests/test_IBD_segment.cpp \
    -std=c++11 \
    -o tests/segment -lboost_system -lboost_thread -lboost_unit_test_framework
tests/segment

g++ IBDmix/IBD_Stack.cpp IBDmix/IBD_Segment.cpp tests/test_IBD_segment_single.cpp \
    -std=c++11 \
    -o tests/single -lboost_system -lboost_thread -lboost_unit_test_framework
tests/single

g++ IBDmix/Genotype_Reader.cpp tests/test_Genotype_Reader.cpp \
    -std=c++11 \
    -o tests/reader -lboost_system -lboost_thread -lboost_unit_test_framework
tests/reader
