module load boost
g++ test_mergeVCF.cpp -std=c++11 -lboost_system -lboost_thread -lboost_unit_test_framework -o utest
./utest
exit
g++ tests/test_FUIntroSeg.cpp -lboost_system -lboost_thread -lboost_unit_test_framework -o utest
./utest
