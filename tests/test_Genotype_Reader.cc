#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include "IBDmix/Genotype_Reader.h"

class SampleGenotype : public ::testing::Test {
    protected:
        void SetUp() {
            std::ofstream genotype;
            genotype.open("test_genotype.txt");
            genotype << "chrom\tpos\tref\talt\tn1\tm1\tm2\tm3\tm4\n"
                << "1\t2\tA\tT\t1\t0\t0\t0\t0\n"
                << "1\t3\tA\tT\t2\t0\t0\t0\t0\n"
                << "1\t4\tA\tT\t1\t0\t1\t1\t1\n"
                << "1\t104\tA\tT\t1\t0\t1\t1\t1\n"
                << "1\t105\tA\tT\t0\t2\t1\t1\t1\n"
                << "2\t125\tA\tT\t0\t2\t2\t1\t1\n"
                << "3\t126\tA\tT\t0\t2\t2\t2\t2\n"
                ;
            genotype.close();
            std::ofstream mask;
            mask.open("test_mask.txt");
            mask << "1 100 120\n"
                << "1 130 140\n"
                << "1 160 161\n"
                << "1 190 200\n"
                << "1 260 281\n"
                << "2 130 140\n"
                << "4 130 140\n"
                << "8 130 140\n"
                ;
            mask.close();
        }
        void TearDown() {
            remove("test_genotype.txt");
            remove("test_mask.txt");
        }
};

TEST_F(SampleGenotype, CanInitialize){
    FILE * file = fopen("test_genotype.txt", "rt");
    // use explicit values in case defaults change
    Genotype_Reader reader = Genotype_Reader(
            file, nullptr, 0.01, 0.002, 2, 1e-200);
    ASSERT_EQ(4, reader.initialize(nullptr));
    rewind(file);
    std::ofstream samples;
    samples.open("test_samples.txt");
    samples << "m2\nm4\n";
    samples.close();
    ASSERT_EQ(2, reader.initialize(fopen("test_samples.txt", "rt")));
    rewind(file);
    samples.open("test_samples.txt");
    samples << "m2\nm4";
    samples.close();
    ASSERT_EQ(2, reader.initialize(fopen("test_samples.txt", "rt")));
    rewind(file);
    samples.open("test_samples.txt");
    samples << "m2";
    samples.close();
    ASSERT_EQ(1, reader.initialize(fopen("test_samples.txt", "rt")));
}

TEST(TestGenotypeReader, CanFindToken){
    ASSERT_EQ(2, find_token("tes", "test\tte\ttes"));
    ASSERT_EQ(-1, find_token("not in", "test\tstring\t"));
    ASSERT_EQ(0, find_token("test", "test\tstring\t"));
    ASSERT_EQ(1, find_token("string", "test\tstring\t"));
    ASSERT_EQ(1, find_token("string", "test\tstring"));
    ASSERT_EQ(1, find_token("string\tthing", "test\tstring"));
    ASSERT_EQ(1, find_token("string\tthing", "test\tstring\tthings"));
    ASSERT_EQ(-1, find_token("strinG", "test\tstring"));
    ASSERT_EQ(-1, find_token("string", "test\tstrinG"));
    ASSERT_EQ(-1, find_token("string1", "test\tstring"));
    ASSERT_EQ(-1, find_token("string", "test\tstring1"));
    ASSERT_EQ(-1, find_token("String", "test\tstring1"));
    ASSERT_EQ(-1, find_token("\0", "test\tstring1"));
}

TEST(TestGenotypeReader, CanFindArchaic){
    Genotype_Reader reader = Genotype_Reader(nullptr);
    reader.samples = new char[100];
    reader.buffer = new char[100];
    strcpy(reader.buffer, "n1\ts1\ts2");

    // test both null, delete first token in samples
    strcpy(reader.samples, "n1\ts1\ts2");
    reader.find_archaic(nullptr, nullptr);
    ASSERT_EQ(0, reader.archaic_index);
    ASSERT_STREQ(reader.samples,  "s1\ts2");

    // test with sample line, delete nothing
    strcpy(reader.samples, "n1\ts1\ts2");
    reader.find_archaic(nullptr, "");
    ASSERT_EQ(0, reader.archaic_index);
    ASSERT_STREQ(reader.samples,  "n1\ts1\ts2");

    // test with archaic
    strcpy(reader.samples, "n1\ts1\ts2");
    reader.find_archaic("s1", nullptr);
    ASSERT_EQ(1, reader.archaic_index);
    ASSERT_STREQ(reader.samples,  "n1\ts2");

    strcpy(reader.samples, "n1\ts1\ts2");
    reader.find_archaic("s2", nullptr);
    ASSERT_EQ(2, reader.archaic_index);
    ASSERT_STREQ(reader.samples,  "n1\ts1\t");

    strcpy(reader.samples, "n1\ts1\ts2\t");
    reader.find_archaic("s2", nullptr);
    ASSERT_EQ(2, reader.archaic_index);
    ASSERT_STREQ(reader.samples,  "n1\ts1\t");

    strcpy(reader.samples, "n1\ts1\ts2\t");
    reader.find_archaic("s2", "");
    ASSERT_EQ(2, reader.archaic_index);
    ASSERT_STREQ(reader.samples,  "n1\ts1\ts2\t");
}

TEST(TestGenotypeReader, CanDetermineSampleMapping){
    Genotype_Reader reader = Genotype_Reader(nullptr);
    reader.buffer = new char[100];
    strcpy(reader.buffer, "n1\ts1\ts2");
    reader.sample_to_index = new int[10];
    reader.num_samples = 10;

    reader.archaic_index = 0;
    reader.determine_sample_mapping(nullptr);
    for(int i = 0; i < reader.num_samples; i++)
        ASSERT_EQ(i+1, reader.sample_to_index[i]);

    reader.archaic_index = 10;
    reader.determine_sample_mapping(nullptr);
    for(int i = 0; i < reader.num_samples; i++)
        ASSERT_EQ(i, reader.sample_to_index[i]);

    reader.archaic_index = 5;
    reader.determine_sample_mapping(nullptr);
    for(int i = 0; i < 5; i++)
        ASSERT_EQ(i, reader.sample_to_index[i]);
    for(int i = 5; i < reader.num_samples; i++)
        ASSERT_EQ(i+1, reader.sample_to_index[i]);

    strcpy(reader.buffer, "n1\ts1\ts2\ts3\ts4\ts5\ts6\ts7");
    reader.num_samples = 3;
    reader.determine_sample_mapping("s6\tn1\ts3");
    int result[] = {6, 0, 3};
    for(int i = 0; i < reader.num_samples; i++)
        ASSERT_EQ(result[i], reader.sample_to_index[i]);

    reader.determine_sample_mapping("s1\ts1\ts1");
    int result2[] = {1, 1, 1};
    for(int i = 0; i < reader.num_samples; i++)
        ASSERT_EQ(result2[i], reader.sample_to_index[i]);
}

TEST(TestGenotypeReader, CanYieldSample){
    Genotype_Reader reader = Genotype_Reader(nullptr);
    reader.num_samples = 3;
    reader.samples = new char[100];
    strcpy(reader.samples, "n1\ts1\ts2");
    char *ptr = nullptr;
    int count = 0;

    ASSERT_TRUE(reader.yield_sample(ptr, count++));
    ASSERT_STREQ(ptr, "n1");

    ASSERT_TRUE(reader.yield_sample(ptr, count++));
    ASSERT_STREQ(ptr, "s1");

    ASSERT_TRUE(reader.yield_sample(ptr, count++));
    ASSERT_STREQ(ptr, "s2");

    ASSERT_TRUE(!reader.yield_sample(ptr, count++));
    ASSERT_STREQ(reader.samples, "n1\ts1\ts2");
}

double original_cal_lod(int source_gt, int target_gt, double pb, double aerr, double merr,
        double minesp){
    double pa = 1-pb;
    if(source_gt == 0 && target_gt == 0)
        return log10(((1-aerr)*(1-merr)+aerr*merr) / pa / (1-aerr*(1-aerr)));
    if(source_gt == 0 && target_gt == 1){
        double temp = 0.5 * ((1-aerr)*merr+aerr*(1-merr)) / pb /
                (1-aerr*(1-aerr)) + 0.5 * ((1-aerr)*(1-merr)+aerr*merr) /
                pa / (1-aerr*(1-aerr));
        return log10(temp);
    }
    if(source_gt == 0 && target_gt == 2)
        return log10(((1-aerr)*merr+aerr*(1-merr) + minesp) /
                pb / (1-aerr*(1-aerr) + minesp));
    if(source_gt == 1 && target_gt == 0)
        return -log10(pa * (1+2*aerr*(1-aerr)));
    if(source_gt == 1 && target_gt == 1)
        return -log10(2*pa*pb*(1+2*aerr*(1-aerr)));
    if(source_gt == 1 && target_gt == 2)
        return -log10(pb * (1+2*aerr*(1-aerr)));
    if(source_gt == 2 && target_gt == 0)
        return log10(((1-aerr)*merr+aerr*(1-merr) + minesp) /
                pa / (1-aerr*(1-aerr) + minesp));
    if(source_gt == 2 && target_gt == 1)
        return log10(0.5*((1-aerr)*(1-merr)+aerr*merr) / pb /
                (1-aerr*(1-aerr)) + 0.5 * ((1-aerr)*merr+aerr*(1-merr)) /
                pa / (1-aerr*(1-aerr)));
    if(source_gt == 2 && target_gt == 2)
        return log10(((1-aerr)*(1-merr)+aerr*merr) / pb / (1-aerr*(1-aerr)));
    return 0;
}

TEST(TestGenotypeReader, CanCalculateLod){
    Genotype_Reader reader = Genotype_Reader(nullptr);

    double pbs[] = {0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.99, 1};
    double aerrs[] = {0.01, 0.02, 0.03};
    double merrs[] = {0.0001, 0.001, 0.002};
    double minesps[] = {1e-200, 1e-190, 1e-100};
    char gts[] = {'0', '1', '2', '9'};
    for(double pb : pbs) for(double aerr : aerrs) for(double merr : merrs)
        for(double minesp : minesps) for(char arch : gts){
            reader.archaic_error = aerr;
            reader.minesp = minesp;
            reader.update_lod_cache(arch, pb, merr);
            for(char mod : gts){
                double read = reader.calculate_lod(mod);
                double orig = original_cal_lod(arch-'0', mod-'0', pb, aerr, merr, minesp);
                if(isinf(read) && isinf(orig))
                    continue;
                ASSERT_TRUE(read == orig || abs((read - orig)/orig) < 0.001);
            }
            reader.update_lod_cache(arch, pb, merr, false); // not selected
            for(char mod : gts){
                double read = reader.calculate_lod(mod);
                double orig = original_cal_lod(arch-'0', mod-'0', pb, aerr, merr, minesp);
                if((mod == '0' && arch == '2') || (mod == '2' && arch == '0')){
                    if(isinf(read) && isinf(orig))
                        continue;
                    //std::cout << "mod=" << mod-'0' << " arch=" << arch-'0' << 
                    //        " pb=" << pb << " aerr=" << aerr << " merr=" << merr
                    //        << " minesp=" << minesp << "\n"
                    //        << "reader=" << read
                    //        << "\norig=  " << orig
                    //        << "\n";
                    ASSERT_DOUBLE_EQ(read, orig);
                }
                else{
                    ASSERT_EQ(0, read);
                }
            }
        }

}

TEST(TestGenotypeReader, CanGetFrequency){
    double result = 0;
    Genotype_Reader reader = Genotype_Reader(nullptr);
    reader.sample_to_index = new int[4];
    for(int i = 0; i < 4; i++)
        reader.sample_to_index[i] = i+1;
    reader.num_samples = 4;
    reader.buffer = new char[100];
    reader.minor_allele_cutoff = 1;

    strcpy(reader.buffer, "0 0 1 0 0");
    // fail since cutoff too low
    ASSERT_TRUE(!reader.get_frequency(result));
    reader.minor_allele_cutoff = 0;
    ASSERT_TRUE(reader.get_frequency(result));
    ASSERT_EQ(0.125, result);

    reader.minor_allele_cutoff = 1;
    strcpy(reader.buffer, "0 0 2 0 0");
    ASSERT_TRUE(reader.get_frequency(result));
    ASSERT_EQ(0.25, result);

    reader.minor_allele_cutoff = 0;
    strcpy(reader.buffer, "0 0 2 9 9");
    ASSERT_TRUE(reader.get_frequency(result));
    ASSERT_EQ(0.5, result);

    strcpy(reader.buffer, "0 9 2 9 9");
    ASSERT_TRUE(!reader.get_frequency(result));
    ASSERT_EQ(1, result);

    strcpy(reader.buffer, "0 9 9 9 9");
    ASSERT_TRUE(!reader.get_frequency(result));
    ASSERT_EQ(0, result);

    strcpy(reader.buffer, "0 0 0 0 0");
    ASSERT_TRUE(!reader.get_frequency(result));
    ASSERT_EQ(0, result);

}

TEST(TestGenotypeReader, CanGetModernError){
    Genotype_Reader reader = Genotype_Reader(nullptr, nullptr, 0.01, 0.1, 2);

    ASSERT_DOUBLE_EQ(0.1, reader.get_modern_error(0.4));
    ASSERT_DOUBLE_EQ(0.1, reader.get_modern_error(0.6));
    ASSERT_DOUBLE_EQ(0.1, reader.get_modern_error(0.1));
    ASSERT_DOUBLE_EQ(0.1, reader.get_modern_error(0.9));
    ASSERT_DOUBLE_EQ(0.02, reader.get_modern_error(0.01));
    ASSERT_TRUE(abs((reader.get_modern_error(0.99) - 0.02)/0.02) < 0.001);
}

TEST(TestGenotypeReader, CanProcessLineBuffer){
    Genotype_Reader reader = Genotype_Reader(nullptr);
    reader.buffer = new char[100];
    reader.num_samples = 4;
    reader.sample_to_index = new int[4];
    for(int i = 0; i < 4; i++)
        reader.sample_to_index[i] = 4-i;
    reader.archaic_index = 0;
    reader.lod_scores = new double[4];
    reader.recover_type = new unsigned char[4];

    strcpy(reader.buffer, "0 0 1 2 9");
    reader.process_line_buffer(true);
    double result[] = {0, -1.617, 0.00432, 0.300};
    for(int i = 0; i < 4; i++)
        ASSERT_TRUE(result[i] == reader.lod_scores[i] ||
                abs((reader.lod_scores[i] - result[i])/result[i]) < 0.001);

    reader.process_line_buffer(false);
    result[0] = 0;
    result[1] = -1.617;
    result[2] = 0;
    result[3] = 0;
    for(int i = 0; i < 4; i++)
        ASSERT_TRUE(result[i] == reader.lod_scores[i] ||
                abs((reader.lod_scores[i] - result[i])/result[i]) < 0.001);

    strcpy(reader.buffer, "9 0 1 2 9");
    reader.process_line_buffer(true);
    result[1] = 0;
    for(int i = 0; i < 4; i++)
        ASSERT_TRUE(result[i] == reader.lod_scores[i] ||
                abs((reader.lod_scores[i] - result[i])/result[i]) < 0.001);

    reader.process_line_buffer(false);
    for(int i = 0; i < 4; i++)
        ASSERT_TRUE(result[i] == reader.lod_scores[i] ||
                abs((reader.lod_scores[i] - result[i])/result[i]) < 0.001);

    strcpy(reader.buffer, "0 0 1 0 9");
    reader.process_line_buffer(true);
    for(int i = 0; i < 4; i++)
        ASSERT_TRUE(result[i] == reader.lod_scores[i] ||
                abs((reader.lod_scores[i] - result[i])/result[i]) < 0.001);

    strcpy(reader.buffer, "1 1 1 2 9");
    reader.process_line_buffer(true);
    result[1] = 0.16758;
    result[2] = 0.34366;
    result[3] = 0.34366;
    for(int i = 0; i < 4; i++)
        ASSERT_TRUE(result[i] == reader.lod_scores[i] ||
                abs((reader.lod_scores[i] - result[i])/result[i]) < 0.001);
        //std::cout << reader.lod_scores[i] << "\n";
}

TEST_F(SampleGenotype, CanTestInMask){
    Genotype_Reader reader = Genotype_Reader(
            nullptr, fopen("test_mask.txt", "rt"));

    reader.chromosome = 1;
    reader.position = 90;
    ASSERT_TRUE(!reader.in_mask());

    // check same position
    ASSERT_TRUE(!reader.in_mask());
    ASSERT_TRUE(!reader.in_mask());
    ASSERT_TRUE(!reader.in_mask());

    reader.position = 100;
    ASSERT_TRUE(!reader.in_mask());

    reader.position = 101;
    ASSERT_TRUE(reader.in_mask());

    reader.position = 102;
    ASSERT_TRUE(reader.in_mask());

    reader.position = 120;
    ASSERT_TRUE(reader.in_mask());

    // skip a range
    reader.position = 161;
    ASSERT_TRUE(reader.in_mask());

    reader.position = 261;
    ASSERT_TRUE(reader.in_mask());
    
    //skip chromosomes
    reader.chromosome = 2;
    reader.position = 130;
    ASSERT_TRUE(!reader.in_mask());

    reader.chromosome = 3;
    reader.position = 130;
    ASSERT_TRUE(!reader.in_mask());

    reader.chromosome = 4;
    reader.position = 131;
    ASSERT_TRUE(reader.in_mask());

    reader.chromosome = 5;
    reader.position = 131;
    ASSERT_TRUE(!reader.in_mask());

    reader.chromosome = 6;
    reader.position = 131;
    ASSERT_TRUE(!reader.in_mask());

    reader.chromosome = 8;
    reader.position = 131;
    ASSERT_TRUE(reader.in_mask());

    reader.chromosome = 9;
    reader.position = 131;
    ASSERT_TRUE(!reader.in_mask());
    ASSERT_TRUE(!reader.in_mask());
    ASSERT_TRUE(!reader.in_mask());
    ASSERT_TRUE(!reader.in_mask());

    reader.chromosome = 8;
    // still retain the last read values (not actual use case)
    ASSERT_TRUE(reader.in_mask());

    // setting mask to null should short the region check
    Genotype_Reader reader2(fopen("test_genotype.txt", "rt"));
    ASSERT_TRUE(!reader2.in_mask());
}

TEST_F(SampleGenotype, CanUpdateDefaults)
{
    Genotype_Reader reader = Genotype_Reader(fopen("test_genotype.txt", "rt"));
    reader.initialize();
    ASSERT_EQ(4, reader.num_samples);

    //"1\t2\tA\tT\t1\t0\t0\t0\t0\n"  fails allele check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(2, reader.position);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[0]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[1]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[2]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[3]);

    // "1\t3\tA\tT\t2\t0\t0\t0\t0\n" fails check, 2/0 override
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(3, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - -1.99568)/-1.99568) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - -1.99568)/-1.99568) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[2] - -1.99568)/-1.99568) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[3] - -1.99568)/-1.99568) < 0.001);

    // "1\t4\tA\tT\t1\t0\t1\t1\t1\n"  passes check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(4, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - 0.195605)/0.195605) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - 0.320544)/0.320544) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[2] - 0.320544)/0.320544) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[3] - 0.320544)/0.320544) < 0.001);

    // "1\t104\tA\tT\t1\t0\t1\t1\t1\n" same as above
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(104, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - 0.195605)/0.195605) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - 0.320544)/0.320544) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[2] - 0.320544)/0.320544) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[3] - 0.320544)/0.320544) < 0.001);

    // "1\t105\tA\tT\t0\t2\t1\t1\t1\n"  passes check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(105, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - -1.713827)/-1.713827) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - 0.127177)/0.127177) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[2] - 0.127177)/0.127177) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[3] - 0.127177)/0.127177) < 0.001);

    // "2\t125\tA\tT\t0\t2\t2\t1\t1\n" change chromosome
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(2, reader.chromosome);
    ASSERT_EQ(125, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - -1.79301)/-1.79301) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - -1.79301)/-1.79301) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[2] - 0.301874)/0.301874) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[3] - 0.301874)/0.301874) < 0.001);

    // "3\t126\tA\tT\t0\t2\t2\t2\t2\n" change chrom, 0/2 override check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(3, reader.chromosome);
    ASSERT_EQ(126, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - -1.99568)/-1.99568) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - -1.99568)/-1.99568) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[2] - -1.99568)/-1.99568) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[3] - -1.99568)/-1.99568) < 0.001);

    // eof
    ASSERT_TRUE(!reader.update());
    ASSERT_TRUE(!reader.update());
    ASSERT_TRUE(!reader.update());
}

TEST_F(SampleGenotype, CanUpdateSamples)
{
    Genotype_Reader reader = Genotype_Reader(fopen("test_genotype.txt", "rt"));
    std::ofstream samples;
    samples.open("test_samples.txt");
    samples << "m4\nm2\n";
    samples.close();
    reader.initialize(fopen("test_samples.txt", "rt"), "m1");
    ASSERT_EQ(2, reader.num_samples);

    //"1\t2\tA\tT\t1\t0\t0\t0\t0\n"  fails allele check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(2, reader.position);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[0]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[1]);

    // "1\t3\tA\tT\t2\t0\t0\t0\t0\n" fails check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(3, reader.position);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[0]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[1]);

    // "1\t4\tA\tT\t1\t0\t1\t1\t1\n"  passes check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(4, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - 0.00432094)/0.00432094) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - 0.00432094)/0.00432094) < 0.001);

    // "1\t104\tA\tT\t1\t0\t1\t1\t1\n" same as above
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(104, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - 0.00432094)/0.00432094) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - 0.00432094)/0.00432094) < 0.001);

    // "1\t105\tA\tT\t0\t2\t1\t1\t1\n"  passes check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(105, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - 0.00432094)/0.00432094) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - 0.00432094)/0.00432094) < 0.001);

    // "2\t125\tA\tT\t0\t2\t2\t1\t1\n" change chromosome
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(2, reader.chromosome);
    ASSERT_EQ(125, reader.position);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[0]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[1]);

    // "3\t126\tA\tT\t0\t2\t2\t2\t2\n" change chrom, 0/2 override check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(3, reader.chromosome);
    ASSERT_EQ(126, reader.position);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[0]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[1]);

    // eof
    ASSERT_TRUE(!reader.update());
    ASSERT_TRUE(!reader.update());
    ASSERT_TRUE(!reader.update());
}

TEST_F(SampleGenotype, CanUpdateMask)
{
    Genotype_Reader reader = Genotype_Reader(
            fopen("test_genotype.txt", "rt"),
            fopen("test_mask.txt", "rt"));
    reader.initialize();
    ASSERT_EQ(4, reader.num_samples);

    //"1\t2\tA\tT\t1\t0\t0\t0\t0\n"  fails allele check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(2, reader.position);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[0]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[1]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[2]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[3]);

    // "1\t3\tA\tT\t2\t0\t0\t0\t0\n" fails check, 2/0 override
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(3, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - -1.99568)/-1.99568) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - -1.99568)/-1.99568) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[2] - -1.99568)/-1.99568) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[3] - -1.99568)/-1.99568) < 0.001);

    // "1\t4\tA\tT\t1\t0\t1\t1\t1\n"  passes check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(4, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - 0.195605)/0.195605) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - 0.320544)/0.320544) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[2] - 0.320544)/0.320544) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[3] - 0.320544)/0.320544) < 0.001);

    // "1\t104\tA\tT\t1\t0\t1\t1\t1\n" in mask
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(104, reader.position);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[0]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[1]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[2]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[3]);

    // "1\t105\tA\tT\t0\t2\t1\t1\t1\n"  in mask, 0, 2 override
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(105, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - -1.713827)/-1.713827) < 0.001);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[1]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[2]);
    ASSERT_DOUBLE_EQ(0, reader.lod_scores[3]);

    // "2\t125\tA\tT\t0\t2\t2\t1\t1\n" change chromosome
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(2, reader.chromosome);
    ASSERT_EQ(125, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - -1.79301)/-1.79301) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - -1.79301)/-1.79301) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[2] - 0.301874)/0.301874) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[3] - 0.301874)/0.301874) < 0.001);

    // "3\t126\tA\tT\t0\t2\t2\t2\t2\n" change chrom, 0/2 override check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(3, reader.chromosome);
    ASSERT_EQ(126, reader.position);
    ASSERT_TRUE(abs((reader.lod_scores[0] - -1.99568)/-1.99568) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[1] - -1.99568)/-1.99568) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[2] - -1.99568)/-1.99568) < 0.001);
    ASSERT_TRUE(abs((reader.lod_scores[3] - -1.99568)/-1.99568) < 0.001);

    // eof
    ASSERT_TRUE(!reader.update());
    ASSERT_TRUE(!reader.update());
    ASSERT_TRUE(!reader.update());
}
