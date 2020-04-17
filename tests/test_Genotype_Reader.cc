#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include "IBDmix/Genotype_Reader.h"

using ::testing::ElementsAre;

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
    std::ifstream file;
    file.open("test_genotype.txt");
    // use explicit values in case defaults change
    Genotype_Reader reader(&file, nullptr, 0.01, 0.002, 2, 1e-200);
    std::istream sample_dummy(nullptr);
    ASSERT_EQ(4, reader.initialize(sample_dummy));
    ASSERT_THAT(reader.get_samples(), ElementsAre("m1", "m2", "m3", "m4"));
    file.clear();
    file.seekg(0);
    std::istringstream samples;

    samples.str("m2\nm4\n");
    ASSERT_EQ(2, reader.initialize(samples));
    ASSERT_THAT(reader.get_samples(), ElementsAre("m2", "m4"));

    file.clear();
    file.seekg(0);
    samples.str("m2\nm4");
    samples.clear();
    ASSERT_EQ(2, reader.initialize(samples));
    ASSERT_THAT(reader.get_samples(), ElementsAre("m2", "m4"));

    file.clear();
    file.seekg(0);
    samples.str("m2");
    samples.clear();
    ASSERT_EQ(1, reader.initialize(samples));
    ASSERT_THAT(reader.get_samples(), ElementsAre("m2"));
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

TEST(GenotypeReader, CanCalculateLod){
    std::ifstream file;
    Genotype_Reader reader(&file);

    double pbs[] = {0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.99, 1};
    double aerrs[] = {0.01, 0.02, 0.03};
    double merrs[] = {0.0001, 0.001, 0.002};
    double minesps[] = {1e-200, 1e-190, 1e-100};
    char gts[] = {'0', '1', '2', '9'};
    for(double pb : pbs) for(double aerr : aerrs) for(double merr : merrs)
        for(double minesp : minesps) for(char arch : gts){
            reader.archaic_error = aerr;
            reader.minesp = minesp;
            reader.update_lod_cache(arch, pb, merr, true);
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
                    ASSERT_DOUBLE_EQ(read, orig);
                }
                else{
                    ASSERT_EQ(0, read);
                }
            }
        }

}

TEST_F(SampleGenotype, CanGetFrequency){
    std::ifstream file;
    file.open("test_genotype.txt");
    Genotype_Reader reader(&file);
    std::istringstream sample("m1\nm2\nm3\nm4");
    reader.initialize(sample);

    double result = 0;
    reader.minor_allele_cutoff = 1;

    reader.buffer = "0 0 1 0 0";
    // fail since cutoff too low
    ASSERT_FALSE(reader.get_frequency(result));
    reader.minor_allele_cutoff = 0;
    ASSERT_TRUE(reader.get_frequency(result));
    ASSERT_EQ(0.125, result);

    reader.minor_allele_cutoff = 1;
    reader.buffer = "0 0 2 0 0";
    ASSERT_TRUE(reader.get_frequency(result));
    ASSERT_EQ(0.25, result);

    reader.minor_allele_cutoff = 0;
    reader.buffer = "0 0 2 9 9";
    ASSERT_TRUE(reader.get_frequency(result));
    ASSERT_EQ(0.5, result);

    reader.buffer = "0 9 2 9 9";
    ASSERT_FALSE(reader.get_frequency(result));
    ASSERT_EQ(1, result);

    reader.buffer = "0 9 9 9 9";
    ASSERT_FALSE(reader.get_frequency(result));
    ASSERT_EQ(0, result);

    reader.buffer = "0 0 0 0 0";
    ASSERT_FALSE(reader.get_frequency(result));
    ASSERT_EQ(0, result);

}

TEST(GenotypeReader, CanGetModernError){
    Genotype_Reader reader(nullptr, nullptr, 0.01, 0.1, 2);

    ASSERT_DOUBLE_EQ(0.1, reader.get_modern_error(0.4));
    ASSERT_DOUBLE_EQ(0.1, reader.get_modern_error(0.6));
    ASSERT_DOUBLE_EQ(0.1, reader.get_modern_error(0.1));
    ASSERT_DOUBLE_EQ(0.1, reader.get_modern_error(0.9));
    ASSERT_DOUBLE_EQ(0.02, reader.get_modern_error(0.01));
    ASSERT_TRUE(abs((reader.get_modern_error(0.99) - 0.02)/0.02) < 0.001);
}

TEST_F(SampleGenotype, CanProcessLineBuffer){
    std::ifstream file;
    file.open("test_genotype.txt");
    Genotype_Reader reader(&file);
    std::istringstream sample("m4\nm3\nm2\nm1");
    reader.initialize(sample);

    reader.lod_scores.reserve(4);
    reader.recover_type.reserve(4);

    reader.buffer = "0 0 1 2 9";
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

    reader.buffer = "9 0 1 2 9";
    reader.process_line_buffer(true);
    result[1] = 0;
    for(int i = 0; i < 4; i++)
        ASSERT_TRUE(result[i] == reader.lod_scores[i] ||
                abs((reader.lod_scores[i] - result[i])/result[i]) < 0.001);

    reader.process_line_buffer(false);
    for(int i = 0; i < 4; i++)
        ASSERT_TRUE(result[i] == reader.lod_scores[i] ||
                abs((reader.lod_scores[i] - result[i])/result[i]) < 0.001);

    reader.buffer = "0 0 1 0 9";
    reader.process_line_buffer(true);
    for(int i = 0; i < 4; i++)
        ASSERT_TRUE(result[i] == reader.lod_scores[i] ||
                abs((reader.lod_scores[i] - result[i])/result[i]) < 0.001);

    reader.buffer = "1 1 1 2 9";
    reader.process_line_buffer(true);
    result[1] = 0.16758;
    result[2] = 0.34366;
    result[3] = 0.34366;
    for(int i = 0; i < 4; i++)
        ASSERT_TRUE(result[i] == reader.lod_scores[i] ||
                abs((reader.lod_scores[i] - result[i])/result[i]) < 0.001);
}

TEST_F(SampleGenotype, CanUpdateDefaults) {
    std::ifstream file;
    file.open("test_genotype.txt");
    Genotype_Reader reader(&file);
    std::istream sample_dummy(nullptr);
    reader.initialize(sample_dummy);
    ASSERT_EQ(4, reader.num_samples());

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
    ASSERT_FALSE(reader.update());
    ASSERT_FALSE(reader.update());
    ASSERT_FALSE(reader.update());
}

TEST_F(SampleGenotype, CanUpdateSamples) {
    std::ifstream file;
    file.open("test_genotype.txt");
    Genotype_Reader reader(&file);
    std::ofstream samples;
    samples.open("test_samples.txt");
    samples << "m4\nm2\n";
    samples.close();
    std::ifstream sample_file;
    sample_file.open("test_samples.txt");
    reader.initialize(sample_file, "m1");
    sample_file.close();
    ASSERT_EQ(2, reader.num_samples());

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
    ASSERT_FALSE(reader.update());
    ASSERT_FALSE(reader.update());
    ASSERT_FALSE(reader.update());
}

TEST_F(SampleGenotype, CanUpdateMask) {
    std::ifstream mask;
    mask.open("test_mask.txt");
    std::ifstream file;
    file.open("test_genotype.txt");
    Genotype_Reader reader(&file, &mask);
    std::istream sample_dummy(nullptr);
    reader.initialize(sample_dummy);
    ASSERT_EQ(4, reader.num_samples());

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
    ASSERT_FALSE(reader.update());
    ASSERT_FALSE(reader.update());
    ASSERT_FALSE(reader.update());
    mask.close();
}

TEST_F(SampleGenotype, CanCheckLineFilter) {
    std::ofstream genotype;
    genotype.open("test_genotype.txt", std::ofstream::out | std::ofstream::app);
    genotype << "4\t136\tA\tT\t0\t2\t2\t2\t2\n"
        << "4\t137\tA\tT\t2\t0\t0\t0\t0\n";
    genotype.close();
    std::ifstream mask;
    mask.open("test_mask.txt");
    std::ifstream file;
    file.open("test_genotype.txt");
    Genotype_Reader reader(&file, &mask);
    std::istream sample_dummy(nullptr);
    reader.initialize(sample_dummy);
    ASSERT_EQ(4, reader.num_samples());

    //"1\t2\tA\tT\t1\t0\t0\t0\t0\n"  fails allele check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(2, reader.position);
    ASSERT_EQ(MAF_LOW, reader.line_filtering);
    ASSERT_EQ(reader.recover_type[0], 0);
    ASSERT_EQ(reader.recover_type[1], 0);
    ASSERT_EQ(reader.recover_type[2], 0);
    ASSERT_EQ(reader.recover_type[3], 0);

    // "1\t3\tA\tT\t2\t0\t0\t0\t0\n" fails check, 2/0 override
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(3, reader.position);
    ASSERT_EQ(MAF_LOW, reader.line_filtering);
    ASSERT_EQ(reader.recover_type[0], RECOVER_2_0);
    ASSERT_EQ(reader.recover_type[1], RECOVER_2_0);
    ASSERT_EQ(reader.recover_type[2], RECOVER_2_0);
    ASSERT_EQ(reader.recover_type[3], RECOVER_2_0);

    // "1\t4\tA\tT\t1\t0\t1\t1\t1\n"  passes check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(4, reader.position);
    ASSERT_EQ(0, reader.line_filtering);
    ASSERT_EQ(reader.recover_type[0], 0);
    ASSERT_EQ(reader.recover_type[1], 0);
    ASSERT_EQ(reader.recover_type[2], 0);
    ASSERT_EQ(reader.recover_type[3], 0);

    // "1\t104\tA\tT\t1\t0\t1\t1\t1\n" in mask
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(104, reader.position);
    ASSERT_EQ(IN_MASK, reader.line_filtering);
    ASSERT_EQ(reader.recover_type[0], 0);
    ASSERT_EQ(reader.recover_type[1], 0);
    ASSERT_EQ(reader.recover_type[2], 0);
    ASSERT_EQ(reader.recover_type[3], 0);

    // "1\t105\tA\tT\t0\t2\t1\t1\t1\n"  in mask, 0, 2 override
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(1, reader.chromosome);
    ASSERT_EQ(105, reader.position);
    ASSERT_EQ(IN_MASK, reader.line_filtering);
    ASSERT_EQ(reader.recover_type[0], RECOVER_0_2);
    ASSERT_EQ(reader.recover_type[1], 0);
    ASSERT_EQ(reader.recover_type[2], 0);
    ASSERT_EQ(reader.recover_type[3], 0);

    // "2\t125\tA\tT\t0\t2\t2\t1\t1\n" change chromosome
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(2, reader.chromosome);
    ASSERT_EQ(125, reader.position);
    ASSERT_EQ(0, reader.line_filtering);
    ASSERT_EQ(reader.recover_type[0], 0);
    ASSERT_EQ(reader.recover_type[1], 0);
    ASSERT_EQ(reader.recover_type[2], 0);
    ASSERT_EQ(reader.recover_type[3], 0);

    // "3\t126\tA\tT\t0\t2\t2\t2\t2\n" change chrom, 0/2 override check
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(3, reader.chromosome);
    ASSERT_EQ(126, reader.position);
    ASSERT_EQ(MAF_HIGH, reader.line_filtering);
    ASSERT_EQ(reader.recover_type[0], RECOVER_0_2);
    ASSERT_EQ(reader.recover_type[1], RECOVER_0_2);
    ASSERT_EQ(reader.recover_type[2], RECOVER_0_2);
    ASSERT_EQ(reader.recover_type[3], RECOVER_0_2);

    // "4\t136\tA\tT\t0\t2\t2\t2\t2\n" in mask and maf high
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(4, reader.chromosome);
    ASSERT_EQ(136, reader.position);
    ASSERT_EQ(IN_MASK | MAF_HIGH, reader.line_filtering);
    ASSERT_EQ(reader.recover_type[0], RECOVER_0_2);
    ASSERT_EQ(reader.recover_type[1], RECOVER_0_2);
    ASSERT_EQ(reader.recover_type[2], RECOVER_0_2);
    ASSERT_EQ(reader.recover_type[3], RECOVER_0_2);

    // "4\t137\tA\tT\t2\t0\t0\t0\t0\n" in mask and maf low
    ASSERT_TRUE(reader.update());
    ASSERT_EQ(4, reader.chromosome);
    ASSERT_EQ(137, reader.position);
    ASSERT_EQ(IN_MASK | MAF_LOW, reader.line_filtering);
    ASSERT_EQ(reader.recover_type[0], RECOVER_2_0);
    ASSERT_EQ(reader.recover_type[1], RECOVER_2_0);
    ASSERT_EQ(reader.recover_type[2], RECOVER_2_0);
    ASSERT_EQ(reader.recover_type[3], RECOVER_2_0);

    // eof
    ASSERT_FALSE(reader.update());
    mask.close();
}
