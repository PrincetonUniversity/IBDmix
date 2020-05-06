#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <stdio.h>

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "IBDmix/Sample_Mapper.h"

using ::testing::ElementsAre;

TEST(SampleMapper, CanFindArchaic) {
    Sample_Mapper mapper;
    std::string arch = "";
    mapper.samples = std::vector<std::string> {"n1", "s1", "s2"};

    // test no archaic, default to first sample
    mapper.find_archaic(arch);
    ASSERT_EQ(0, mapper.archaic_index);
    ASSERT_EQ(mapper.samples.size(), 3);

    // test with archaic
    arch = "s1";
    mapper.find_archaic(arch);
    ASSERT_EQ(1, mapper.archaic_index);

    arch = "s2";
    mapper.find_archaic(arch);
    ASSERT_EQ(2, mapper.archaic_index);

    arch = "s3";
    ASSERT_THROW(mapper.find_archaic(arch), std::invalid_argument);
}

TEST(SampleMapper, CanMap) {
    Sample_Mapper mapper;
    mapper.samples = std::vector<std::string> {"n1", "s1", "s2", "s3", "s4"};
    std::vector<std::string> requested {};

    mapper.archaic_index = 0;
    mapper.map(requested);
    for (int i = 0; i < 4; i++)
        ASSERT_EQ(i+1, mapper.sample_to_index[i]);
    ASSERT_EQ(mapper.sample_to_index.size(), 4);
    mapper.sample_to_index.clear();

    mapper.archaic_index = 4;
    mapper.samples = std::vector<std::string> {"s1", "s2", "s3", "s4", "n1"};
    mapper.map(requested);
    for (int i = 0; i < 4; i++)
        ASSERT_EQ(i, mapper.sample_to_index[i]);
    ASSERT_EQ(mapper.sample_to_index.size(), 4);
    mapper.sample_to_index.clear();

    mapper.archaic_index = 2;
    mapper.samples = std::vector<std::string> {"s1", "s2", "n1", "s3", "s4"};
    mapper.map(requested);
    for (int i = 0; i < 2; i++)
        ASSERT_EQ(i, mapper.sample_to_index[i]);
    for (int i = 2; i < 4; i++)
        ASSERT_EQ(i+1, mapper.sample_to_index[i]);
    ASSERT_EQ(mapper.sample_to_index.size(), 4);
    mapper.sample_to_index.clear();

    mapper.samples = std::vector<std::string> {
        "n1", "s1", "s2", "s3", "s4", "s5", "s6", "s7"};
    requested = std::vector<std::string> {
        "s6", "n1", "s3"};
    mapper.map(requested);
    int result[] = {6, 0, 3};
    for (int i = 0; i < 3; i++)
        ASSERT_EQ(result[i], mapper.sample_to_index[i]);
    ASSERT_EQ(mapper.sample_to_index.size(), 3);
    mapper.sample_to_index.clear();

    // duplicate samples
    requested = std::vector<std::string> {
        "s1", "s1", "s1"};
    mapper.samples = std::vector<std::string> {"n1", "s1"};
    mapper.map(requested);
    int result2[] = {1, 1, 1};
    for (int i = 0; i < 3; i++)
        ASSERT_EQ(result2[i], mapper.sample_to_index[i]);
    ASSERT_EQ(mapper.sample_to_index.size(), 3);
    mapper.sample_to_index.clear();

    // can't find a sample
    requested = std::vector<std::string> {
        "s1", "S1", "s1"};
    ASSERT_THROW(mapper.map(requested), std::invalid_argument);
}

TEST(SampleMapper, CanInitialize) {
    Sample_Mapper mapper;
    std::string genotype_contents(
            "chrom\tpos\tref\talt\tn1\tm1\tm2\tm3\tm4\n"
            "1\t2\tA\tT\t1\t0\t0\t0\t0\n");
    std::istringstream genotype(genotype_contents);
    std::istream sample_dummy(nullptr);
    // just genotype provided
    int result = mapper.initialize(genotype, sample_dummy);
    ASSERT_EQ(result, 4);
    ASSERT_EQ(mapper.archaic_index, 0);
    ASSERT_THAT(mapper.samples, ElementsAre("m1", "m2", "m3", "m4"));

    genotype.clear();
    genotype.str(genotype_contents);
    result = mapper.initialize(genotype, sample_dummy, "m3");
    ASSERT_EQ(result, 4);
    ASSERT_EQ(mapper.archaic_index, 3);
    ASSERT_THAT(mapper.samples, ElementsAre("n1", "m1", "m2", "m4"));

    genotype.str(genotype_contents);
    genotype.clear();
    std::istringstream samples("m1\nm2\nm3");
    result = mapper.initialize(genotype, samples, "n1");
    ASSERT_EQ(result, 3);
    ASSERT_EQ(mapper.archaic_index, 0);
    ASSERT_THAT(mapper.samples, ElementsAre("m1", "m2", "m3"));

    genotype.str(genotype_contents);
    genotype.clear();
    samples.str("m1\nm1\nm2");
    samples.clear();
    result = mapper.initialize(genotype, samples, "m3");
    ASSERT_EQ(result, 3);
    ASSERT_EQ(mapper.archaic_index, 3);
    ASSERT_THAT(mapper.samples, ElementsAre("m1", "m1", "m2"));

    // setting archaic equal to a sample is legal, but probably not useful
    genotype.str(genotype_contents);
    genotype.clear();
    samples.str("m1\nm1\nm2");
    samples.clear();
    result = mapper.initialize(genotype, samples, "m1");
    ASSERT_EQ(result, 3);
    ASSERT_EQ(mapper.archaic_index, 1);
    ASSERT_THAT(mapper.samples, ElementsAre("m1", "m1", "m2"));
}
