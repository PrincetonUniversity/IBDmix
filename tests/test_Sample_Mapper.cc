#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <stdio.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "IBDmix/Sample_Mapper.h"

using ::testing::ElementsAre;

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
  ASSERT_EQ(mapper.getArchaicIndex(), 0);
  ASSERT_THAT(mapper.getSamples(), ElementsAre("m1", "m2", "m3", "m4"));

  genotype.clear();
  genotype.str(genotype_contents);
  result = mapper.initialize(genotype, sample_dummy, "m3");
  ASSERT_EQ(result, 4);
  ASSERT_EQ(mapper.getArchaicIndex(), 3);
  ASSERT_THAT(mapper.getSamples(), ElementsAre("n1", "m1", "m2", "m4"));

  genotype.str(genotype_contents);
  genotype.clear();
  std::istringstream samples("m1\nm2\nm3");
  result = mapper.initialize(genotype, samples, "n1");
  ASSERT_EQ(result, 3);
  ASSERT_EQ(mapper.getArchaicIndex(), 0);
  ASSERT_THAT(mapper.getSamples(), ElementsAre("m1", "m2", "m3"));

  genotype.str(genotype_contents);
  genotype.clear();
  samples.str("m1\nm1\nm2");
  samples.clear();
  result = mapper.initialize(genotype, samples, "m3");
  ASSERT_EQ(result, 3);
  ASSERT_EQ(mapper.getArchaicIndex(), 3);
  ASSERT_THAT(mapper.getSamples(), ElementsAre("m1", "m1", "m2"));

  // setting archaic equal to a sample is legal, but probably not useful
  genotype.str(genotype_contents);
  genotype.clear();
  samples.str("m1\nm1\nm2");
  samples.clear();
  result = mapper.initialize(genotype, samples, "m1");
  ASSERT_EQ(result, 3);
  ASSERT_EQ(mapper.getArchaicIndex(), 1);
  ASSERT_THAT(mapper.getSamples(), ElementsAre("m1", "m1", "m2"));
}
