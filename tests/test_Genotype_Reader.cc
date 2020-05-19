#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>

#include "IBDmix/Genotype_Reader.h"

using ::testing::ElementsAre;

class SampleGenotype : public ::testing::Test {
 protected:
  void SetUp() {
    genotype.str(
        "chrom\tpos\tref\talt\tn1\tm1\tm2\tm3\tm4\n"
        "1\t2\tA\tT\t1\t0\t0\t0\t0\n"
        "1\t3\tA\tT\t2\t0\t0\t0\t0\n"
        "1\t4\tA\tT\t1\t0\t1\t1\t1\n"
        "1\t104\tA\tT\t1\t0\t1\t1\t1\n"
        "1\t105\tA\tT\t0\t2\t1\t1\t1\n"
        "2\t125\tA\tT\t0\t2\t2\t1\t1\n"
        "3\t126\tA\tT\t0\t2\t2\t2\t2\n");
    mask.str(
        "1 100 120\n"
        "1 130 140\n"
        "1 160 161\n"
        "1 190 200\n"
        "1 260 281\n"
        "2 130 140\n"
        "4 130 140\n"
        "8 130 140\n");
  }
  std::istringstream genotype;
  std::istringstream mask;
};

TEST_F(SampleGenotype, CanInitialize) {
  // use explicit values in case defaults change
  Genotype_Reader reader(&genotype, nullptr, 0.01, 0.002, 2, 1e-200);
  std::istream sample_dummy(nullptr);
  ASSERT_EQ(4, reader.initialize(sample_dummy));
  ASSERT_THAT(reader.get_samples(), ElementsAre("m1", "m2", "m3", "m4"));
  genotype.clear();
  genotype.seekg(0);
  std::istringstream samples;

  samples.str("m2\nm4\n");
  ASSERT_EQ(2, reader.initialize(samples));
  ASSERT_THAT(reader.get_samples(), ElementsAre("m2", "m4"));

  genotype.clear();
  genotype.seekg(0);
  samples.str("m2\nm4");
  samples.clear();
  ASSERT_EQ(2, reader.initialize(samples));
  ASSERT_THAT(reader.get_samples(), ElementsAre("m2", "m4"));

  genotype.clear();
  genotype.seekg(0);
  samples.str("m2");
  samples.clear();
  ASSERT_EQ(1, reader.initialize(samples));
  ASSERT_THAT(reader.get_samples(), ElementsAre("m2"));
}

// check allele frequency and line_filtering
TEST_F(SampleGenotype, CanGetFrequency) {
  std::istringstream tempGen(
      "chrom\tpos\tref\talt\tn1\tm1\tm2\tm3\tm4\n"
      "1\t2\tA\tT\t0\t0\t1\t0\t0\n"
      "1\t3\tA\tT\t0\t0\t2\t0\t0\n"
      "1\t4\tA\tT\t0\t2\t2\t2\t2\n"
      "1\t5\tA\tT\t0\t0\t0\t0\t0\n"
      "1\t6\tA\tT\t0\t0\t2\t9\t9\n"
      "1\t7\tA\tT\t0\t9\t2\t9\t9\n"
      "1\t8\tA\tT\t0\t9\t9\t9\t9\n");
  Genotype_Reader reader(&tempGen, nullptr, 0.01, 0.002, 2, 1e-200, 0);
  std::istringstream sample("m1\nm2\nm3\nm4");
  reader.initialize(sample);

  ASSERT_TRUE(reader.update());
  ASSERT_EQ(0, reader.getLineFilter());
  ASSERT_EQ(0.125, reader.getAlleleFrequency());

  ASSERT_TRUE(reader.update());
  ASSERT_EQ(0, reader.getLineFilter());
  ASSERT_EQ(0.25, reader.getAlleleFrequency());

  ASSERT_TRUE(reader.update());
  ASSERT_EQ(MAF_HIGH, reader.getLineFilter());
  ASSERT_EQ(1, reader.getAlleleFrequency());

  ASSERT_TRUE(reader.update());
  ASSERT_EQ(MAF_LOW, reader.getLineFilter());
  ASSERT_EQ(0, reader.getAlleleFrequency());

  ASSERT_TRUE(reader.update());
  ASSERT_EQ(0, reader.getLineFilter());
  ASSERT_EQ(0.5, reader.getAlleleFrequency());

  ASSERT_TRUE(reader.update());
  ASSERT_EQ(MAF_HIGH, reader.getLineFilter());
  ASSERT_EQ(1, reader.getAlleleFrequency());

  // all 9's, maf is too low and high (no counts!)
  ASSERT_TRUE(reader.update());
  ASSERT_EQ(MAF_LOW | MAF_HIGH, reader.getLineFilter());
  ASSERT_EQ(0, reader.getAlleleFrequency());

  ASSERT_FALSE(reader.update());

  // with MAC = 1
  std::istringstream tempGen2(
      "chrom\tpos\tref\talt\tn1\tm1\tm2\tm3\tm4\n"
      "1\t2\tA\tT\t0\t0\t1\t0\t0\n"
      "1\t3\tA\tT\t0\t0\t2\t0\t0\n");
  Genotype_Reader reader2(&tempGen2, nullptr, 0.01, 0.002, 2, 1e-200, 1);
  std::istringstream sample2("m1\nm2\nm3\nm4");
  reader2.initialize(sample2);

  ASSERT_TRUE(reader2.update());
  ASSERT_EQ(MAF_LOW, reader2.getLineFilter());
  ASSERT_EQ(0.125, reader2.getAlleleFrequency());

  ASSERT_TRUE(reader2.update());
  ASSERT_EQ(0, reader2.getLineFilter());
  ASSERT_EQ(0.25, reader2.getAlleleFrequency());
}

TEST_F(SampleGenotype, CanUpdateDefaults) {
  Genotype_Reader reader(&genotype);
  std::istream sample_dummy(nullptr);
  reader.initialize(sample_dummy);
  ASSERT_EQ(4, reader.num_samples());

  // "1\t2\tA\tT\t1\t0\t0\t0\t0\n"  fails allele check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(2, reader.getPosition());
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(0));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(1));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(2));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(3));

  // "1\t3\tA\tT\t2\t0\t0\t0\t0\n" fails check, 2/0 override
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(3, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - -1.99568) / -1.99568) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - -1.99568) / -1.99568) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(2) - -1.99568) / -1.99568) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(3) - -1.99568) / -1.99568) < 0.001);

  // "1\t4\tA\tT\t1\t0\t1\t1\t1\n"  passes check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(4, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - 0.195605) / 0.195605) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - 0.320544) / 0.320544) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(2) - 0.320544) / 0.320544) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(3) - 0.320544) / 0.320544) < 0.001);

  // "1\t104\tA\tT\t1\t0\t1\t1\t1\n" same as above
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(104, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - 0.195605) / 0.195605) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - 0.320544) / 0.320544) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(2) - 0.320544) / 0.320544) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(3) - 0.320544) / 0.320544) < 0.001);

  // "1\t105\tA\tT\t0\t2\t1\t1\t1\n"  passes check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(105, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - -1.713827) / -1.713827) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - 0.127177) / 0.127177) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(2) - 0.127177) / 0.127177) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(3) - 0.127177) / 0.127177) < 0.001);

  // "2\t125\tA\tT\t0\t2\t2\t1\t1\n" change chromosome
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("2", reader.getChromosome());
  ASSERT_EQ(125, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - -1.79301) / -1.79301) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - -1.79301) / -1.79301) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(2) - 0.301874) / 0.301874) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(3) - 0.301874) / 0.301874) < 0.001);

  // "3\t126\tA\tT\t0\t2\t2\t2\t2\n" change chrom, 0/2 override check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("3", reader.getChromosome());
  ASSERT_EQ(126, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - -1.99568) / -1.99568) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - -1.99568) / -1.99568) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(2) - -1.99568) / -1.99568) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(3) - -1.99568) / -1.99568) < 0.001);

  // eof
  ASSERT_FALSE(reader.update());
  ASSERT_FALSE(reader.update());
  ASSERT_FALSE(reader.update());
}

TEST_F(SampleGenotype, CanUpdateSamples) {
  Genotype_Reader reader(&genotype);
  std::istringstream samples;
  samples.str("m4\nm2\n");
  reader.initialize(samples, "m1");
  ASSERT_EQ(2, reader.num_samples());

  // "1\t2\tA\tT\t1\t0\t0\t0\t0\n"  fails allele check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(2, reader.getPosition());
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(0));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(1));

  // "1\t3\tA\tT\t2\t0\t0\t0\t0\n" fails check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(3, reader.getPosition());
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(0));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(1));

  // "1\t4\tA\tT\t1\t0\t1\t1\t1\n"  passes check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(4, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - 0.00432094) / 0.00432094) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - 0.00432094) / 0.00432094) < 0.001);

  // "1\t104\tA\tT\t1\t0\t1\t1\t1\n" same as above
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(104, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - 0.00432094) / 0.00432094) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - 0.00432094) / 0.00432094) < 0.001);

  // "1\t105\tA\tT\t0\t2\t1\t1\t1\n"  passes check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(105, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - 0.00432094) / 0.00432094) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - 0.00432094) / 0.00432094) < 0.001);

  // "2\t125\tA\tT\t0\t2\t2\t1\t1\n" change chromosome
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("2", reader.getChromosome());
  ASSERT_EQ(125, reader.getPosition());
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(0));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(1));

  // "3\t126\tA\tT\t0\t2\t2\t2\t2\n" change chrom, 0/2 override check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("3", reader.getChromosome());
  ASSERT_EQ(126, reader.getPosition());
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(0));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(1));

  // eof
  ASSERT_FALSE(reader.update());
  ASSERT_FALSE(reader.update());
  ASSERT_FALSE(reader.update());
}

TEST_F(SampleGenotype, CanUpdateMask) {
  Genotype_Reader reader(&genotype, &mask);
  std::istream sample_dummy(nullptr);
  reader.initialize(sample_dummy);
  ASSERT_EQ(4, reader.num_samples());

  // "1\t2\tA\tT\t1\t0\t0\t0\t0\n"  fails allele check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(2, reader.getPosition());
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(0));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(1));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(2));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(3));

  // "1\t3\tA\tT\t2\t0\t0\t0\t0\n" fails check, 2/0 override
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(3, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - -1.99568) / -1.99568) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - -1.99568) / -1.99568) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(2) - -1.99568) / -1.99568) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(3) - -1.99568) / -1.99568) < 0.001);

  // "1\t4\tA\tT\t1\t0\t1\t1\t1\n"  passes check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(4, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - 0.195605) / 0.195605) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - 0.320544) / 0.320544) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(2) - 0.320544) / 0.320544) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(3) - 0.320544) / 0.320544) < 0.001);

  // "1\t104\tA\tT\t1\t0\t1\t1\t1\n" in mask
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(104, reader.getPosition());
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(0));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(1));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(2));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(3));

  // "1\t105\tA\tT\t0\t2\t1\t1\t1\n"  in mask, 0, 2 override
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(105, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - -1.713827) / -1.713827) < 0.001);
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(1));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(2));
  ASSERT_DOUBLE_EQ(0, reader.getLodScore(3));

  // "2\t125\tA\tT\t0\t2\t2\t1\t1\n" change chromosome
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("2", reader.getChromosome());
  ASSERT_EQ(125, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - -1.79301) / -1.79301) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - -1.79301) / -1.79301) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(2) - 0.301874) / 0.301874) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(3) - 0.301874) / 0.301874) < 0.001);

  // "3\t126\tA\tT\t0\t2\t2\t2\t2\n" change chrom, 0/2 override check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("3", reader.getChromosome());
  ASSERT_EQ(126, reader.getPosition());
  ASSERT_TRUE(abs((reader.getLodScore(0) - -1.99568) / -1.99568) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(1) - -1.99568) / -1.99568) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(2) - -1.99568) / -1.99568) < 0.001);
  ASSERT_TRUE(abs((reader.getLodScore(3) - -1.99568) / -1.99568) < 0.001);

  // eof
  ASSERT_FALSE(reader.update());
  ASSERT_FALSE(reader.update());
  ASSERT_FALSE(reader.update());
}

TEST_F(SampleGenotype, CanCheckLineFilter) {
  genotype.str(genotype.str() +
               "4\t136\tA\tT\t0\t2\t2\t2\t2\n"
               "4\t137\tA\tT\t2\t0\t0\t0\t0\n");
  Genotype_Reader reader(&genotype, &mask);
  std::istream sample_dummy(nullptr);
  reader.initialize(sample_dummy);
  ASSERT_EQ(4, reader.num_samples());

  // "1\t2\tA\tT\t1\t0\t0\t0\t0\n"  fails allele check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(2, reader.getPosition());
  ASSERT_EQ(MAF_LOW, reader.getLineFilter());
  ASSERT_EQ(reader.getRecoverType(0), 0);
  ASSERT_EQ(reader.getRecoverType(1), 0);
  ASSERT_EQ(reader.getRecoverType(2), 0);
  ASSERT_EQ(reader.getRecoverType(3), 0);

  // "1\t3\tA\tT\t2\t0\t0\t0\t0\n" fails check, 2/0 override
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(3, reader.getPosition());
  ASSERT_EQ(MAF_LOW, reader.getLineFilter());
  ASSERT_EQ(reader.getRecoverType(0), RECOVER_2_0);
  ASSERT_EQ(reader.getRecoverType(1), RECOVER_2_0);
  ASSERT_EQ(reader.getRecoverType(2), RECOVER_2_0);
  ASSERT_EQ(reader.getRecoverType(3), RECOVER_2_0);

  // "1\t4\tA\tT\t1\t0\t1\t1\t1\n"  passes check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(4, reader.getPosition());
  ASSERT_EQ(0, reader.getLineFilter());
  ASSERT_EQ(reader.getRecoverType(0), 0);
  ASSERT_EQ(reader.getRecoverType(1), 0);
  ASSERT_EQ(reader.getRecoverType(2), 0);
  ASSERT_EQ(reader.getRecoverType(3), 0);

  // "1\t104\tA\tT\t1\t0\t1\t1\t1\n" in mask
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(104, reader.getPosition());
  ASSERT_EQ(IN_MASK, reader.getLineFilter());
  ASSERT_EQ(reader.getRecoverType(0), 0);
  ASSERT_EQ(reader.getRecoverType(1), 0);
  ASSERT_EQ(reader.getRecoverType(2), 0);
  ASSERT_EQ(reader.getRecoverType(3), 0);

  // "1\t105\tA\tT\t0\t2\t1\t1\t1\n"  in mask, 0, 2 override
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("1", reader.getChromosome());
  ASSERT_EQ(105, reader.getPosition());
  ASSERT_EQ(IN_MASK, reader.getLineFilter());
  ASSERT_EQ(reader.getRecoverType(0), RECOVER_0_2);
  ASSERT_EQ(reader.getRecoverType(1), 0);
  ASSERT_EQ(reader.getRecoverType(2), 0);
  ASSERT_EQ(reader.getRecoverType(3), 0);

  // "2\t125\tA\tT\t0\t2\t2\t1\t1\n" change chromosome
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("2", reader.getChromosome());
  ASSERT_EQ(125, reader.getPosition());
  ASSERT_EQ(0, reader.getLineFilter());
  ASSERT_EQ(reader.getRecoverType(0), 0);
  ASSERT_EQ(reader.getRecoverType(1), 0);
  ASSERT_EQ(reader.getRecoverType(2), 0);
  ASSERT_EQ(reader.getRecoverType(3), 0);

  // "3\t126\tA\tT\t0\t2\t2\t2\t2\n" change chrom, 0/2 override check
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("3", reader.getChromosome());
  ASSERT_EQ(126, reader.getPosition());
  ASSERT_EQ(MAF_HIGH, reader.getLineFilter());
  ASSERT_EQ(reader.getRecoverType(0), RECOVER_0_2);
  ASSERT_EQ(reader.getRecoverType(1), RECOVER_0_2);
  ASSERT_EQ(reader.getRecoverType(2), RECOVER_0_2);
  ASSERT_EQ(reader.getRecoverType(3), RECOVER_0_2);

  // "4\t136\tA\tT\t0\t2\t2\t2\t2\n" in mask and maf high
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("4", reader.getChromosome());
  ASSERT_EQ(136, reader.getPosition());
  ASSERT_EQ(IN_MASK | MAF_HIGH, reader.getLineFilter());
  ASSERT_EQ(reader.getRecoverType(0), RECOVER_0_2);
  ASSERT_EQ(reader.getRecoverType(1), RECOVER_0_2);
  ASSERT_EQ(reader.getRecoverType(2), RECOVER_0_2);
  ASSERT_EQ(reader.getRecoverType(3), RECOVER_0_2);

  // "4\t137\tA\tT\t2\t0\t0\t0\t0\n" in mask and maf low
  ASSERT_TRUE(reader.update());
  ASSERT_EQ("4", reader.getChromosome());
  ASSERT_EQ(137, reader.getPosition());
  ASSERT_EQ(IN_MASK | MAF_LOW, reader.getLineFilter());
  ASSERT_EQ(reader.getRecoverType(0), RECOVER_2_0);
  ASSERT_EQ(reader.getRecoverType(1), RECOVER_2_0);
  ASSERT_EQ(reader.getRecoverType(2), RECOVER_2_0);
  ASSERT_EQ(reader.getRecoverType(3), RECOVER_2_0);

  // eof
  ASSERT_FALSE(reader.update());
}
