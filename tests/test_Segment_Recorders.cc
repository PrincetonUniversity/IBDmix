#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>
#include <sstream>

#include "IBDmix/Genotype_Reader.h"
#include "IBDmix/IBD_Stack.h"
#include "IBDmix/Segment_Recorders.h"

TEST(CountRecorder, CanWriteHeader) {
  CountRecorder counter;
  std::ostringstream oss;
  counter.writeHeader(oss);
  ASSERT_STREQ(oss.str().c_str(),
               "\tsites\tpositive_lods\tnegative_lods\tmask_and_maf\tin_mask\t"
               "maf_low\tmaf_high\trec_2_0\trec_0_2");
}

TEST(CountRecorder, CanRecord) {
  CountRecorder counter;
  std::ostringstream oss;
  IBD_Pool pool(5);
  IBD_Node *node = pool.get_node(1, 0, 0);

  counter.initializeSegment();
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t0\t0\t0\t0\t0\t0\t0\t0\t0");
  oss.str("");
  oss.clear();

  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t1\t0\t0\t0\t0\t0\t0\t0\t0");
  oss.str("");
  oss.clear();

  node->lod = 1;
  node->bitmask = IN_MASK | RECOVER_2_0;
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t2\t1\t0\t0\t1\t0\t0\t1\t0");
  oss.str("");
  oss.clear();

  node->lod = -1;
  node->bitmask = MAF_LOW | RECOVER_0_2;
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t3\t1\t1\t0\t1\t1\t0\t1\t1");
  oss.str("");
  oss.clear();

  node->bitmask = MAF_HIGH;
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t4\t1\t2\t0\t1\t1\t1\t1\t1");
  oss.str("");
  oss.clear();

  node->bitmask = MAF_HIGH | MAF_LOW;  // impossible but valid
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t5\t1\t3\t0\t1\t2\t2\t1\t1");
  oss.str("");
  oss.clear();

  node->bitmask = IN_MASK | MAF_HIGH | MAF_LOW;  // impossible but valid
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t6\t1\t4\t1\t1\t2\t2\t1\t1");
  oss.str("");
  oss.clear();

  node->bitmask = IN_MASK | MAF_HIGH;
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t7\t1\t5\t2\t1\t2\t2\t1\t1");
  oss.str("");
  oss.clear();

  node->bitmask = IN_MASK | MAF_LOW;
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t8\t1\t6\t3\t1\t2\t2\t1\t1");
  oss.str("");
  oss.clear();

  counter.initializeSegment();
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t0\t0\t0\t0\t0\t0\t0\t0\t0");
}

TEST(SiteRecorder, CanWriteHeader) {
  SiteRecorder counter;
  std::ostringstream oss;
  counter.writeHeader(oss);
  ASSERT_STREQ(oss.str().c_str(), "\tSNPs");
}

TEST(SiteRecorder, CanRecord) {
  SiteRecorder counter;
  std::ostringstream oss;
  IBD_Pool pool(5);
  IBD_Node *node = pool.get_node(1, 0, 0);

  counter.initializeSegment();
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t");
  oss.str("");
  oss.clear();

  node->lod = 1;
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t1");
  oss.str("");
  oss.clear();

  node->position = 2;
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t1,2");
  oss.str("");
  oss.clear();

  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t1,2,2");
  oss.str("");
  oss.clear();

  node->lod = -1;  // ignored
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t1,2,2");
  oss.str("");
  oss.clear();
}

TEST(LODRecorder, CanWriteHeader) {
  LODRecorder counter;
  std::ostringstream oss;
  counter.writeHeader(oss);
  ASSERT_STREQ(oss.str().c_str(), "\tLODs");
}

TEST(LODRecorder, CanRecord) {
  LODRecorder counter;
  std::ostringstream oss;
  IBD_Pool pool(5);
  IBD_Node *node = pool.get_node(1, 0, 0);

  counter.initializeSegment();
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t");
  oss.str("");
  oss.clear();

  node->lod = 1;
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t1");
  oss.str("");
  oss.clear();

  node->lod = 2;
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t1,2");
  oss.str("");
  oss.clear();

  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t1,2,2");
  oss.str("");
  oss.clear();

  node->lod = -1;  // ignored
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t1,2,2");
  oss.str("");
  oss.clear();

  node->lod = 0.123456;  // ignored
  counter.record(node);
  counter.report(oss);
  ASSERT_STREQ(oss.str().c_str(), "\t1,2,2,0.1235");
  oss.str("");
  oss.clear();

  // confirm precision reset to default 6
  oss << 0.123456789;
  ASSERT_STREQ(oss.str().c_str(), "0.123457");
}
