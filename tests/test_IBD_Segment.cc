#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <memory>

#include "IBDmix/IBD_Segment.h"
#include "IBDmix/Segment_Recorders.h"
#include "IBDmix/IBD_Stack.h"
#include "IBDmix/Genotype_Reader.h"

unsigned char none = '\0';

TEST(IBDSegment, CanConstruct) {
    IBD_Pool pool(5);
    IBD_Segment s1("test", 0, &pool);
    ASSERT_EQ(s1.size(), 0);
}

TEST(IBDSegment, CanAddBasicLOD) {
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg("test", 0, &pool);
    ASSERT_EQ(seg.size(), 0);

    seg.add_lod(1, 1, -1, none, output);
    ASSERT_EQ(seg.size(), 0);

    seg.add_lod(1, 2, 0.5, none, output);
    ASSERT_EQ(seg.size(), 1);

    // add some decreasing values to keep end at start
    seg.add_lod(1, 3, -0.1, none, output);
    seg.add_lod(1, 4, -0.1, none, output);
    seg.add_lod(1, 5, -0.1, none, output);
    ASSERT_EQ(seg.size(), 4);

    // add an increasing, back to 0.4
    seg.add_lod(1, 6, 0.2, none, output);
    ASSERT_EQ(seg.size(), 5);

    // add new maxes (equal, then more)
    seg.add_lod(1, 7, 0.1, none, output);
    ASSERT_EQ(seg.size(), 2);

    seg.add_lod(1, 9, -0.1, none, output);
    seg.add_lod(1, 10, -0.1, none, output);
    seg.add_lod(1, 11, 0.3, none, output);
    ASSERT_EQ(seg.size(), 2);
}

TEST(IBDSegment, CanAddLODOutput) {
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg("test", 0, &pool);

    // add some positions
    seg.add_lod(2, 1, 0.1, none, output);
    seg.add_lod(2, 2, 0.1, none, output);
    seg.add_lod(2, 3, 0.1, none, output);
    seg.add_lod(2, 4, 0.1, none, output);
    ASSERT_EQ(seg.size(), 2);

    // make cumsum 0
    seg.add_lod(2, 50, -0.4, none, output);
    ASSERT_EQ(seg.size(), 3);
    // and less than 0
    seg.add_lod(2, 60, -0.1, none, output);
    ASSERT_EQ(output.str(), "test\t2\t1\t50\t0.4\n");
    ASSERT_EQ(seg.size(), 0);
    ASSERT_EQ(pool.size(), 5);

    // have start = end
    output.str("");
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -3, none, output);
    ASSERT_EQ(output.str(), "test\t2\t1\t2\t2\n");
    ASSERT_EQ(seg.size(), 0);
    ASSERT_EQ(pool.size(), 5);

    // generate multiple outputs at once
    output.str("");
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -1, none, output);
    seg.add_lod(2, 3, 0.5, none, output);
    seg.add_lod(2, 4, -1, none, output);
    seg.add_lod(2, 5, 0.7, none, output);
    seg.add_lod(2, 6, -2, none, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\ntest\t2\t3\t4\t0.5\ntest\t2\t5\t6\t0.7\n");
    ASSERT_EQ(seg.size(), 0);
    ASSERT_EQ(pool.size(), 10);

    // split last into two positions
    output.str("");
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -1, none, output);
    seg.add_lod(2, 3, 0.5, none, output);
    seg.add_lod(2, 4, -1, none, output);
    seg.add_lod(2, 5, 0.3, none, output);
    seg.add_lod(2, 6, 0.3, none, output);
    seg.add_lod(2, 7, -2, none, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\ntest\t2\t3\t4\t0.5\ntest\t2\t5\t7\t0.6\n");
    ASSERT_EQ(seg.size(), 0);
    ASSERT_EQ(pool.size(), 10);

    // trigger reversal twice
    output.str("");
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -1, none, output);
    seg.add_lod(2, 3, 0.5, none, output);
    seg.add_lod(2, 4, -1, none, output);
    seg.add_lod(2, 5, 0.5, none, output);
    seg.add_lod(2, 6, -1, none, output);
    seg.add_lod(2, 7, 1.9, none, output);
    seg.add_lod(2, 8, -1, none, output);
    seg.add_lod(2, 9, 0.5, none, output);
    seg.add_lod(2, 10, -3, none, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\n"
            "test\t2\t3\t4\t0.5\n"
            "test\t2\t5\t6\t0.5\n"
            "test\t2\t7\t8\t1.9\n"
            "test\t2\t9\t10\t0.5\n");
    ASSERT_EQ(seg.size(), 0);
    ASSERT_EQ(pool.size(), 10);
}

TEST(IBDSegment, CanPurge) {
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg("test", 0, &pool);
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -1, none, output);
    seg.add_lod(2, 3, 0.5, none, output);
    seg.add_lod(2, 4, -1, none, output);
    seg.add_lod(2, 5, 0.5, none, output);
    seg.add_lod(2, 6, -1, none, output);
    seg.add_lod(2, 7, 1.9, none, output);
    seg.add_lod(2, 8, -1, none, output);
    seg.add_lod(2, 9, 0.5, none, output);
    ASSERT_EQ('\0', output.str().c_str()[0]);
    seg.purge(output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\n"
            "test\t2\t3\t4\t0.5\n"
            "test\t2\t5\t6\t0.5\n"
            "test\t2\t7\t8\t1.9\n"
            "test\t2\t9\t9\t0.5\n");
}

TEST(IBDSegment, CanPurgeInclusive) {
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg("test", 0, &pool, false);
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -1, none, output);
    seg.add_lod(2, 3, 0.5, none, output);
    seg.add_lod(2, 4, -1, none, output);
    seg.add_lod(2, 5, 0.5, none, output);
    seg.add_lod(2, 6, -1, none, output);
    seg.add_lod(2, 7, 1.9, none, output);
    seg.add_lod(2, 8, -1, none, output);
    seg.add_lod(2, 9, 0.5, none, output);
    ASSERT_EQ('\0', output.str().c_str()[0]);
    seg.purge(output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t1\t2\n"
            "test\t2\t3\t3\t0.5\n"
            "test\t2\t5\t5\t0.5\n"
            "test\t2\t7\t7\t1.9\n"
            "test\t2\t9\t9\t0.5\n");
}

TEST(IBDSegment, CanPrint) {
    IBD_Pool pool(5);
    std::ostringstream output, print_out;
    IBD_Segment seg("test", 0, &pool);
    print_out << seg;
    ASSERT_STREQ(print_out.str().c_str(),
                    "--- test ---\n");
    print_out.str("");
    print_out.clear();

    seg.add_lod(2, 1, 2, none, output);
    print_out << seg;
    ASSERT_STREQ(print_out.str().c_str(),
                    "--- test ---\n"
                    "1\t2\t2 <- top <- start <- end\n");
    print_out.clear();
    print_out.str("");

    seg.add_lod(2, 2, -1, none, output);
    print_out << seg;
    ASSERT_STREQ(print_out.str().c_str(),
                    "--- test ---\n"
                    "2\t-1\t1 <- top\n"
                    "1\t2\t2 <- start <- end\n");
    print_out.clear();
    print_out.str("");

    seg.add_lod(2, 3, 0.5, none, output);
    seg.add_lod(2, 4, -1, none, output);
    seg.add_lod(2, 5, 0.5, none, output);
    seg.add_lod(2, 6, -1, none, output);
    print_out << seg;
    ASSERT_STREQ(print_out.str().c_str(),
                    "--- test ---\n"
                    "6\t-1\t0 <- top\n"
                    "5\t0.5\t1\n"
                    "4\t-1\t0.5\n"
                    "3\t0.5\t1.5\n"
                    "2\t-1\t1\n"
                    "1\t2\t2 <- start <- end\n");
    print_out.clear();
    print_out.str("");

    seg.add_lod(2, 7, 2.1, none, output);
    seg.add_lod(2, 8, -1, none, output);
    // collapse and move end
    print_out << seg;
    ASSERT_STREQ(print_out.str().c_str(),
                    "--- test ---\n"
                    "8\t-1\t1.1 <- top\n"
                    "7\t2.1\t2.1 <- end\n"
                    "1\t2\t2 <- start\n");
}

TEST(IBDSegment, CanRecordStats) {
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg("test", 0, &pool, true);
    seg.add_recorder(std::make_shared<CountRecorder>());
    seg.add_lod(2, 1, 1, IN_MASK, output);
    seg.add_lod(2, 2, 1, IN_MASK | MAF_LOW, output);
    seg.add_lod(2, 3, 1, IN_MASK | MAF_HIGH, output);
    seg.add_lod(2, 4, 1, MAF_LOW, output);
    seg.add_lod(2, 5, 1, MAF_HIGH, output);
    seg.add_lod(2, 6, 1, RECOVER_2_0, output);
    seg.add_lod(2, 7, 1, RECOVER_0_2, output);
    seg.add_lod(2, 8, 1, RECOVER_0_2 | IN_MASK, output);
    seg.add_lod(2, 9, 1, RECOVER_0_2 | MAF_LOW, output);
    ASSERT_EQ('\0', output.str().c_str()[0]);
    seg.purge(output);
    ASSERT_EQ(output.str(),
            // starting at last 9, sites positve_lod both
            // in_mask maf_low maf_high 2_0 0_2
            "test\t2\t1\t9\t9\t9\t9\t2\t2\t2\t1\t1\t3\n");
    output.str("");
    output.clear();

    // add some positions
    seg.add_lod(2, 1, 0.1, IN_MASK, output);
    seg.add_lod(2, 2, 0.1, IN_MASK, output);
    seg.add_lod(2, 3, 0.1, IN_MASK, output);
    seg.add_lod(2, 4, 0.1, IN_MASK, output);
    ASSERT_EQ(seg.size(), 2);

    // make cumsum 0
    seg.add_lod(2, 50, -0.4, MAF_LOW, output);
    ASSERT_EQ(seg.size(), 3);
    // and less than 0
    seg.add_lod(2, 60, -0.1, MAF_HIGH, output);
    // last site (50) doesn't apply as it's excluded
    ASSERT_EQ(output.str(),
            "test\t2\t1\t50\t0.4\t4\t4\t0\t4\t0\t0\t0\t0\n");
    ASSERT_EQ(seg.size(), 0);

    // have start = end
    output.str("");
    seg.add_lod(2, 1, 2, MAF_LOW, output);
    seg.add_lod(2, 2, -3, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\t1\t1\t0\t0\t1\t0\t0\t0\n");
    ASSERT_EQ(seg.size(), 0);

    // generate multiple outputs at once
    output.str("");
    seg.add_lod(2, 1, 2, IN_MASK, output);
    seg.add_lod(2, 2, -1, IN_MASK, output);
    seg.add_lod(2, 3, 0.5, MAF_LOW, output);
    seg.add_lod(2, 4, -1, IN_MASK, output);
    seg.add_lod(2, 5, 0.7, MAF_HIGH, output);
    seg.add_lod(2, 6, -2, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\t1\t1\t0\t1\t0\t0\t0\t0\n"
            "test\t2\t3\t4\t0.5\t1\t1\t0\t0\t1\t0\t0\t0\n"
            "test\t2\t5\t6\t0.7\t1\t1\t0\t0\t0\t1\t0\t0\n");
    ASSERT_EQ(seg.size(), 0);

    // split into more positions
    output.str("");
    seg.add_lod(2, 1, 2, IN_MASK, output);
    seg.add_lod(2, 2, -1, IN_MASK, output);
    seg.add_lod(2, 3, 0.1, MAF_LOW, output);
    seg.add_lod(2, 4, 0.1, MAF_LOW, output);
    seg.add_lod(2, 5, 0.1, MAF_LOW, output);
    seg.add_lod(2, 6, 0.1, MAF_LOW, output);
    seg.add_lod(2, 7, 0.1, MAF_LOW, output);
    seg.add_lod(2, 8, -1, IN_MASK, output);
    seg.add_lod(2, 9, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 10, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 11, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 12, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 13, -2, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\t1\t1\t0\t1\t0\t0\t0\t0\n"
            "test\t2\t3\t8\t0.5\t5\t5\t0\t0\t5\t0\t0\t0\n"
            "test\t2\t9\t13\t1.2\t4\t4\t0\t0\t0\t4\t0\t0\n");
    ASSERT_EQ(seg.size(), 0);

    // trigger reversal twice
    output.str("");
    seg.add_lod(2, 1, 2, IN_MASK, output);
    seg.add_lod(2, 2, -1, IN_MASK, output);
    seg.add_lod(2, 3, 0.5, MAF_LOW, output);
    seg.add_lod(2, 4, -1, IN_MASK, output);
    seg.add_lod(2, 5, 0.5, MAF_HIGH, output);
    seg.add_lod(2, 6, -1, IN_MASK, output);
    seg.add_lod(2, 7, 1.9, RECOVER_0_2, output);
    seg.add_lod(2, 8, -1, IN_MASK, output);
    seg.add_lod(2, 9, 0.5, RECOVER_2_0, output);
    seg.add_lod(2, 10, -3, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\t1\t1\t0\t1\t0\t0\t0\t0\n"
            "test\t2\t3\t4\t0.5\t1\t1\t0\t0\t1\t0\t0\t0\n"
            "test\t2\t5\t6\t0.5\t1\t1\t0\t0\t0\t1\t0\t0\n"
            "test\t2\t7\t8\t1.9\t1\t1\t0\t0\t0\t0\t0\t1\n"
            "test\t2\t9\t10\t0.5\t1\t1\t0\t0\t0\t0\t1\t0\n");
    ASSERT_EQ(seg.size(), 0);

    // one output with a late max
    output.str("");
    seg.add_lod(2, 1, 5, IN_MASK, output);
    seg.add_lod(2, 2, -1, MAF_LOW, output);
    seg.add_lod(2, 3, -1, MAF_LOW, output);
    seg.add_lod(2, 4, -1, MAF_LOW, output);
    seg.add_lod(2, 5, -0.5, MAF_LOW, output);
    seg.add_lod(2, 6, -0.5, MAF_LOW, output);
    seg.add_lod(2, 7, 20, IN_MASK, output);
    seg.add_lod(2, 8, 5, MAF_HIGH, output);
    seg.add_lod(2, 9, -30, RECOVER_0_2, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t9\t26\t8\t3\t0\t2\t5\t1\t0\t0\n");
    ASSERT_EQ(seg.size(), 0);
}

TEST(IBDSegment, CanRecordStatsInclusive) {
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg("test", 0, &pool, false);
    seg.add_recorder(std::make_shared<CountRecorder>());
    seg.add_lod(2, 1, 1, IN_MASK, output);
    seg.add_lod(2, 2, 1, IN_MASK | MAF_LOW, output);
    seg.add_lod(2, 3, 1, IN_MASK | MAF_HIGH, output);
    seg.add_lod(2, 4, 1, MAF_LOW, output);
    seg.add_lod(2, 5, 1, MAF_HIGH, output);
    seg.add_lod(2, 6, 1, RECOVER_2_0, output);
    seg.add_lod(2, 7, 1, RECOVER_0_2, output);
    seg.add_lod(2, 8, 1, RECOVER_0_2 | IN_MASK, output);
    seg.add_lod(2, 9, 1, RECOVER_0_2 | MAF_LOW, output);
    ASSERT_EQ('\0', output.str().c_str()[0]);
    seg.purge(output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t9\t9\t9\t9\t2\t2\t2\t1\t1\t3\n");
    output.str("");
    output.clear();

    // add some positions
    seg.add_lod(2, 1, 0.1, IN_MASK, output);
    seg.add_lod(2, 2, 0.1, IN_MASK, output);
    seg.add_lod(2, 3, 0.1, IN_MASK, output);
    seg.add_lod(2, 4, 0.1, IN_MASK, output);
    ASSERT_EQ(seg.size(), 2);

    // make cumsum 0
    seg.add_lod(2, 50, -0.4, MAF_LOW, output);
    ASSERT_EQ(seg.size(), 3);
    // and less than 0
    seg.add_lod(2, 60, -0.1, MAF_HIGH, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t4\t0.4\t4\t4\t0\t4\t0\t0\t0\t0\n");
    ASSERT_EQ(seg.size(), 0);

    // have start = end
    output.str("");
    seg.add_lod(2, 1, 2, MAF_LOW, output);
    seg.add_lod(2, 2, -3, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t1\t2\t1\t1\t0\t0\t1\t0\t0\t0\n");
    ASSERT_EQ(seg.size(), 0);

    // generate multiple outputs at once
    output.str("");
    seg.add_lod(2, 1, 2, IN_MASK, output);
    seg.add_lod(2, 2, -1, IN_MASK, output);
    seg.add_lod(2, 3, 0.5, MAF_LOW, output);
    seg.add_lod(2, 4, -1, IN_MASK, output);
    seg.add_lod(2, 5, 0.7, MAF_HIGH, output);
    seg.add_lod(2, 6, -2, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t1\t2\t1\t1\t0\t1\t0\t0\t0\t0\n"
            "test\t2\t3\t3\t0.5\t1\t1\t0\t0\t1\t0\t0\t0\n"
            "test\t2\t5\t5\t0.7\t1\t1\t0\t0\t0\t1\t0\t0\n");
    ASSERT_EQ(seg.size(), 0);

    // split into more positions
    output.str("");
    seg.add_lod(2, 1, 2, IN_MASK, output);
    seg.add_lod(2, 2, -1, IN_MASK, output);
    seg.add_lod(2, 3, 0.1, MAF_LOW, output);
    seg.add_lod(2, 4, 0.1, MAF_LOW, output);
    seg.add_lod(2, 5, 0.1, MAF_LOW, output);
    seg.add_lod(2, 6, 0.1, MAF_LOW, output);
    seg.add_lod(2, 7, 0.1, MAF_LOW, output);
    seg.add_lod(2, 8, -1, IN_MASK, output);
    seg.add_lod(2, 9, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 10, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 11, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 12, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 13, -2, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t1\t2\t1\t1\t0\t1\t0\t0\t0\t0\n"
            "test\t2\t3\t7\t0.5\t5\t5\t0\t0\t5\t0\t0\t0\n"
            "test\t2\t9\t12\t1.2\t4\t4\t0\t0\t0\t4\t0\t0\n");
    ASSERT_EQ(seg.size(), 0);

    // trigger reversal twice
    output.str("");
    seg.add_lod(2, 1, 2, IN_MASK, output);
    seg.add_lod(2, 2, -1, IN_MASK, output);
    seg.add_lod(2, 3, 0.5, MAF_LOW, output);
    seg.add_lod(2, 4, -1, IN_MASK, output);
    seg.add_lod(2, 5, 0.5, MAF_HIGH, output);
    seg.add_lod(2, 6, -1, IN_MASK, output);
    seg.add_lod(2, 7, 1.9, RECOVER_0_2, output);
    seg.add_lod(2, 8, -1, IN_MASK, output);
    seg.add_lod(2, 9, 0.5, RECOVER_2_0, output);
    seg.add_lod(2, 10, -3, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t1\t2\t1\t1\t0\t1\t0\t0\t0\t0\n"
            "test\t2\t3\t3\t0.5\t1\t1\t0\t0\t1\t0\t0\t0\n"
            "test\t2\t5\t5\t0.5\t1\t1\t0\t0\t0\t1\t0\t0\n"
            "test\t2\t7\t7\t1.9\t1\t1\t0\t0\t0\t0\t0\t1\n"
            "test\t2\t9\t9\t0.5\t1\t1\t0\t0\t0\t0\t1\t0\n");
    ASSERT_EQ(seg.size(), 0);
}

TEST(IBDSegmentSites, CanAddLODOutput) {
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg("test", 0, &pool);
    seg.add_recorder(std::make_shared<SiteRecorder>());

    // add some positions
    seg.add_lod(2, 1, 0.1, none, output);
    seg.add_lod(2, 2, 0.1, none, output);
    seg.add_lod(2, 3, 0.1, none, output);
    seg.add_lod(2, 4, 0.1, none, output);
    ASSERT_EQ(seg.size(), 2);

    // make cumsum 0
    seg.add_lod(2, 50, -0.4, none, output);
    ASSERT_EQ(seg.size(), 3);
    // and less than 0
    seg.add_lod(2, 60, -0.1, none, output);
    ASSERT_EQ(output.str(), "test\t2\t1\t50\t0.4\t1,2,3,4\n");
    ASSERT_EQ(seg.size(), 0);
    ASSERT_EQ(pool.size(), 5);

    // have start = end
    output.str("");
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -3, none, output);
    ASSERT_EQ(output.str(), "test\t2\t1\t2\t2\t1\n");
    ASSERT_EQ(seg.size(), 0);
    ASSERT_EQ(pool.size(), 5);

    // generate multiple outputs at once
    output.str("");
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -1, none, output);
    seg.add_lod(2, 3, 0.5, none, output);
    seg.add_lod(2, 4, -1, none, output);
    seg.add_lod(2, 5, 0.7, none, output);
    seg.add_lod(2, 6, -2, none, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\t1\n"
            "test\t2\t3\t4\t0.5\t3\n"
            "test\t2\t5\t6\t0.7\t5\n");
    ASSERT_EQ(seg.size(), 0);
    ASSERT_EQ(pool.size(), 10);

    // split last into two positions
    output.str("");
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -1, none, output);
    seg.add_lod(2, 3, 0.5, none, output);
    seg.add_lod(2, 4, -1, none, output);
    seg.add_lod(2, 5, 0.3, none, output);
    seg.add_lod(2, 6, 0.3, none, output);
    seg.add_lod(2, 7, -2, none, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\t1\n"
            "test\t2\t3\t4\t0.5\t3\n"
            "test\t2\t5\t7\t0.6\t5,6\n");
    ASSERT_EQ(seg.size(), 0);
    ASSERT_EQ(pool.size(), 10);

    // trigger reversal twice
    output.str("");
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -1, none, output);
    seg.add_lod(2, 3, 0.5, none, output);
    seg.add_lod(2, 4, -1, none, output);
    seg.add_lod(2, 5, 0.5, none, output);
    seg.add_lod(2, 6, -1, none, output);
    seg.add_lod(2, 7, 1.9, none, output);
    seg.add_lod(2, 8, -1, none, output);
    seg.add_lod(2, 9, 0.3, none, output);
    seg.add_lod(2, 10, -0.1, none, output);
    seg.add_lod(2, 11, 0.3, none, output);
    seg.add_lod(2, 12, -3, none, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\t1\n"
            "test\t2\t3\t4\t0.5\t3\n"
            "test\t2\t5\t6\t0.5\t5\n"
            "test\t2\t7\t8\t1.9\t7\n"
            "test\t2\t9\t12\t0.5\t9,11\n");
    ASSERT_EQ(seg.size(), 0);
    ASSERT_EQ(pool.size(), 20);

    // one output with a late max
    output.str("");
    seg.add_lod(2, 1, 5, IN_MASK, output);
    seg.add_lod(2, 2, -1, MAF_LOW, output);
    seg.add_lod(2, 3, -1, MAF_LOW, output);
    seg.add_lod(2, 4, -1, MAF_LOW, output);
    seg.add_lod(2, 5, -0.5, MAF_LOW, output);
    seg.add_lod(2, 6, -0.5, MAF_LOW, output);
    seg.add_lod(2, 7, 0.1, MAF_LOW, output);
    seg.add_lod(2, 8, 0.1, MAF_LOW, output);
    seg.add_lod(2, 9, 20, IN_MASK, output);
    seg.add_lod(2, 10, 5, MAF_HIGH, output);
    seg.add_lod(2, 11, -30, RECOVER_0_2, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t11\t26.2\t1,7,8,9,10\n");
    ASSERT_EQ(seg.size(), 0);
}

TEST(IBDSegmentSites, CanPurge) {
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg("test", 0, &pool);
    seg.add_recorder(std::make_shared<SiteRecorder>());
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -1, none, output);
    seg.add_lod(2, 3, 0.5, none, output);
    seg.add_lod(2, 4, -1, none, output);
    seg.add_lod(2, 5, 0.5, none, output);
    seg.add_lod(2, 6, -1, none, output);
    seg.add_lod(2, 7, 1.9, none, output);
    seg.add_lod(2, 8, -1, none, output);
    seg.add_lod(2, 9, 0.5, none, output);
    seg.add_lod(2, 10, 0, none, output);
    seg.add_lod(2, 11, 0, none, output);
    ASSERT_EQ('\0', output.str().c_str()[0]);
    seg.purge(output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\t1\n"
            "test\t2\t3\t4\t0.5\t3\n"
            "test\t2\t5\t6\t0.5\t5\n"
            "test\t2\t7\t8\t1.9\t7\n"
            "test\t2\t9\t11\t0.5\t9\n");
    // though the last segment goes to 11, 10 and 11 have lod == 0
}

TEST(IBDSegmentSites, CanPurgeInclusive) {
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg("test", 0, &pool, false);
    seg.add_recorder(std::make_shared<SiteRecorder>());
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -1, none, output);
    seg.add_lod(2, 3, 0.5, none, output);
    seg.add_lod(2, 4, -1, none, output);
    seg.add_lod(2, 5, 0.5, none, output);
    seg.add_lod(2, 6, -1, none, output);
    seg.add_lod(2, 7, 1.9, none, output);
    seg.add_lod(2, 8, -1, none, output);
    seg.add_lod(2, 9, 0.5, none, output);
    ASSERT_EQ('\0', output.str().c_str()[0]);
    seg.purge(output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t1\t2\t1\n"
            "test\t2\t3\t3\t0.5\t3\n"
            "test\t2\t5\t5\t0.5\t5\n"
            "test\t2\t7\t7\t1.9\t7\n"
            "test\t2\t9\t9\t0.5\t9\n");
}

TEST(IBDSegmentSites, CanRecordStats) {
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg("test", 0, &pool);
    seg.add_recorder(std::make_shared<CountRecorder>());
    seg.add_recorder(std::make_shared<SiteRecorder>());
    seg.add_lod(2, 1, 1, IN_MASK, output);
    seg.add_lod(2, 2, 1, IN_MASK | MAF_LOW, output);
    seg.add_lod(2, 3, 1, IN_MASK | MAF_HIGH, output);
    seg.add_lod(2, 4, 1, MAF_LOW, output);
    seg.add_lod(2, 5, 1, MAF_HIGH, output);
    seg.add_lod(2, 6, 1, RECOVER_2_0, output);
    seg.add_lod(2, 7, 1, RECOVER_0_2, output);
    seg.add_lod(2, 8, 1, RECOVER_0_2 | IN_MASK, output);
    seg.add_lod(2, 9, 1, RECOVER_0_2 | MAF_LOW, output);
    ASSERT_EQ('\0', output.str().c_str()[0]);
    seg.purge(output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t9\t9\t9\t9\t2\t2\t2\t1\t1\t3\t1,2,3,4,5,6,7,8,9\n");
    output.str("");
    output.clear();

    // add some positions
    seg.add_lod(2, 1, 0.1, IN_MASK, output);
    seg.add_lod(2, 2, 0.1, IN_MASK, output);
    seg.add_lod(2, 3, 0.1, IN_MASK, output);
    seg.add_lod(2, 4, 0.1, IN_MASK, output);
    ASSERT_EQ(seg.size(), 2);

    // make cumsum 0
    seg.add_lod(2, 50, -0.4, MAF_LOW, output);
    ASSERT_EQ(seg.size(), 3);
    // and less than 0
    seg.add_lod(2, 60, -0.1, MAF_HIGH, output);
    // last site (50) doesn't apply as it's excluded
    ASSERT_EQ(output.str(),
            "test\t2\t1\t50\t0.4\t4\t4\t0\t4\t0\t0\t0\t0\t1,2,3,4\n");
    ASSERT_EQ(seg.size(), 0);

    // have start = end
    output.str("");
    seg.add_lod(2, 1, 2, MAF_LOW, output);
    seg.add_lod(2, 2, -3, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\t1\t1\t0\t0\t1\t0\t0\t0\t1\n");
    ASSERT_EQ(seg.size(), 0);

    // generate multiple outputs at once
    output.str("");
    seg.add_lod(2, 1, 2, IN_MASK, output);
    seg.add_lod(2, 2, -1, IN_MASK, output);
    seg.add_lod(2, 3, 0.5, MAF_LOW, output);
    seg.add_lod(2, 4, -1, IN_MASK, output);
    seg.add_lod(2, 5, 0.7, MAF_HIGH, output);
    seg.add_lod(2, 6, -2, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\t1\t1\t0\t1\t0\t0\t0\t0\t1\n"
            "test\t2\t3\t4\t0.5\t1\t1\t0\t0\t1\t0\t0\t0\t3\n"
            "test\t2\t5\t6\t0.7\t1\t1\t0\t0\t0\t1\t0\t0\t5\n");
    ASSERT_EQ(seg.size(), 0);

    // split into more positions
    output.str("");
    seg.add_lod(2, 1, 2, IN_MASK, output);
    seg.add_lod(2, 2, -1, IN_MASK, output);
    seg.add_lod(2, 3, 0.1, MAF_LOW, output);
    seg.add_lod(2, 4, 0.1, MAF_LOW, output);
    seg.add_lod(2, 5, 0.1, MAF_LOW, output);
    seg.add_lod(2, 6, 0.1, MAF_LOW, output);
    seg.add_lod(2, 7, 0.1, MAF_LOW, output);
    seg.add_lod(2, 8, -1, IN_MASK, output);
    seg.add_lod(2, 9, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 10, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 11, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 12, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 13, -2, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\t1\t1\t0\t1\t0\t0\t0\t0\t1\n"
            "test\t2\t3\t8\t0.5\t5\t5\t0\t0\t5\t0\t0\t0\t3,4,5,6,7\n"
            "test\t2\t9\t13\t1.2\t4\t4\t0\t0\t0\t4\t0\t0\t9,10,11,12\n");
    ASSERT_EQ(seg.size(), 0);

    // trigger reversal twice
    output.str("");
    seg.add_lod(2, 1, 2, IN_MASK, output);
    seg.add_lod(2, 2, -1, IN_MASK, output);
    seg.add_lod(2, 3, 0.5, MAF_LOW, output);
    seg.add_lod(2, 4, -1, IN_MASK, output);
    seg.add_lod(2, 5, 0.5, MAF_HIGH, output);
    seg.add_lod(2, 6, -1, IN_MASK, output);
    seg.add_lod(2, 7, 1.9, RECOVER_0_2, output);
    seg.add_lod(2, 8, -1, IN_MASK, output);
    seg.add_lod(2, 9, 0.5, RECOVER_2_0, output);
    seg.add_lod(2, 10, -3, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\t1\t1\t0\t1\t0\t0\t0\t0\t1\n"
            "test\t2\t3\t4\t0.5\t1\t1\t0\t0\t1\t0\t0\t0\t3\n"
            "test\t2\t5\t6\t0.5\t1\t1\t0\t0\t0\t1\t0\t0\t5\n"
            "test\t2\t7\t8\t1.9\t1\t1\t0\t0\t0\t0\t0\t1\t7\n"
            "test\t2\t9\t10\t0.5\t1\t1\t0\t0\t0\t0\t1\t0\t9\n");
    ASSERT_EQ(seg.size(), 0);
}

TEST(IBDSegmentSites, CanRecordStatsInclusive) {
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg("test", 0, &pool, false);
    seg.add_recorder(std::make_shared<CountRecorder>());
    seg.add_recorder(std::make_shared<SiteRecorder>());
    seg.add_lod(2, 1, 1, IN_MASK, output);
    seg.add_lod(2, 2, 1, IN_MASK | MAF_LOW, output);
    seg.add_lod(2, 3, 1, IN_MASK | MAF_HIGH, output);
    seg.add_lod(2, 4, 1, MAF_LOW, output);
    seg.add_lod(2, 5, 1, MAF_HIGH, output);
    seg.add_lod(2, 6, 1, RECOVER_2_0, output);
    seg.add_lod(2, 7, 1, RECOVER_0_2, output);
    seg.add_lod(2, 8, 1, RECOVER_0_2 | IN_MASK, output);
    seg.add_lod(2, 9, 1, RECOVER_0_2 | MAF_LOW, output);
    ASSERT_EQ('\0', output.str().c_str()[0]);
    seg.purge(output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t9\t9\t9\t9\t2\t2\t2\t1\t1\t3\t1,2,3,4,5,6,7,8,9\n");
    output.str("");
    output.clear();

    // add some positions
    seg.add_lod(2, 1, 0.1, IN_MASK, output);
    seg.add_lod(2, 2, 0.1, IN_MASK, output);
    seg.add_lod(2, 3, 0.1, IN_MASK, output);
    seg.add_lod(2, 4, 0.1, IN_MASK, output);
    ASSERT_EQ(seg.size(), 2);

    // make cumsum 0
    seg.add_lod(2, 50, -0.4, MAF_LOW, output);
    ASSERT_EQ(seg.size(), 3);
    // and less than 0
    seg.add_lod(2, 60, -0.1, MAF_HIGH, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t4\t0.4\t4\t4\t0\t4\t0\t0\t0\t0\t1,2,3,4\n");
    ASSERT_EQ(seg.size(), 0);

    // have start = end
    output.str("");
    seg.add_lod(2, 1, 2, MAF_LOW, output);
    seg.add_lod(2, 2, -3, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t1\t2\t1\t1\t0\t0\t1\t0\t0\t0\t1\n");
    ASSERT_EQ(seg.size(), 0);

    // generate multiple outputs at once
    output.str("");
    seg.add_lod(2, 1, 2, IN_MASK, output);
    seg.add_lod(2, 2, -1, IN_MASK, output);
    seg.add_lod(2, 3, 0.5, MAF_LOW, output);
    seg.add_lod(2, 4, -1, IN_MASK, output);
    seg.add_lod(2, 5, 0.7, MAF_HIGH, output);
    seg.add_lod(2, 6, -2, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t1\t2\t1\t1\t0\t1\t0\t0\t0\t0\t1\n"
            "test\t2\t3\t3\t0.5\t1\t1\t0\t0\t1\t0\t0\t0\t3\n"
            "test\t2\t5\t5\t0.7\t1\t1\t0\t0\t0\t1\t0\t0\t5\n");
    ASSERT_EQ(seg.size(), 0);

    // split into more positions
    output.str("");
    seg.add_lod(2, 1, 2, IN_MASK, output);
    seg.add_lod(2, 2, -1, IN_MASK, output);
    seg.add_lod(2, 3, 0.1, MAF_LOW, output);
    seg.add_lod(2, 4, 0.1, MAF_LOW, output);
    seg.add_lod(2, 5, 0.1, MAF_LOW, output);
    seg.add_lod(2, 6, 0.1, MAF_LOW, output);
    seg.add_lod(2, 7, 0.1, MAF_LOW, output);
    seg.add_lod(2, 8, -1, IN_MASK, output);
    seg.add_lod(2, 9, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 10, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 11, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 12, 0.3, MAF_HIGH, output);
    seg.add_lod(2, 13, -2, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t1\t2\t1\t1\t0\t1\t0\t0\t0\t0\t1\n"
            "test\t2\t3\t7\t0.5\t5\t5\t0\t0\t5\t0\t0\t0\t3,4,5,6,7\n"
            "test\t2\t9\t12\t1.2\t4\t4\t0\t0\t0\t4\t0\t0\t9,10,11,12\n");
    ASSERT_EQ(seg.size(), 0);

    // trigger reversal twice
    output.str("");
    seg.add_lod(2, 1, 2, IN_MASK, output);
    seg.add_lod(2, 2, -1, IN_MASK, output);
    seg.add_lod(2, 3, 0.5, MAF_LOW, output);
    seg.add_lod(2, 4, -1, IN_MASK, output);
    seg.add_lod(2, 5, 0.5, MAF_HIGH, output);
    seg.add_lod(2, 6, -1, IN_MASK, output);
    seg.add_lod(2, 7, 1.9, RECOVER_0_2, output);
    seg.add_lod(2, 8, -1, IN_MASK, output);
    seg.add_lod(2, 9, 0.5, RECOVER_2_0, output);
    seg.add_lod(2, 10, -3, IN_MASK, output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t1\t2\t1\t1\t0\t1\t0\t0\t0\t0\t1\n"
            "test\t2\t3\t3\t0.5\t1\t1\t0\t0\t1\t0\t0\t0\t3\n"
            "test\t2\t5\t5\t0.5\t1\t1\t0\t0\t0\t1\t0\t0\t5\n"
            "test\t2\t7\t7\t1.9\t1\t1\t0\t0\t0\t0\t0\t1\t7\n"
            "test\t2\t9\t9\t0.5\t1\t1\t0\t0\t0\t0\t1\t0\t9\n");
    ASSERT_EQ(seg.size(), 0);
}
