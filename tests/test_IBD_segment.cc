#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "IBDmix/IBD_Segment.h"
#include "IBDmix/IBD_Stack.h"

unsigned char none = '\0';

TEST(IBDSegment, CanConstruct){
    ASSERT_EQ(pool_length(), 0);
    IBD_Segment s1 = IBD_Segment("test", 0);
    ASSERT_EQ(s1.length(), 0);
}

TEST(IBDSegment, CanAddBasicLOD){
    // get at least one allocation
    reclaim_node(get_node(0));
    int pool = pool_length();
    std::ostringstream output;
    IBD_Segment seg = IBD_Segment("test", 0);
    ASSERT_EQ(seg.length(), 0);

    seg.add_lod(1, 1, -1, output, none);
    ASSERT_EQ(seg.length(), 0);

    seg.add_lod(1, 2, 0.5, output, none);
    ASSERT_EQ(seg.length(), 1);

    // add some decreasing values to keep end at start
    seg.add_lod(1, 3, -0.1, output, none);
    seg.add_lod(1, 4, -0.1, output, none);
    seg.add_lod(1, 5, -0.1, output, none);
    ASSERT_EQ(seg.length(), 4);

    // add an increasing, back to 0.4
    seg.add_lod(1, 6, 0.2, output, none);
    ASSERT_EQ(seg.length(), 5);
    
    //add new maxes (equal, then more)
    seg.add_lod(1, 7, 0.1, output, none);
    ASSERT_EQ(seg.length(), 2);

    seg.add_lod(1, 9, -0.1, output, none);
    seg.add_lod(1, 10, -0.1, output, none);
    seg.add_lod(1, 11, 0.3, output, none);
    ASSERT_EQ(seg.length(), 2);
    ASSERT_EQ(pool_length(), pool -2);
}

TEST(IBDSegment, CanAddLODOutput){
    // get at least one allocation
    reclaim_node(get_node(0));
    int pool = pool_length();
    std::ostringstream output;
    IBD_Segment seg = IBD_Segment("test", 0);

    // add some positions
    seg.add_lod(2, 1, 0.1, output, none);
    seg.add_lod(2, 2, 0.1, output, none);
    seg.add_lod(2, 3, 0.1, output, none);
    seg.add_lod(2, 4, 0.1, output, none);
    ASSERT_EQ(seg.length(), 2);

    // make cumsum 0
    seg.add_lod(2, 50, -0.4, output, none);
    ASSERT_EQ(seg.length(), 3);
    // and less than 0
    seg.add_lod(2, 60, -0.1, output, none);
    ASSERT_EQ(output.str(), "test\t2\t1\t50\t0.4\n");
    ASSERT_EQ(seg.length(), 0);
    ASSERT_EQ(pool_length(), pool);

    //have start = end
    output.str("");
    seg.add_lod(2, 1, 2, output, none);
    seg.add_lod(2, 2, -3, output, none);
    ASSERT_EQ(output.str(), "test\t2\t1\t2\t2\n");
    ASSERT_EQ(seg.length(), 0);
    ASSERT_EQ(pool_length(), pool);

    //generate multiple outputs at once
    output.str("");
    seg.add_lod(2, 1, 2, output, none);
    seg.add_lod(2, 2, -1, output, none);
    seg.add_lod(2, 3, 0.5, output, none);
    seg.add_lod(2, 4, -1, output, none);
    seg.add_lod(2, 5, 0.7, output, none);
    seg.add_lod(2, 6, -2, output, none);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\ntest\t2\t3\t4\t0.5\ntest\t2\t5\t6\t0.7\n");
    ASSERT_EQ(seg.length(), 0);
    ASSERT_EQ(pool_length(), pool);

    //split last into two positions
    output.str("");
    seg.add_lod(2, 1, 2, output, none);
    seg.add_lod(2, 2, -1, output, none);
    seg.add_lod(2, 3, 0.5, output, none);
    seg.add_lod(2, 4, -1, output, none);
    seg.add_lod(2, 5, 0.3, output, none);
    seg.add_lod(2, 6, 0.3, output, none);
    seg.add_lod(2, 7, -2, output, none);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\ntest\t2\t3\t4\t0.5\ntest\t2\t5\t7\t0.6\n");
    ASSERT_EQ(seg.length(), 0);
    ASSERT_EQ(pool_length(), pool);

    //trigger reversal twice
    output.str("");
    seg.add_lod(2, 1, 2, output, none);
    seg.add_lod(2, 2, -1, output, none);
    seg.add_lod(2, 3, 0.5, output, none);
    seg.add_lod(2, 4, -1, output, none);
    seg.add_lod(2, 5, 0.5, output, none);
    seg.add_lod(2, 6, -1, output, none);
    seg.add_lod(2, 7, 1.9, output, none);
    seg.add_lod(2, 8, -1, output, none);
    seg.add_lod(2, 9, 0.5, output, none);
    seg.add_lod(2, 10, -3, output, none);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\ntest\t2\t3\t4\t0.5\n"
            "test\t2\t5\t6\t0.5\ntest\t2\t7\t8\t1.9\n"
            "test\t2\t9\t10\t0.5\n");
    ASSERT_EQ(seg.length(), 0);
    ASSERT_EQ(pool_length(), pool);
}

TEST(IBDSegment, CanPurge){
    std::ostringstream output;
    IBD_Segment seg = IBD_Segment("test", 0);
    seg.add_lod(2, 1, 2, output, none);
    seg.add_lod(2, 2, -1, output, none);
    seg.add_lod(2, 3, 0.5, output, none);
    seg.add_lod(2, 4, -1, output, none);
    seg.add_lod(2, 5, 0.5, output, none);
    seg.add_lod(2, 6, -1, output, none);
    seg.add_lod(2, 7, 1.9, output, none);
    seg.add_lod(2, 8, -1, output, none);
    seg.add_lod(2, 9, 0.5, output, none);
    ASSERT_EQ('\0', output.str().c_str()[0]);
    seg.purge(output);
    ASSERT_EQ(output.str(),
            "test\t2\t1\t2\t2\ntest\t2\t3\t4\t0.5\n"
            "test\t2\t5\t6\t0.5\ntest\t2\t7\t8\t1.9\n"
            "test\t2\t9\t9\t0.5\n");
}

// TODO add more tests for exclusive end and more stats, display  function should use stream properly... 
