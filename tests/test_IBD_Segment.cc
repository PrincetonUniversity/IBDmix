#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "IBDmix/IBD_Segment.h"
#include "IBDmix/IBD_Stack.h"

unsigned char none = '\0';

TEST(IBDSegment, CanConstruct){
    IBD_Pool pool(5);
    IBD_Segment s1 = IBD_Segment("test", 0, &pool);
    ASSERT_EQ(s1.size(), 0);
}

TEST(IBDSegment, CanAddBasicLOD){
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg = IBD_Segment("test", 0, &pool);
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
    
    //add new maxes (equal, then more)
    seg.add_lod(1, 7, 0.1, none, output);
    ASSERT_EQ(seg.size(), 2);

    seg.add_lod(1, 9, -0.1, none, output);
    seg.add_lod(1, 10, -0.1, none, output);
    seg.add_lod(1, 11, 0.3, none, output);
    ASSERT_EQ(seg.size(), 2);
}

TEST(IBDSegment, CanAddLODOutput){
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg = IBD_Segment("test", 0, &pool);

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

    //have start = end
    output.str("");
    seg.add_lod(2, 1, 2, none, output);
    seg.add_lod(2, 2, -3, none, output);
    ASSERT_EQ(output.str(), "test\t2\t1\t2\t2\n");
    ASSERT_EQ(seg.size(), 0);
    ASSERT_EQ(pool.size(), 5);

    //generate multiple outputs at once
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

    //split last into two positions
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

    //trigger reversal twice
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
            "test\t2\t1\t2\t2\ntest\t2\t3\t4\t0.5\n"
            "test\t2\t5\t6\t0.5\ntest\t2\t7\t8\t1.9\n"
            "test\t2\t9\t10\t0.5\n");
    ASSERT_EQ(seg.size(), 0);
    ASSERT_EQ(pool.size(), 10);
}

TEST(IBDSegment, CanPurge){
    IBD_Pool pool(5);
    std::ostringstream output;
    IBD_Segment seg = IBD_Segment("test", 0, &pool);
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

// TODO add more tests for exclusive end and more stats, display  function should use stream properly... 
