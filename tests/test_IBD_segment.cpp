#define BOOST_TEST_MODULE IBD_Segment
#include <boost/test/included/unit_test.hpp>
#include "../IBDmix/IBD_Segment.hpp"
#include "../IBDmix/IBD_Stack.hpp"

BOOST_AUTO_TEST_CASE(test_constructor){
    BOOST_REQUIRE_EQUAL(pool_length(), 0);
    IBD_Segment s1 = IBD_Segment("test", 1, 0);
    BOOST_REQUIRE_EQUAL(s1.name, "test");
    BOOST_REQUIRE_EQUAL(s1.length(), 0);
}

BOOST_AUTO_TEST_CASE(test_add_lod_basic){
    // get at least one allocation
    reclaim_node(get_node(0));
    int pool = pool_length();
    char output[100];
    IBD_Segment seg = IBD_Segment("test", 1, 0);
    BOOST_REQUIRE_EQUAL(seg.length(), 0);

    seg.add_lod(1, -1, output);
    BOOST_REQUIRE_EQUAL(seg.length(), 0);

    seg.add_lod(2, 0.5, output);
    BOOST_REQUIRE_EQUAL(seg.length(), 1);
    // add some decreasing values to keep end at start
    seg.add_lod(3, -0.1, output);
    seg.add_lod(4, -0.1, output);
    seg.add_lod(5, -0.1, output);
    BOOST_REQUIRE_EQUAL(seg.length(), 4);

    // add an increasing, back to 0.4
    seg.add_lod(6, 0.2, output);
    BOOST_REQUIRE_EQUAL(seg.length(), 5);
    
    //add new maxes (equal, then more)
    seg.add_lod(7, 0.1, output);
    BOOST_REQUIRE_EQUAL(seg.length(), 2);

    seg.add_lod(9, -0.1, output);
    seg.add_lod(10, -0.1, output);
    seg.add_lod(11, 0.3, output);
    BOOST_REQUIRE_EQUAL(seg.length(), 2);
    BOOST_REQUIRE_EQUAL(pool_length(), pool -2);
}

BOOST_AUTO_TEST_CASE(test_add_lod_output){
    // get at least one allocation
    reclaim_node(get_node(0));
    int pool = pool_length();
    char output[100];
    IBD_Segment seg = IBD_Segment("test", 2, 0);

    // add some positions
    seg.add_lod(1, 0.1, output);
    seg.add_lod(2, 0.1, output);
    seg.add_lod(3, 0.1, output);
    seg.add_lod(4, 0.1, output);
    BOOST_REQUIRE_EQUAL(seg.length(), 2);

    // make cumsum 0
    seg.add_lod(5, -0.4, output);
    BOOST_REQUIRE_EQUAL(seg.length(), 3);
    // and less than 0
    seg.add_lod(6, -0.1, output);
    BOOST_REQUIRE_EQUAL(output, "test\t2\t1\t4\t0.4\n");
    BOOST_REQUIRE_EQUAL(seg.length(), 0);
    BOOST_REQUIRE_EQUAL(pool_length(), pool);

    //have start = end
    seg.add_lod(1, 2, output);
    seg.add_lod(2, -3, output);
    BOOST_REQUIRE_EQUAL(output, "test\t2\t1\t1\t2\n");
    BOOST_REQUIRE_EQUAL(seg.length(), 0);
    BOOST_REQUIRE_EQUAL(pool_length(), pool);

    //generate multiple outputs at once
    seg.add_lod(1, 2, output);
    seg.add_lod(2, -1, output);
    seg.add_lod(3, 0.5, output);
    seg.add_lod(4, -1, output);
    seg.add_lod(5, 0.7, output);
    seg.add_lod(6, -2, output);
    BOOST_REQUIRE_EQUAL(output,
            "test\t2\t1\t1\t2\ntest\t2\t3\t3\t0.5\ntest\t2\t5\t5\t0.7\n");
    BOOST_REQUIRE_EQUAL(seg.length(), 0);
    BOOST_REQUIRE_EQUAL(pool_length(), pool);

    //split last into two positions
    seg.add_lod(1, 2, output);
    seg.add_lod(2, -1, output);
    seg.add_lod(3, 0.5, output);
    seg.add_lod(4, -1, output);
    seg.add_lod(5, 0.3, output);
    seg.add_lod(6, 0.3, output);
    seg.add_lod(7, -2, output);
    BOOST_REQUIRE_EQUAL(output,
            "test\t2\t1\t1\t2\ntest\t2\t3\t3\t0.5\ntest\t2\t5\t6\t0.6\n");
    BOOST_REQUIRE_EQUAL(seg.length(), 0);
    BOOST_REQUIRE_EQUAL(pool_length(), pool);

    //trigger reversal twice
    seg.add_lod(1, 2, output);
    seg.add_lod(2, -1, output);
    seg.add_lod(3, 0.5, output);
    seg.add_lod(4, -1, output);
    seg.add_lod(5, 0.5, output);
    seg.add_lod(6, -1, output);
    seg.add_lod(7, 1.9, output);
    seg.add_lod(8, -1, output);
    seg.add_lod(9, 0.5, output);
    seg.add_lod(10, -3, output);
    BOOST_REQUIRE_EQUAL(output,
            "test\t2\t1\t1\t2\ntest\t2\t3\t3\t0.5\n"
            "test\t2\t5\t5\t0.5\ntest\t2\t7\t7\t1.9\n"
            "test\t2\t9\t9\t0.5\n");
    BOOST_REQUIRE_EQUAL(seg.length(), 0);
    BOOST_REQUIRE_EQUAL(pool_length(), pool);
}
