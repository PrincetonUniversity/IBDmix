#define BOOST_TEST_MODULE FUIntroSegTest
#include <boost/test/included/unit_test.hpp>
#include "../FUIntroSeg_v2.1.cpp"
#include <math.h>

BOOST_AUTO_TEST_CASE(test_cal_lod)
{
    double result;
    double pb = 0.1;
    double aerr = 0.01;
    double merr = 0.02;
    double pa = 1-pb;
    double minesp = 1e-200;

    // 0, 0
    result = log10(((1-aerr)*(1-merr)+aerr*merr) / pa / (1-aerr*(1-aerr)));
    BOOST_REQUIRE_EQUAL(cal_lod(0, 0, pb, aerr, merr), result);

    // 0, 1
    result = log10(0.5 * ((1-aerr)*merr+aerr*(1-merr)) / pb / (1-aerr*(1-aerr)) +
            0.5 * ((1-aerr)*(1-merr)+aerr*merr) / pa / (1-aerr*(1-aerr)));
    BOOST_REQUIRE_EQUAL(cal_lod(0, 1, pb, aerr, merr), result);

//    if(source_gt == 0 && target_gt == 2)
//        result = log10(((1-aerr)*merr+aerr*(1-merr) + minesp) / pb / (1-aerr*(1-aerr) + minesp));
//    if(source_gt == 1 && target_gt == 0)
//        result = -log10(pa * (1+2*aerr*(1-aerr)));
//    if(source_gt == 1 && target_gt == 1)
//        result = -log10(2*pa*pb*(1+2*aerr*(1-aerr)));
//    if(source_gt == 1 && target_gt == 2)
//        result = -log10(pb * (1+2*aerr*(1-aerr)));
//    if(source_gt == 2 && target_gt == 0)
//        result = log10(((1-aerr)*merr+aerr*(1-merr) + minesp) / pa / (1-aerr*(1-aerr) + minesp));
//    if(source_gt == 2 && target_gt == 1)
//        result = log10(0.5*((1-aerr)*(1-merr)+aerr*merr) / pb / (1-aerr*(1-aerr)) + 0.5 * ((1-aerr)*merr+aerr*(1-merr)) / pa / (1-aerr*(1-aerr)));
//    if(source_gt == 2 && target_gt == 2)
//        result = log10(((1-aerr)*(1-merr)+aerr*merr) / pb / (1-aerr*(1-aerr)));
}
