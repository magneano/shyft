#include "test_pch.h"
#define _USE_MATH_DEFINES
#include "core/time_series_dd.h"

using shyft::time_series::dd::apoint_ts;
using shyft::time_series::dd::ats_vector;
using time_axis=shyft::time_series::dd::gta_t;
using shyft::core::utctime;
using shyft::core::utcperiod;
using shyft::core::calendar;
using shyft::core::deltahours;
using shyft::time_series::ts_point_fx;
using std::vector;

using shyft::time_series::dd::clip_to_period;


TEST_SUITE("time_series") {
    TEST_CASE("ts_clip_empty") {
        apoint_ts a;
        auto b=clip_to_period(a,utcperiod(utctime(-100),utctime(0)));// no overlap -> empty ts.
        FAST_CHECK_EQ(a.size(),0u);
        FAST_CHECK_EQ(b.size(),0u);
    }
    TEST_CASE("ts_clip_to_period_stair_case") {
        apoint_ts a(time_axis(utctime(0),utctime(10),5),vector<double>{1,2,3,4.,5.},ts_point_fx::POINT_AVERAGE_VALUE);
        auto b=clip_to_period(a,utcperiod(utctime(0),utctime(10)*5));// equal or full overlap
        FAST_CHECK_EQ(b.time_axis(),a.time_axis());
        b=clip_to_period(a,utcperiod(utctime(-100),utctime(0)));// no overlap -> empty ts.
        FAST_CHECK_EQ(b.time_axis().size(),0u);
        // partial overlap, ts have a part that ends into period p:
        b=clip_to_period(a,utcperiod(utctime(-100),utctime(1)));// barely one first value.
        FAST_REQUIRE_EQ(b.time_axis().size(),1u);// one point
        FAST_CHECK_EQ(b.value(0),doctest::Approx(1.0));
        FAST_CHECK_EQ(b.time_axis(),time_axis{utctime(0),utctime(10),1});
        // partial overlap, ts begins in period, then extend after
        b=clip_to_period(a,utcperiod(utctime(3*10+1),utctime(6*10)));// two last values
        FAST_REQUIRE_EQ(b.time_axis().size(),2u);
        FAST_CHECK_EQ(b.value(0),doctest::Approx(4.0));
        FAST_CHECK_EQ(b.value(1),doctest::Approx(5.0));
        FAST_CHECK_EQ(b.time_axis(),time_axis{utctime(30),utctime(10),2});
        // last case, ts cliped in the middle
        b=clip_to_period(a,utcperiod(utctime(1*10),utctime(2*10)));// mid value
        FAST_REQUIRE_EQ(b.time_axis().size(),1u);
        FAST_CHECK_EQ(b.value(0),doctest::Approx(2.0));
        FAST_CHECK_EQ(b.time_axis(),time_axis{utctime(10),utctime(10),1});
        
    }
    TEST_CASE("ts_clip_to_period_linear") {
        apoint_ts a(time_axis(utctime(0),utctime(10),5),vector<double>{1,2,3,4.,5.},ts_point_fx::POINT_INSTANT_VALUE);
        auto b=clip_to_period(a,utcperiod(utctime(0),utctime(10)*5));// equal or full overlap
        FAST_CHECK_EQ(b.time_axis(),a.time_axis());
        b=clip_to_period(a,utcperiod(utctime(-100),utctime(0)));// no overlap -> empty ts.
        FAST_CHECK_EQ(b.time_axis().size(),0u);
        // partial overlap, ts have a part that ends into period p:
        b=clip_to_period(a,utcperiod(utctime(-100),utctime(1)));// barely one first value.
        FAST_REQUIRE_EQ(b.time_axis().size(),2u);// two points, because it needs those to eval. correct value
        FAST_CHECK_EQ(b.value(0),doctest::Approx(1.0));
        FAST_CHECK_EQ(b.value(1),doctest::Approx(2.0));
        
        FAST_CHECK_EQ(b.time_axis(),time_axis{utctime(0),utctime(10),2});
        // partial overlap, ts begins in period, then extend after
        b=clip_to_period(a,utcperiod(utctime(3*10-1),utctime(6*10)));// three(barely )wo last values
        FAST_REQUIRE_EQ(b.time_axis().size(),3u);
        FAST_CHECK_EQ(b.value(0),doctest::Approx(3.0));
        FAST_CHECK_EQ(b.value(1),doctest::Approx(4.0));
        FAST_CHECK_EQ(b.value(2),doctest::Approx(5.0));
        FAST_CHECK_EQ(b.time_axis(),time_axis{utctime(20),utctime(10),3});
        // last case, ts cliped in the middle
        b=clip_to_period(a,utcperiod(utctime(1*10+1),utctime(2*10)));// two mid values
        FAST_REQUIRE_EQ(b.time_axis().size(),2u);
        FAST_CHECK_EQ(b.value(0),doctest::Approx(2.0));
        FAST_CHECK_EQ(b.value(1),doctest::Approx(3.0));
        FAST_CHECK_EQ(b.time_axis(),time_axis{utctime(10),utctime(10),2});
    }
    TEST_CASE("ts_clip_to_period_linear_bp") {
        apoint_ts a(time_axis( // also test with break-point ts, because they are kind of exotic when it comes to logic.
                vector<utctime>{utctime(0),utctime(10),utctime(20),utctime(30),utctime(40)},utctime(50)),
                vector<double>{1,2,3,4.,5.},
                ts_point_fx::POINT_INSTANT_VALUE);
        // notice that time-axis compare are 'semantic' not exact, so one time-axis is equal to the other if they contains same periods
        auto b=clip_to_period(a,utcperiod(utctime(0),utctime(10)*5));// equal or full overlap
        FAST_CHECK_EQ(b.time_axis(),a.time_axis());
        b=clip_to_period(a,utcperiod(utctime(-100),utctime(0)));// no overlap -> empty ts.
        FAST_CHECK_EQ(b.time_axis().size(),0u);
        // partial overlap, ts have a part that ends into period p:
        b=clip_to_period(a,utcperiod(utctime(-100),utctime(1)));// barely one first value.
        FAST_REQUIRE_EQ(b.time_axis().size(),2u);// two points, because it needs those to eval. correct value
        FAST_CHECK_EQ(b.value(0),doctest::Approx(1.0));
        FAST_CHECK_EQ(b.value(1),doctest::Approx(2.0));
        
        FAST_CHECK_EQ(b.time_axis(),time_axis{utctime(0),utctime(10),2});
        // partial overlap, ts begins in period, then extend after
        b=clip_to_period(a,utcperiod(utctime(3*10-1),utctime(6*10)));// three(barely )wo last values
        FAST_REQUIRE_EQ(b.time_axis().size(),3u);
        FAST_CHECK_EQ(b.value(0),doctest::Approx(3.0));
        FAST_CHECK_EQ(b.value(1),doctest::Approx(4.0));
        FAST_CHECK_EQ(b.value(2),doctest::Approx(5.0));
        FAST_CHECK_EQ(b.time_axis(),time_axis{utctime(20),utctime(10),3});
        // last case, ts cliped in the middle
        b=clip_to_period(a,utcperiod(utctime(1*10+1),utctime(2*10)));// two last values(notice we tick one up at lower bound)
        FAST_REQUIRE_EQ(b.time_axis().size(),2u);
        FAST_CHECK_EQ(b.value(0),doctest::Approx(2.0));
        FAST_CHECK_EQ(b.value(1),doctest::Approx(3.0));
        FAST_CHECK_EQ(b.time_axis(),time_axis{utctime(10),utctime(10),2});
    }
    
    TEST_CASE("tsv_clip_to_period") {
        // we trust the tests above for the logic, here we only do vector-related sanity-checks
        ats_vector a;
        a.emplace_back(time_axis( // also test with break-point ts, because they are kind of exotic when it comes to logic.
                vector<utctime>{utctime(0),utctime(10),utctime(20),utctime(30),utctime(40)},utctime(50)),
                vector<double>{1,2,3,4.,5.},
                ts_point_fx::POINT_INSTANT_VALUE
        );
        a.emplace_back(apoint_ts{});
        a.emplace_back(time_axis(utctime(0),utctime(10),5),vector<double>{1,2,3,4.,5.},ts_point_fx::POINT_AVERAGE_VALUE);
        auto b = clip_to_period(a,utcperiod(utctime(1*10+1),utctime(2*10)));
        FAST_REQUIRE_EQ(b.size(),a.size());
        for(size_t i=0;i<b.size();++i)
            FAST_CHECK_LE(b[i].size(),a[i].size());// verify at least shorter time-axis
        ats_vector z;
        auto zb=clip_to_period(z,utcperiod(utctime(1*10+1),utctime(2*10)));
        FAST_REQUIRE_EQ(z.size(),zb.size());
    }   
}
