#include "test_pch.h"

#include <cmath>
#include <vector>

#include "core/utctime_utilities.h"
#include "core/time_axis.h"
#include "core/time_series_dd.h"

using shyft::core::no_utctime;
using shyft::core::utctime;
using shyft::core::seconds;
using shyft::core::utctime_0;
using std::numeric_limits;
const double eps = numeric_limits<double>::epsilon();
using shyft::time_series::dd::apoint_ts;
using shyft::time_series::dd::derivative_method;
using shyft::time_axis::generic_dt;
using shyft::time_series::ts_point_fx;
using std::vector;
using std::make_shared;
using std::isfinite;
//Helper
static inline utctime _t(int64_t t1970s) {return utctime{seconds(t1970s)};}

static void tc_ts_linear_derivative(generic_dt ta) {
    apoint_ts f(ta,vector<double>{1,1,2,3,-1.0,5.0},ts_point_fx::POINT_INSTANT_VALUE);
    auto d_f=f.derivative();
    FAST_CHECK_EQ(d_f.size(),f.size());
    FAST_CHECK_EQ(d_f.value(0),doctest::Approx(0.0));
    FAST_CHECK_EQ(d_f.value(1),doctest::Approx(0.1));
    FAST_CHECK_EQ(d_f.value(2),doctest::Approx(0.1));
    FAST_CHECK_EQ(d_f.value(3),doctest::Approx(-0.4));
    FAST_CHECK_EQ(d_f(f.time(3)+seconds(5)),doctest::Approx(-0.4));
    FAST_CHECK_EQ(d_f.value(4),doctest::Approx(0.6));
    FAST_CHECK_EQ(false,isfinite(d_f.value(5)));
    FAST_CHECK_EQ(false,isfinite(d_f(f.time(5))));
    FAST_CHECK_EQ(true,isfinite(d_f(f.time(5)-seconds(1) )));
    auto v=d_f.values();
    FAST_CHECK_EQ(v.size(),f.size());
    FAST_CHECK_EQ(v[0],doctest::Approx(0.0));
    FAST_CHECK_EQ(v[1],doctest::Approx(0.1));
    FAST_CHECK_EQ(v[2],doctest::Approx(0.1));
    FAST_CHECK_EQ(v[3],doctest::Approx(-0.4));
    FAST_CHECK_EQ(v[4],doctest::Approx(0.6));
    FAST_CHECK_EQ(false,isfinite(v[5]));
}

static void tc_ts_stair_case_derivative(generic_dt ta) {
    apoint_ts f(ta,vector<double>{1,1,2,shyft::nan,-1.0,5.0},ts_point_fx::POINT_AVERAGE_VALUE);
    auto d_f=f.derivative();
    FAST_CHECK_EQ(d_f.size(),f.size());
    FAST_CHECK_EQ(d_f.value(0),doctest::Approx(0.0));
    FAST_CHECK_EQ(d_f.value(1),doctest::Approx((2-1)/20.0));
    FAST_CHECK_EQ(d_f.value(2),doctest::Approx((2-1)/20.0));
    FAST_CHECK_EQ(false,isfinite(d_f.value(3)));
    FAST_CHECK_EQ(false,isfinite(d_f(f.time(3)+seconds(5) )));
    FAST_CHECK_EQ(d_f.value(4),doctest::Approx((5.0-(-1.0))/20.0));
    FAST_CHECK_EQ(d_f.value(5),doctest::Approx((5.0-(-1.0))/20.0));
    auto v=d_f.values();
    FAST_CHECK_EQ(v[0],doctest::Approx(0.0));
    FAST_CHECK_EQ(v[1],doctest::Approx((2-1)/20.0));
    FAST_CHECK_EQ(v[2],doctest::Approx((2-1)/20.0));
    FAST_CHECK_EQ(false,isfinite(v[3]));
    FAST_CHECK_EQ(v[4],doctest::Approx((5.0-(-1.0))/20.0));
    FAST_CHECK_EQ(v[5],doctest::Approx((5.0-(-1.0))/20.0));
}

static void tc_ts_stair_case_fwd_derivative(generic_dt ta) {
                                                            //0 1 2      3      4   5
    apoint_ts f(ta,vector<double>{1,1,2,shyft::nan,-1.0,5.0},ts_point_fx::POINT_AVERAGE_VALUE);
    auto d_f=f.derivative(derivative_method::forward_diff);
    FAST_CHECK_EQ(d_f.size(),f.size());
    FAST_CHECK_EQ(d_f.value(0),doctest::Approx(0.0));
    FAST_CHECK_EQ(d_f.value(1),doctest::Approx((2-1)/10.0));
    FAST_CHECK_EQ(d_f.value(2),doctest::Approx((2-2)/10.0));
    FAST_CHECK_EQ(false,isfinite(d_f.value(3)));
    FAST_CHECK_EQ(false,isfinite(d_f(f.time(3)+ seconds(5) )));
    FAST_CHECK_EQ(d_f.value(4),doctest::Approx((5.0-(-1.0))/10.0));
    FAST_CHECK_EQ(d_f.value(5),doctest::Approx((5.0-(5.0))/10.0));
    auto v=d_f.values();
    FAST_CHECK_EQ(v[0],doctest::Approx(0.0));
    FAST_CHECK_EQ(v[1],doctest::Approx((2-1)/10.0));
    FAST_CHECK_EQ(v[2],doctest::Approx((2-2)/10.0));
    FAST_CHECK_EQ(false,isfinite(v[3]));
    FAST_CHECK_EQ(v[4],doctest::Approx((5.0-(-1.0))/10.0));
    FAST_CHECK_EQ(v[5],doctest::Approx((5.0-(5.0))/10.0));
}

static void tc_ts_stair_case_bwd_derivative(generic_dt ta) {
                                                             //0 1 2      3      4   5
    apoint_ts f(ta,vector<double>{1,1,2,shyft::nan,-1.0,5.0},ts_point_fx::POINT_AVERAGE_VALUE);
    auto d_f=f.derivative(derivative_method::backward_diff);
    FAST_CHECK_EQ(d_f.size(),f.size());
    FAST_CHECK_EQ(d_f.value(0),doctest::Approx(0.0));
    FAST_CHECK_EQ(d_f.value(1),doctest::Approx((1-1)/10.0));
    FAST_CHECK_EQ(d_f.value(2),doctest::Approx((2-1)/10.0));
    FAST_CHECK_EQ(false,isfinite(d_f.value(3)));
    FAST_CHECK_EQ(false,isfinite(d_f(f.time(3)+seconds(5) )));
    FAST_CHECK_EQ(d_f.value(4),doctest::Approx(0.0));
    FAST_CHECK_EQ(d_f.value(5),doctest::Approx((5.0-(-1.0))/10.0));
    auto v=d_f.values();
    FAST_CHECK_EQ(v[0],doctest::Approx(0.0));
    FAST_CHECK_EQ(v[1],doctest::Approx((1-1)/10.0));
    FAST_CHECK_EQ(v[2],doctest::Approx((2-1)/10.0));
    FAST_CHECK_EQ(false,isfinite(v[3]));
    FAST_CHECK_EQ(v[4],doctest::Approx(0.0));
    FAST_CHECK_EQ(v[5],doctest::Approx((5.0-(-1.0))/10.0));
}

static void tc_ts_stair_case_one_value(derivative_method dm) {
    apoint_ts f1(generic_dt(utctime_0,seconds(10),1),vector<double>{0},ts_point_fx::POINT_AVERAGE_VALUE);
    apoint_ts f2(generic_dt(utctime_0,seconds(10),1),vector<double>{shyft::nan},ts_point_fx::POINT_AVERAGE_VALUE);
    auto d1=f1.derivative(dm);
    auto d2=f2.derivative(dm);
    auto v1=d1.values();
    auto v2=d2.values();
    FAST_CHECK_EQ(d1.value(0),doctest::Approx(0.0));
    FAST_CHECK_EQ(v1[0],doctest::Approx(0.0));
    FAST_CHECK_EQ(false,isfinite(d2.value(0)));
    FAST_CHECK_EQ(false,isfinite(v2[0]));
}

static void tc_ts_linear_one_segment(derivative_method dm) {
    apoint_ts f1(generic_dt(utctime_0,seconds(1),2),vector<double>{0,1.0},ts_point_fx::POINT_AVERAGE_VALUE);
    apoint_ts f2(generic_dt(utctime_0,seconds(1),2),vector<double>{shyft::nan,1.0},ts_point_fx::POINT_AVERAGE_VALUE);
    auto d1=f1.derivative(dm);
    auto d2=f2.derivative(dm);
    auto v1=d1.values();
    auto v2=d2.values();
    FAST_CHECK_EQ(d1.value(0),doctest::Approx(1.0/2));
    FAST_CHECK_EQ(v1[0],doctest::Approx(1.0/2));
    FAST_CHECK_EQ(false,isfinite(d2.value(0)));
    FAST_CHECK_EQ(false,isfinite(v2[0]));
}

TEST_SUITE("time_series") {

    generic_dt ta1(utctime_0,seconds(10),6);
    generic_dt ta2(vector<utctime>{_t(0),_t(10),_t(20),_t(30),_t(40),_t(50)},_t(6*10));

    TEST_CASE("ts_linear_derivative") {
        tc_ts_linear_derivative(ta1);
        tc_ts_linear_derivative(ta2);
    }

    TEST_CASE("ts_stair_case_derivative") {
        tc_ts_stair_case_derivative(ta1);
        tc_ts_stair_case_derivative(ta2);
    }
    TEST_CASE("ts_stair_case_fwd_derivative") {
        tc_ts_stair_case_fwd_derivative(ta1);
        tc_ts_stair_case_fwd_derivative(ta2);
    }
    TEST_CASE("ts_stair_case_bwd_derivative") {
        tc_ts_stair_case_bwd_derivative(ta1);
        tc_ts_stair_case_bwd_derivative(ta2);
    }
    TEST_CASE("ts_stair_case_one_value") {
        tc_ts_stair_case_one_value(derivative_method::default_diff);
        tc_ts_stair_case_one_value(derivative_method::forward_diff);
        tc_ts_stair_case_one_value(derivative_method::backward_diff);
    }
    TEST_CASE("ts_linear_one_segment") {
        tc_ts_linear_one_segment(derivative_method::default_diff);
    }
    TEST_CASE("derivative_ts_index_of") {
        using shyft::nan;
        apoint_ts a(ta1,1.0,ts_point_fx::POINT_AVERAGE_VALUE);
        auto b = a.derivative();
        for(size_t i=0;i<ta1.size();++i) {
            FAST_CHECK_EQ(b.index_of(ta1.time(i)),a.index_of(ta1.time(i)));
            FAST_CHECK_EQ(b.time(i),b.time(i));
        }
    }
    TEST_CASE("derivative_inside_issue") {
        /* https://github.com/statkraft/shyft/issues/352
         * The following code reproduces the issue. 
         * As seen from the last statement, 
         * orig+orig_derivative_inside_inf contains values [2. 2. 2. 2.], 
         * but 
         * (orig+orig_derivative_inside_inf).inside evaluates all values as NaN. 
         * Why?
         * Ans: There was a bug in derivative.index_of, which is now fixed, and covered with the above test.
         *      The original test-code, in c++ is shown below, with correct wanted behaviour
         */
        using shyft::nan;
        apoint_ts orig(ta1,1.0,ts_point_fx::POINT_AVERAGE_VALUE);
        auto a = orig.derivative().inside(nan,nan,nan,1.0,0.0);
        auto b = (orig + orig.derivative().inside(nan,nan,nan,1.0,0.0)).inside(nan,nan,nan,1.0,0.0);
        auto a_v = a.values();
        auto b_v = b.values();
        FAST_CHECK_EQ(a_v.size(),ta1.size());
        for(auto&v:b_v)
            FAST_CHECK_EQ(v,doctest::Approx(1.0));
    }
}
