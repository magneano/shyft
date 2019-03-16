#include "test_pch.h"
#define _USE_MATH_DEFINES
#include "core/time_series_dd.h"

using shyft::time_series::dd::apoint_ts;
using time_axis=shyft::time_series::dd::gta_t;
using shyft::core::utctime;
using shyft::core::calendar;
using shyft::core::deltahours;
using shyft::time_series::ts_point_fx;
using std::vector;
using v_=vector<double>;
using t_=vector<utctime>;
using shyft::core::no_utctime;

using shyft::time_series::dd::extend_ts_fill_policy;
using shyft::time_series::dd::extend_ts_split_policy;
using std::isfinite;

 static inline utctime t(double x){return shyft::core::from_seconds(x);};

namespace test::time_series {
    void is_expected_result (apoint_ts const&ts,time_axis const& expect_ta,v_ const&expect_v) {
        FAST_CHECK_EQ(ts.time_axis(),expect_ta);
        auto v=ts.values();
        FAST_REQUIRE_EQ(v.size(),expect_v.size());
        for(size_t i=0;i<v.size();++i) {
            if(isfinite(expect_v[i])) {
                FAST_CHECK_EQ(v[i],doctest::Approx(expect_v[i]));
                FAST_CHECK_EQ(ts.value(i),doctest::Approx(expect_v[i]));
                FAST_CHECK_EQ(ts(ts.time(i)),doctest::Approx(expect_v[i]));
            } else {
                FAST_CHECK_EQ(isfinite(v[i]),false);
                FAST_CHECK_EQ(isfinite(ts.value(i)),false);
                FAST_CHECK_EQ(isfinite(ts(ts.time(i))),false);
            }
        }
    }
    
    #define TST_CASE(x)
    /** test to verify results, given two time-series a and b, formulated using different type of time-axis */
    void verify_case_with_gap (apoint_ts const &a, apoint_ts const &b) {
        TST_CASE("ts_extend_gap_lhs_last") {
                is_expected_result (
                    a.extend(b,extend_ts_split_policy::EPS_LHS_LAST,extend_ts_fill_policy::EPF_NAN,no_utctime,shyft::nan),
                    time_axis{t_{t(0),t(1),t(2),t(3),t(4),t(5)}}, v_{0,1,shyft::nan,3,4}
                );
                is_expected_result (
                    a.extend(b,extend_ts_split_policy::EPS_LHS_LAST,extend_ts_fill_policy::EPF_LAST,no_utctime,shyft::nan),
                    time_axis{t_{t(0),t(1),t(2),t(3),t(4),t(5)}}, v_{0,1,1,3,4}
                );
                is_expected_result (
                    a.extend(b,extend_ts_split_policy::EPS_LHS_LAST,extend_ts_fill_policy::EPF_FILL,no_utctime,9.0),
                    time_axis{t_{t(0),t(1),t(2),t(3),t(4),t(5)}}, v_{0,1,9.0,3,4}
                );
        }
        TST_CASE("ts_extend_gap_rhs_first") {
                is_expected_result (
                    a.extend(b,extend_ts_split_policy::EPS_RHS_FIRST,extend_ts_fill_policy::EPF_NAN,no_utctime,shyft::nan),
                    time_axis{t_{t(0),t(1),t(2),t(3),t(4),t(5)}}, v_{0,1,shyft::nan,3,4}
                );
                is_expected_result (
                    a.extend(b,extend_ts_split_policy::EPS_RHS_FIRST,extend_ts_fill_policy::EPF_LAST,no_utctime,shyft::nan),
                    time_axis{t_{t(0),t(1),t(2),t(3),t(4),t(5)}}, v_{0,1,1,3,4}
                );
                is_expected_result (
                    a.extend(b,extend_ts_split_policy::EPS_RHS_FIRST,extend_ts_fill_policy::EPF_FILL,no_utctime,9.0),
                    time_axis{t_{t(0),t(1),t(2),t(3),t(4),t(5)}}, v_{0,1,9.0,3,4}
                );
        }
        TST_CASE("ts_extend_gap_value") {
                is_expected_result (
                    a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_NAN,t(3),shyft::nan),
                    time_axis{t_{t(0),t(1),t(2),t(3),t(4),t(5)}}, v_{0,1,shyft::nan,3,4}
                );
                is_expected_result (
                    a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_LAST,t(3),shyft::nan),
                    time_axis{t_{t(0),t(1),t(2),t(3),t(4),t(5)}}, v_{0,1,1,3,4}
                );
                is_expected_result (
                    a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_FILL,t(3),9.0),
                    time_axis{t_{t(0),t(1),t(2),t(3),t(4),t(5)}}, v_{0,1,9.0,3,4}
                );
        }
        TST_CASE("ts_extend_t_value_at_end_of_lhs") {
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_NAN,t(2),shyft::nan),
                time_axis{t_{t(0),t(1),t(2),t(3),t(4),t(5)}}, v_{0,1,shyft::nan,3,4}
            );
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_LAST,t(2),shyft::nan),
            time_axis{t_{t(0),t(1),t(2),t(3),t(4),t(5)}}, v_{0,1,1,3,4}
            );
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_FILL,t(2),9.0),
                time_axis{t_{t(0),t(1),t(2),t(3),t(4),t(5)}}, v_{0,1.0,9.0,3,4}
            );
        }
        TST_CASE("ts_extend_t_value_inside_lhs") {
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_NAN,t(1.9),shyft::nan),
                time_axis{t_{t(0),t(1),t(3),t(4),t(5)}}, v_{0,1,3,4}
            );
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_LAST,t(1.9),shyft::nan),
            time_axis{t_{t(0),t(1),t(3),t(4),t(5)}}, v_{0,1,3,4}
            );
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_FILL,t(1.9),9.0),
                time_axis{t_{t(0),t(1),t(3),t(4),t(5)}}, v_{0,1.0,3,4}
            );
        }        
    }
    
    void verify_case_with_overlap (apoint_ts const &a, apoint_ts const &b) {
        TST_CASE("ts_extend_overlap_lhs_last") {
            //apoint_ts a(time_axis{t_{t(0),t(1),t(2)}},v_{0,1},ts_point_fx::POINT_AVERAGE_VALUE);
            //apoint_ts b(time_axis{t_{t(1),t(2),t(3)}},v_{1.1,2.1},ts_point_fx::POINT_AVERAGE_VALUE);
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_LHS_LAST,extend_ts_fill_policy::EPF_NAN,no_utctime,shyft::nan),
                time_axis{t_{t(0),t(1),t(2),t(3)}}, v_{0,1,2.1}
            );
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_LHS_LAST,extend_ts_fill_policy::EPF_LAST,no_utctime,shyft::nan),
                time_axis{t_{t(0),t(1),t(2),t(3)}}, v_{0,1,2.1}
            );
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_LHS_LAST,extend_ts_fill_policy::EPF_FILL,no_utctime,9.0),
                time_axis{t_{t(0),t(1),t(2),t(3)}}, v_{0,1,2.1}
            );
        }
        TST_CASE("ts_extend_overlap_rhs_first") {
            //apoint_ts a(time_axis{t_{t(0),t(1),t(2)}},v_{0,1},ts_point_fx::POINT_AVERAGE_VALUE);
            //apoint_ts b(time_axis{t_{t(1),t(2),t(3)}},v_{1.1,2.1},ts_point_fx::POINT_AVERAGE_VALUE);
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_RHS_FIRST,extend_ts_fill_policy::EPF_NAN,no_utctime,shyft::nan),
                time_axis{t_{t(0),t(1),t(2),t(3)}}, v_{0,1.1,2.1}
            );
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_RHS_FIRST,extend_ts_fill_policy::EPF_LAST,no_utctime,shyft::nan),
                time_axis{t_{t(0),t(1),t(2),t(3)}}, v_{0,1.1,2.1}
            );
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_RHS_FIRST,extend_ts_fill_policy::EPF_FILL,no_utctime,9.0),
                time_axis{t_{t(0),t(1),t(2),t(3)}}, v_{0,1.1,2.1}
            );
        }
        TST_CASE("ts_extend_overlap_t_value_at_end_of_lhs") {
            //apoint_ts a(time_axis{t_{t(0),t(1),t(2)}},v_{0,1},ts_point_fx::POINT_AVERAGE_VALUE);
            //apoint_ts b(time_axis{t_{t(1),t(2),t(3)}},v_{1.1,2.1},ts_point_fx::POINT_AVERAGE_VALUE);
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_NAN,t(2),shyft::nan),
                time_axis{t_{t(0),t(1),t(2),t(3)}}, v_{0,1,2.1}
            );
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_LAST,t(2),shyft::nan),
                time_axis{t_{t(0),t(1),t(2),t(3)}}, v_{0,1,2.1}
            );
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_FILL,t(2),9.0),
                time_axis{t_{t(0),t(1),t(2),t(3)}}, v_{0,1,2.1}
            );
        }
        TST_CASE("ts_extend_overlap_t_value_inside_lhs") {
            //apoint_ts a(time_axis{t_{t(0),t(1),t(2)}},v_{0,1},ts_point_fx::POINT_AVERAGE_VALUE);
            //apoint_ts b(time_axis{t_{t(1),t(2),t(3)}},v_{1.1,2.1},ts_point_fx::POINT_AVERAGE_VALUE);
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_NAN,t(1.9),shyft::nan),
                time_axis{t_{t(0),t(1),t(2),t(3)}}, v_{0,1.1,2.1}
            );
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_LAST,t(1.9),shyft::nan),
                time_axis{t_{t(0),t(1),t(2),t(3)}}, v_{0,1.1,2.1}
            );
            is_expected_result (
                a.extend(b,extend_ts_split_policy::EPS_VALUE,extend_ts_fill_policy::EPF_FILL,t(1.9),9.0),
                time_axis{t_{t(0),t(1),t(2),t(3)}}, v_{0,1.1,2.1}
            );
        }          
    }
}
using namespace test::time_series;

TEST_SUITE("time_series_extend") {
    TEST_CASE("ts_extend_gap_f_f") {
        apoint_ts a(time_axis{t(0),t(1),2},v_{0,1},ts_point_fx::POINT_AVERAGE_VALUE);
        apoint_ts b(time_axis{t(3),t(1),2},v_{3.0,4.0},ts_point_fx::POINT_AVERAGE_VALUE);
        verify_case_with_gap (a,b);
    }
    TEST_CASE("ts_extend_gap_f_p") {
        apoint_ts a(time_axis{t(0),t(1),2},v_{0,1},ts_point_fx::POINT_AVERAGE_VALUE);
        apoint_ts b(time_axis{t_{t(3),t(4),t(5)}},v_{3.0,4.0},ts_point_fx::POINT_AVERAGE_VALUE);
        verify_case_with_gap (a,b);
    }
    TEST_CASE("ts_extend_gap_f") {
        apoint_ts a(time_axis{t_{t(0),t(1),t(2)}},v_{0,1},ts_point_fx::POINT_AVERAGE_VALUE);
        apoint_ts b(time_axis{t(3),t(1),2},v_{3.0,4.0},ts_point_fx::POINT_AVERAGE_VALUE);
        verify_case_with_gap (a,b);
    }
    TEST_CASE("ts_extend_gap_p_p") {
        apoint_ts a(time_axis{t_{t(0),t(1),t(2)}},v_{0,1},ts_point_fx::POINT_AVERAGE_VALUE);
        apoint_ts b(time_axis{t_{t(3),t(4),t(5)}},v_{3.0,4.0},ts_point_fx::POINT_AVERAGE_VALUE);
        verify_case_with_gap (a,b);
    }
    TEST_CASE("ts_extend_overlap_p_p"){
        apoint_ts a(time_axis{t_{t(0),t(1),t(2)}},v_{0,1},ts_point_fx::POINT_AVERAGE_VALUE);
        apoint_ts b(time_axis{t_{t(1),t(2),t(3)}},v_{1.1,2.1},ts_point_fx::POINT_AVERAGE_VALUE);
        verify_case_with_overlap(a,b);
    }
    TEST_CASE("ts_extend_overlap_p_f"){
        apoint_ts a(time_axis{t_{t(0),t(1),t(2)}},v_{0,1},ts_point_fx::POINT_AVERAGE_VALUE);
        apoint_ts b(time_axis{t(1),t(1),2},v_{1.1,2.1},ts_point_fx::POINT_AVERAGE_VALUE);
        verify_case_with_overlap(a,b);
    }
    TEST_CASE("ts_extend_overlap_f_p"){
        apoint_ts a(time_axis{t(0),t(1),2},v_{0,1},ts_point_fx::POINT_AVERAGE_VALUE);
        apoint_ts b(time_axis{t_{t(1),t(2),t(3)}},v_{1.1,2.1},ts_point_fx::POINT_AVERAGE_VALUE);
        verify_case_with_overlap(a,b);
    }
    TEST_CASE("ts_extend_overlap_f_f"){
        apoint_ts a(time_axis{t(0),t(1),2},v_{0,1},ts_point_fx::POINT_AVERAGE_VALUE);
        apoint_ts b(time_axis{t(1),t(1),2},v_{1.1,2.1},ts_point_fx::POINT_AVERAGE_VALUE);
        verify_case_with_overlap(a,b);
    }
    
}
