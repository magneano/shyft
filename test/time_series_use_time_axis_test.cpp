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


TEST_SUITE("time_series") {
    TEST_CASE("ts_use_time_axis") {
        auto t0=calendar().time(2010,1,1);
        vector<double> v{0.0,1.0,2.0,3.0,4.0,5.0};
        auto dt=deltahours(1);
        auto n=v.size();
        time_axis ta{t0,dt,n};
        time_axis to{t0,2*dt,n/2};
        apoint_ts a(ta,v,ts_point_fx::POINT_AVERAGE_VALUE);
        apoint_ts o(to,shyft::nan,ts_point_fx::POINT_AVERAGE_VALUE);
        auto r = a.use_time_axis_from(o);
        FAST_CHECK_EQ(r.time_axis(),o.time_axis());
        auto r_value=r.values();
        FAST_CHECK_EQ(r_value[0],doctest::Approx(v[0]));
        FAST_CHECK_EQ(r_value[1],doctest::Approx(v[2]));
        FAST_CHECK_EQ(r_value[2],doctest::Approx(v[4]));
        
        FAST_CHECK_EQ(r.value(0),doctest::Approx(v[0]));
        FAST_CHECK_EQ(r.value(1),doctest::Approx(v[2]));
        FAST_CHECK_EQ(r.value(2),doctest::Approx(v[4]));
        
        FAST_CHECK_EQ(r(to.time(0)),doctest::Approx(v[0]));
        FAST_CHECK_EQ(r(to.time(1)),doctest::Approx(v[2]));
        FAST_CHECK_EQ(r(to.time(2)),doctest::Approx(v[4]));
        
        FAST_CHECK_EQ(r.point_interpretation(),a.point_interpretation());
        
        auto r2 = a.use_time_axis_from(a);
        FAST_CHECK_EQ(r2,a);
    }
   
}
