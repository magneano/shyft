#include "test_pch.h"
#define _USE_MATH_DEFINES
#include "core/time_series_dd.h"

using shyft::time_series::dd::apoint_ts;
using time_axis=shyft::time_series::dd::gta_t;
using shyft::core::utctime;
using shyft::core::calendar;
using shyft::core::deltahours;
using shyft::core::seconds;
using shyft::time_series::accumulate_value;
using shyft::time_series::ts_point_fx;
using std::vector;

vector<double> accumulate_values(apoint_ts const&ts,time_axis const &ta) {
    utctime tsum;
    size_t ix_hint = 0;
    double a=0.0;
    vector<double> r; r.reserve(ta.size());
    auto lin= ts.point_interpretation() == ts_point_fx::POINT_INSTANT_VALUE;
    for(size_t i=0;i<ta.size();++i) {
        r.push_back(a);
        auto x=accumulate_value(ts, ta.period(i), ix_hint, tsum,lin);
        if(isfinite(x)) a+=x;
    }
    return r;
}

TEST_SUITE("time_series") {
    TEST_CASE("accumulate_ts") {
        auto t0=calendar().time(2010,1,1);
        vector<double> v{shyft::nan,1.0,shyft::nan,1.0,1.0,1.0};
        auto dt=seconds(1);
        auto n=v.size();
        time_axis ta{t0,dt,n};
        apoint_ts a(ta,v,ts_point_fx::POINT_AVERAGE_VALUE);
        auto r = a.accumulate(ta);
        auto r_value=r.values();
        FAST_CHECK_EQ(r_value[0],doctest::Approx(0.0));
        FAST_CHECK_EQ(r_value[1],doctest::Approx(0.0));
        FAST_CHECK_EQ(r_value[2],doctest::Approx(1.0));
        FAST_CHECK_EQ(r_value[3],doctest::Approx(1.0));
        FAST_CHECK_EQ(r_value[4],doctest::Approx(2.0));
        FAST_CHECK_EQ(r.value(0),doctest::Approx(0.0));
        FAST_CHECK_EQ(r.value(1),doctest::Approx(0.0));
        FAST_CHECK_EQ(r.value(2),doctest::Approx(1.0));
        FAST_CHECK_EQ(r.value(3),doctest::Approx(1.0));
        FAST_CHECK_EQ(r.value(4),doctest::Approx(2.0));
        
        FAST_CHECK_EQ(r.point_interpretation(),ts_point_fx::POINT_INSTANT_VALUE);
        
    }
    TEST_CASE("accumulate_ts_lin") {
        auto dt=seconds(1);
        vector<double> v{shyft::nan,1.0,1.0,shyft::nan,1.0,1.0};
        vector<utctime> t{0*dt,1*dt,2*dt,3*dt,4*dt,5*dt};        
        time_axis ta{t,t.back()+utctime{1}};
        apoint_ts a(ta,v,ts_point_fx::POINT_INSTANT_VALUE);
        auto r = a.accumulate(ta);
        auto r_value=r.values();
        FAST_CHECK_EQ(r_value[0],doctest::Approx(0.0));
        FAST_CHECK_EQ(r_value[1],doctest::Approx(0.0));
        FAST_CHECK_EQ(r_value[2],doctest::Approx(1.0));
        FAST_CHECK_EQ(r_value[3],doctest::Approx(1.0));
        FAST_CHECK_EQ(r_value[4],doctest::Approx(1.0));
        FAST_CHECK_EQ(r_value[5],doctest::Approx(2.0));
        
        FAST_CHECK_EQ(r.point_interpretation(),ts_point_fx::POINT_INSTANT_VALUE);
        
    }   
}
