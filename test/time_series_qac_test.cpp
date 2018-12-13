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
    TEST_CASE("qac_speed") {
            // purpose to verify speed of qac(ts..).values
            // to be as fast as ts
            //--
            auto t0=calendar().time(2010,1,1);
            auto n= 24*365*10;
            auto dt=deltahours(1);
            apoint_ts a_1(time_axis(t0,dt,n),1.0,ts_point_fx::POINT_AVERAGE_VALUE);
            apoint_ts a_nan(time_axis(t0,dt,n),shyft::nan,ts_point_fx::POINT_AVERAGE_VALUE);
            apoint_ts z(time_axis(t0,dt*n,1),0.0,ts_point_fx::POINT_AVERAGE_VALUE);
            
            
            auto q_1=a_1.min_max_check_ts_fill(shyft::nan,shyft::nan,dt,z);
            auto t_0 =timing::now();
            auto v_1=q_1.values();
            auto t_1 =elapsed_us(t_0,timing::now());
            auto q_nan=a_nan.min_max_check_ts_fill(shyft::nan,shyft::nan,dt,z);
            auto t_2= timing::now();
            auto v_nan=q_nan.values();
            auto t_3= elapsed_us(t_2,timing::now());
            
            MESSAGE("No check time:"<<t_1<< ", fill_time:"<<t_3<<" [us]");
            FAST_CHECK_GT(t_3,t_1);// no fill should be faster than with-fill (approx 2x )
    }
    
    TEST_CASE("qac_ts_fill_algo") {
        auto   ta = time_axis(0, 1, 5);
        auto ts_src = apoint_ts(ta, vector<double>{1.0, -1.0, 2.0,shyft::nan, 4.0},ts_point_fx::POINT_AVERAGE_VALUE);
        auto  cts = apoint_ts(ta,vector<double>{1.0, 1.8, 2.0, 2.0, 4.0}, ts_point_fx::POINT_AVERAGE_VALUE);
        auto  ts_qac = ts_src.min_max_check_ts_fill(-10.0,10.0,300, cts);
        FAST_CHECK_EQ(ts_qac.value(3), doctest::Approx(2.0));
        ts_qac = ts_src.min_max_check_ts_fill( 0.0,10.0,300,cts);
        FAST_CHECK_EQ(ts_qac.value(1), doctest::Approx(1.8));  // -1 out, replaced with linear between
        FAST_CHECK_EQ(ts_qac.value(3), doctest::Approx(2.0));
        auto v=ts_qac.values();// ensure these are equal
        for(size_t i=0;i<ts_qac.size();++i) {
            FAST_CHECK_EQ(v[i],doctest::Approx(ts_qac.value(i)));
        }
        
    }
}
