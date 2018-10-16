

#include "test_pch.h"
#include "core/pt_hps_k.h"
#include "core/cell_model.h"
#include "core/pt_hps_k_cell_model.h"
#include "mocks.h"
#include "core/time_series.h"
#include "core/utctime_utilities.h"

// Some typedefs for clarity
using namespace shyft::core;
using namespace shyft::time_series;
using namespace shyft::core::pt_hps_k;

using namespace shyfttest::mock;
using namespace shyfttest;

namespace pt = shyft::core::priestley_taylor;
namespace hps = shyft::core::hbv_physical_snow;
namespace gm = shyft::core::glacier_melt;
namespace kr = shyft::core::kirchner;
namespace ae = shyft::core::actual_evapotranspiration;
namespace pc = shyft::core::precipitation_correction;
namespace ta = shyft::time_axis;
typedef TSPointTarget<ta::point_dt> catchment_t;

namespace shyfttest {
    namespace mock {
        // need specialization for pthpsk_response_t above
        template<> template<>
        void ResponseCollector<ta::fixed_dt>::collect<response>(size_t idx, const response& response) {
            _snow_output.set(idx, response.hps.outflow);
        }
        template <> template <>
        void DischargeCollector<ta::fixed_dt>::collect<response>(size_t idx, const response& response) {
            // q_avg is given in mm, so compute the totals
            avg_discharge.set(idx, destination_area*response.kirchner.q_avg / 1000.0 / 3600.0);
        }
    };
}; // End namespace shyfttest
TEST_SUITE("pt_hps_k") {
TEST_CASE("test_call_stack") {
    xpts_t temp;
    xpts_t prec;
    xpts_t rel_hum;
    xpts_t wind_speed;
    xpts_t radiation;

    calendar cal;
    utctime t0 = cal.time(YMDhms(2014, 8, 1, 0, 0, 0));
    size_t n_ts_points = 3*24;
    utctimespan dt  = deltahours(1);
    utctime t1 = t0 + n_ts_points*dt;
    shyfttest::create_time_series(temp, prec, rel_hum, wind_speed, radiation, t0, dt, n_ts_points);

    auto model_dt = deltahours(24);
    vector<utctime> times;
    for (utctime i=t0; i <= t1; i += model_dt)
        times.emplace_back(i);
    ta::fixed_dt time_axis(t0, dt, n_ts_points);
	ta::fixed_dt state_time_axis(t0, dt, n_ts_points + 1);
    // Initialize parameters
    std::vector<double> s = {1.0, 1.0, 1.0, 1.0, 1.0}; // Zero cv distribution of snow (i.e. even)
    std::vector<double> a = {0.0, 0.25, 0.5, 0.75, 1.0};

    pt::parameter pt_param;
    hps::parameter hps_param(s, a);
    ae::parameter ae_param;
    kr::parameter k_param;
    pc::parameter p_corr_param;

    // Initialize the state vectors
    kr::state kirchner_state {5.0};
    hps::state hps_state(vector<double>(5, 0.4), vector<double>(5, 0.0),30000.0, 10.0, 0.5);

    // Initialize response
    response run_response;

    // Initialize collectors
    shyfttest::mock::ResponseCollector<ta::fixed_dt> response_collector(1000*1000, time_axis);
    shyfttest::mock::StateCollector<ta::fixed_dt> state_collector(state_time_axis);

    state state {hps_state, kirchner_state};
    parameter parameter(pt_param, hps_param, ae_param, k_param, p_corr_param);
    geo_cell_data geo_cell_data;
    pt_hps_k::run<direct_accessor, response>(geo_cell_data, parameter, time_axis,0,0, temp,
                                              prec, wind_speed, rel_hum, radiation, state,
                                              state_collector, response_collector);

    auto snow_swe = response_collector.snow_swe();
    for (size_t i = 0; i < snow_swe.size(); ++i)
        TS_ASSERT(std::isfinite(snow_swe.get(i).v) && snow_swe.get(i).v >= 0);
}

TEST_CASE("pt_hps_k_lake_reservoir_response") {
    calendar cal;
    utctime t0 = cal.time(2014, 8, 1, 0, 0, 0);
    utctimespan dt=deltahours(1);
    const int n=50;// need to run some steps to observe kirchner response
    ta::fixed_dt tax(t0,dt,n);
	ta::fixed_dt tax_state(t0, dt, n + 1);
    pt::parameter pt_param;
    hps::parameter gs_param;
    ae::parameter ae_param;
    kr::parameter k_param;
    pc::parameter p_corr_param;
    parameter parameter{pt_param, gs_param, ae_param, k_param, p_corr_param};

    pts_t temp(tax,-15.0,POINT_AVERAGE_VALUE);// freezing cold
    pts_t prec(tax,3.0,POINT_AVERAGE_VALUE);prec.set(0,0.0);// rain except 1.step we use to get initial kirchner response
    pts_t rel_hum(tax,0.8,POINT_AVERAGE_VALUE);
    pts_t wind_speed(tax,2.0,POINT_AVERAGE_VALUE);
    pts_t radiation(tax,300.0,POINT_AVERAGE_VALUE);

    kr::state kirchner_state{1};// 1 mm 
    hps::state gs_state;
    //gs_state.lwc=100.0;
    //gs_state.acc_melt=100.0; 
    
    state s0{gs_state, kirchner_state};// need a equal state for the second run
    state s1{gs_state, kirchner_state};
    
    pt_hps_k::state_collector sc;
    pt_hps_k::all_response_collector rc;
    const double cell_area=1000*1000;
    sc.collect_state=true;
    sc.initialize(tax_state,0,0,cell_area);
    rc.initialize(tax,0,0,cell_area);
    geo_cell_data gcd(geo_point(1000,1000,100));
    land_type_fractions ltf(0.0,0.2,0.3,0.0,0.5);// 0.2 lake, 0.3 reservoir , 0.5 unspec
    gcd.set_land_type_fractions(ltf);
    
    parameter.msp.reservoir_direct_response_fraction=0.0;// all rsv goes to kirchner
    pt_hps_k::run<direct_accessor,pt_hps_k::response>(gcd,parameter,tax,0,0,temp,prec,wind_speed,rel_hum,radiation,s0,sc,rc);
    FAST_CHECK_EQ(rc.avg_discharge.value(0), doctest::Approx(0.266).epsilon(0.01)); // first with 0 precip, nothing should happen
    FAST_CHECK_EQ(rc.avg_discharge.value(n-1), doctest::Approx( 0.5*mmh_to_m3s(prec.value(1),cell_area)).epsilon(0.01));
    
    parameter.msp.reservoir_direct_response_fraction=1.0;// all rsv goes directly to output, lake goes to kirchner
    sc.initialize(tax_state,0,0,cell_area);
    rc.initialize(tax,0,0,cell_area);
    pt_hps_k::run<direct_accessor,pt_hps_k::response>(gcd,parameter,tax,0,0,temp,prec,wind_speed,rel_hum,radiation,s1,sc,rc);
    FAST_CHECK_EQ(rc.avg_discharge.value(0), doctest::Approx(0.266*0.7).epsilon(0.01)); // first with 0 precip, nothing should happen
    auto expected_1 = 0.266+0.3*0.5*mmh_to_m3s(prec.value(1),cell_area);
    FAST_CHECK_EQ(rc.avg_discharge.value(1), doctest::Approx(expected_1).epsilon(0.05)); // precip on rsv direct effect
    auto expected_2 = 0.2*mmh_to_m3s(prec.value(1),cell_area)*(1.0-0.3)+0.3*mmh_to_m3s(prec.value(1),cell_area);
    FAST_CHECK_EQ(rc.avg_discharge.value(n-1), doctest::Approx(expected_2 ).epsilon(0.01));
}

}
