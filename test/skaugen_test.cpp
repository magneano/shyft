#include "test_pch.h"
#include "core/skaugen.h"


using namespace shyft::core::skaugen;
using namespace shyft::core;
typedef calculator<parameter, state, response> SkaugenModel;
TEST_SUITE("skaugen") {
TEST_CASE("test_accumulation") {
    // Model parameters
    const double d_range = 113.0;
    const double unit_size = 0.1;
    const double alpha_0 = 40.77;
    const double max_water_fraction = 0.1;
    const double tx = 0.16;
    const double cx = 2.50;
    const double ts = 0.14;
    const double cfr = 0.01;
    parameter p(alpha_0, d_range, unit_size, max_water_fraction, tx, cx, ts, cfr);

    // Model state variables
    const double alpha = alpha_0;
    const double nu = alpha_0*unit_size;
    const double sca = 0.0;
    const double swe = 0.0;
    const double free_water = 0.0;
    const double residual = 0.0;
    const unsigned long nnn = 0;
    state s(nu, alpha, sca, swe, free_water, residual, nnn);

    // Model input
    utctimespan dt = seconds(60*60);
    double temp = -10.0;
    double prec = 10.0;
    double radiation = 0.0;
    double wind_speed = 0.0;
    std::vector<std::pair<double, double>> tp(10, std::pair<double, double>(temp, prec));

    // Accumulate snow
    SkaugenModel model;
    response r;
	for_each(tp.begin(), tp.end(), [&dt, &p, &model, &radiation, &wind_speed, &s, &r](std::pair<double, double> pair) {
            model.step(dt, p, pair.first, pair.second, radiation, wind_speed, s, r);
        });
    TS_ASSERT_DELTA(s.swe*s.sca, prec*tp.size(), 1.0e-6);
    TS_ASSERT_DELTA(s.sca, 1.0, 1.0e-6);
    TS_ASSERT(s.nu < alpha_0*unit_size);
}


TEST_CASE("test_melt") {
    // Model parameters
    const double d_range = 113.0;
    const double unit_size = 0.1;
    const double alpha_0 = 40.77;
    const double max_water_fraction = 0.1;
    const double tx = 0.16;
    const double cx = 2.50;
    const double ts = 0.14;
    const double cfr = 0.01;
    parameter p(alpha_0, d_range, unit_size, max_water_fraction, tx, cx, ts, cfr);

    // Model state variables
    const double alpha = alpha_0;
    const double nu = alpha_0*unit_size;
    const double sca = 0.0;
    const double swe = 0.0;
    const double free_water = 0.0;
    const double residual = 0.0;
    const unsigned long nnn = 0;
    state s(nu, alpha, sca, swe, free_water, residual, nnn);

    // Model input
    utctimespan dt = seconds(24*60*60);
    double temp = -10.0;
    double prec = 10.0;
    double radiation = 0.0;
    double wind_speed = 0.0;
    std::vector<std::pair<double, double>> tp(10, std::pair<double, double>(temp, prec));

    // Accumulate snow
    SkaugenModel model;
    response r;
	for_each(tp.begin(), tp.end(), [&dt, &p, &model, &radiation, &wind_speed, &s, &r](std::pair<double, double> pair) {
            model.step(dt, p, pair.first, pair.second, radiation, wind_speed, s, r);
        });

    const double total_water = s.swe*s.sca;
    double agg_outflow{0.0}; // For checking mass balance

    // Single melt event
    tp = std::vector<std::pair<double, double>>(1, std::pair<double, double>{10.0, 0.0}); // No precip, but 10.0 degrees for one day
	for_each(tp.begin(), tp.end(), [&dt, &p, &model, &radiation, &wind_speed, &s, &r, &agg_outflow](std::pair<double, double> pair) {
            model.step(dt, p, pair.first, pair.second, radiation, wind_speed, s, r);
            agg_outflow += r.outflow;
        });
    const double total_water_after_melt = s.sca*(s.swe + s.free_water);
    TS_ASSERT_LESS_THAN(total_water_after_melt , total_water); // Less water after melt due to runoff
    TS_ASSERT(r.outflow + s.free_water >= 1.0); // Some runoff or free water in snow
    TS_ASSERT_DELTA(r.outflow + s.sca*(s.free_water + s.swe), total_water, 1.0e-6);

    // One hundred melt events, that should melt everything
    tp = std::vector<std::pair<double, double>>(100, std::pair<double, double>(10.0, 0.0));
    for_each(tp.begin(), tp.end(), [&dt , &p, &model, &radiation, &wind_speed, &s, &r, &agg_outflow] (std::pair<double, double> pair) {
            model.step(dt, p, pair.first, pair.second, radiation, wind_speed, s, r);
            agg_outflow += r.outflow;
        });

    TS_ASSERT_DELTA(s.sca, 0.0, 1.0e-6);
    TS_ASSERT_DELTA(s.swe, 0.0, 1.0e-6);
    TS_ASSERT_DELTA(agg_outflow, total_water, 1.0e-10);
    TS_ASSERT_DELTA(s.alpha, alpha_0, 1.0e-6);
    TS_ASSERT_DELTA(s.nu, alpha_0*unit_size, 1.0e-6);
}


TEST_CASE("test_lwc") {
    // Model parameters
    const double d_range = 113.0;
    const double unit_size = 0.1;
    const double alpha_0 = 40.77;
    const double max_water_fraction = 0.1;
    const double tx = 0.16;
    const double cx = 2.50;
    const double ts = 0.14;
    const double cfr = 0.01;
    parameter p(alpha_0, d_range, unit_size, max_water_fraction, tx, cx, ts, cfr);

    // Model state variables
    const double alpha = alpha_0;
    const double nu = alpha_0*unit_size;
    const double sca = 0.0;
    const double swe = 0.0;
    const double free_water = 0.0;
    const double residual = 0.0;
    const unsigned long nnn = 0;
    state s(nu, alpha, sca, swe, free_water, residual, nnn);

    // Model input
    utctimespan dt = seconds(24*60*60);
    double temp = -10.0;
    double prec = 10.0;
    double radiation = 0.0;
    double wind_speed = 0.0;
    std::vector<std::pair<double, double>> tp(10, std::pair<double, double>(temp, prec));

    // Accumulate snow
    SkaugenModel model;
    response r;
    for_each(tp.begin(), tp.end(), [&dt , &p, &model, &radiation, &wind_speed, &s, &r] (std::pair<double, double> pair) {
            model.step(dt, p, pair.first, pair.second, radiation, wind_speed, s, r);
        });

    TS_ASSERT_DELTA(s.free_water, 0.0, 1.0e-6);  // No free water when dry snow precip
    model.step(dt, p, 10.0, 0.0, radiation, wind_speed, s, r);
    TS_ASSERT(s.free_water <= s.swe*max_water_fraction);  // Can not have more free water in the snow than the capacity of the snowpack
    tp = std::vector<std::pair<double, double>>(5, std::pair<double, double>(10.0, 0.0));
    for_each(tp.begin(), tp.end(), [&dt , &p, &model, &radiation, &wind_speed, &s, &r] (std::pair<double, double> pair) {
            model.step(dt, p, pair.first, pair.second, radiation, wind_speed, s, r);
        });
    TS_ASSERT_DELTA(s.free_water, s.swe*max_water_fraction, 1.0e-6);  // Test that snow is saturated with free water

}
TEST_CASE("skaugen_meltdown") {
#if 0
    Out[32]:

[('2017-07-17T18:00:00Z', 32.2),

('2017-07-17T19:00:00Z', 32.2),

('2017-07-17T20:00:00Z', 32.2),

('2017-07-17T21:00:00Z', 32.2),

('2017-07-17T22:00:00Z', 32.1),

('2017-07-17T23:00:00Z', 0.0),

('2017-07-18T00:00:00Z', 0.0),

('2017-07-18T01:00:00Z', 0.0),

('2017-07-18T02:00:00Z', 0.0),

('2017-07-18T03:00:00Z', 0.0)]



 State at the start of the interval where SWE drops to zero



 [{'swe': 32.1},

{'sca': 0.005033599471562574},

{'alpha': 0.127852277312898},

{'free_water': 3.2100000000000004},

{'nu': 0.012785227731289801},

{'residual': 0.0},

{'num_units': 321}]



 State at the end of the interval where SWE drops to zero



 [{'swe': 0.0},

{'sca': 0.0},

{'alpha': 40.05250550299988},

{'free_water': 0.0},

{'nu': 4.005250550299988},

{'residual': 0.0},

{'num_units': 0}]



 Input at the timestep where SWE drops to zero



[{'precipitation': 0.0010356738461072955},

{'temperature': 4.891358376624782},

{'radiation': 0.0042811137986584116},

{'rel_hum': 0.9769881287683342},

{'wind_speed': 6.411016839448434}]

#endif
    // Model parameters
    const double d_range = 113.0;
    const double unit_size = 0.1;
    const double alpha_0 = 40.77;
    const double max_water_fraction = 0.1;
    const double tx = 0.16;
    const double cx = 2.50;
    const double ts = 0.14;
    const double cfr = 0.01;
    parameter p(alpha_0, d_range, unit_size, max_water_fraction, tx, cx, ts, cfr);

    // Model state variables before melt-down
    const double alpha = 0.127852277312898;
    const double nu = 0.012785227731289801;
    const double sca = 0.005033599471562574;
    const double swe = 32.1;
    const double free_water = 3.21;
    const double residual = 0.0;
    const unsigned long nnn = 321;
    state s(nu, alpha, sca, swe, free_water, residual, nnn);

    // Model input
    utctimespan dt = seconds(3*3600);
    double temp = 4.891358376624782;
    double prec = 0.0010356738461072955;
    double radiation = 0.0042811137986584116;
    double wind_speed = 6.411016839448434;
    std::vector<std::pair<double, double>> tp(2, std::pair<double, double>(temp, prec));

    // melt it
    SkaugenModel model;
    response r;
    model.step(dt, p, temp, prec, radiation, wind_speed, s, r);
    TS_ASSERT(true);

}
}
