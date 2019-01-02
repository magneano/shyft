#include "test_pch.h"

#include <cmath>
#include <limits>
#include <vector>

#include "core/utctime_utilities.h"
#include "core/time_axis.h"
#include "core/time_series_dd.h"

static inline shyft::core::utctime _t(int64_t t1970s) {return shyft::core::utctime{shyft::core::seconds(t1970s)};}


TEST_SUITE("time_series") {

    using shyft::core::no_utctime;
    using shyft::core::seconds;
    using std::numeric_limits;
    const double eps = numeric_limits<double>::epsilon();
    using shyft::time_series::dd::apoint_ts;
    using shyft::time_axis::generic_dt;
    using shyft::time_series::ts_point_fx;
    using std::vector;
    using std::make_shared;
    using std::isfinite;

    using shyft::time_series::dd::qac_parameter;
    using shyft::time_series::dd::qac_ts;

    TEST_CASE("qac_parameter") {

        qac_parameter q;

        SUBCASE("no limits set, allow all values, except nan") {
            FAST_CHECK_EQ(q.is_ok_min_max(shyft::nan),false);
            FAST_CHECK_EQ(q.is_ok_min_max(1.0),true);
        }

        SUBCASE("min/max abs limits") {
            q.max_x=1.0;
            FAST_CHECK_EQ(q.is_ok_min_max(1.0),true);
            FAST_CHECK_EQ(q.is_ok_min_max(1.0+eps),false);
            q.min_x=-1.0;
            FAST_CHECK_EQ(q.is_ok_min_max(-1.0),true);
            FAST_CHECK_EQ(q.is_ok_min_max(-1.0-eps),false);
        }
    }

    TEST_CASE("qac_ts") {
        generic_dt ta{_t(0),seconds(10),5};
        //                 0    1       2     3    4
        vector<double> v {0.0,1.0,shyft::nan,3.0,-20.1};
        vector<double>ev {0.0,1.0,       2.0,3.0,-20.1};
        vector<double>cv {1.0,0.0,      -1.0,3.0,-20.1};
        apoint_ts src(ta,v,ts_point_fx::POINT_AVERAGE_VALUE);
        apoint_ts cts{ta,cv,ts_point_fx::POINT_AVERAGE_VALUE};
        qac_parameter qp = qac_parameter::create_min_max_linear_interpolation_parameters(true, true, shyft::nan, shyft::nan, shyft::core::max_utctime);
        auto ts = make_shared<qac_ts>(src,qp);


        // verify simple min-max limit cases
        FAST_CHECK_UNARY(ts.get()!=nullptr);
        FAST_CHECK_EQ(ts->value(2),doctest::Approx(2.0));
        FAST_CHECK_EQ(ts->value_at(ts->time(2)),doctest::Approx(2.0));
        FAST_CHECK_EQ(ts->value_at(ts->time(2)+seconds(1)),doctest::Approx(2.0));
        src.set_point_interpretation(ts_point_fx::POINT_INSTANT_VALUE);
        FAST_CHECK_EQ(ts->value_at(ts->time(2)+ seconds(1) ),doctest::Approx(2.1));
        FAST_CHECK_EQ(ts->value_at(ts->time(2)- seconds(1) ),doctest::Approx(1.9));
        ts->p.min_x = 0.0;
        FAST_CHECK_UNARY(!isfinite(ts->value_at(ts->time(3)+ seconds(1) )));
        ts->p.min_x = -40.0;
        FAST_CHECK_UNARY(isfinite(ts->value_at(ts->time(3)+ seconds(1) )));
        ts->p.max_x = 2.9; // clip 3.0 out of range
        FAST_CHECK_EQ(ts->value(2),
                      doctest::Approx(
                         v[1] +10*(v[4]-v[1])/(40-10)
                      )
        );

        FAST_CHECK_EQ(ts->value(3),doctest::Approx(
                         v[1] +20*(v[4]-v[1])/(40-10)
                      )
        );

        ts->p.max_x = shyft::nan;
        auto qv=ts->values();
        for(size_t i=0;i<ts->size();++i)
            FAST_CHECK_EQ(qv[i], doctest::Approx(ev[i]));
        src.set(0,shyft::nan);
        src.set(1,shyft::nan);
        FAST_CHECK_UNARY(!isfinite(ts->value(0)));
        FAST_CHECK_UNARY(!isfinite(ts->value(1)));
        FAST_CHECK_UNARY(!isfinite(ts->value(2)));
        FAST_CHECK_UNARY(isfinite(ts->value(3)));
        FAST_CHECK_UNARY(isfinite(ts->value(4)));

        ts->cts=cts.ts;// now set in replacement values as time-series
        FAST_CHECK_EQ(ts->value(0),doctest::Approx(cts.value(0)));
        FAST_CHECK_EQ(ts->value(1),doctest::Approx(cts.value(1)));
        FAST_CHECK_EQ(ts->value(2),doctest::Approx(cts.value(2)));
        FAST_CHECK_EQ(ts->value(3),doctest::Approx(v[3]));// own value
        FAST_CHECK_EQ(ts->value(4),doctest::Approx(v[4]));// own value
        ts->p.min_x = 0.0;// clip out neg values
        FAST_CHECK_EQ(ts->value(4),doctest::Approx(cts.value(4)));// own value replaces -20.1
    }

    TEST_CASE("qac_ts_repeating") {

        SUBCASE("No repetitions") {
            /* Test that in a series without any repetitions no values are replaced.
             */

            std::vector<double> values{0., 1., 3., 2., 4., 2., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
            generic_dt ta{ _t(0), seconds(10), values.size() };

            std::vector<double> expected = values;
        
            apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

            qac_parameter params = qac_parameter::create_repeating_no_fill_parameters(
                false, false, seconds(60), 1e-2
            );
            auto ts = qac_ts(data_ts, params);
            auto result = ts.values();

            for ( std::size_t i = 0; i < values.size(); ++i ) {
                FAST_CHECK_EQ(result[i], expected[i]);
            }
        }

        SUBCASE("One repetition") {

            SUBCASE("Exact detection period span") {
                /* Test that a single repetition sequence is detected when the detection period is exactly as long as the repetition.
                 */

                //                                         *     *     *
                std::vector<double> values  {0., 1., 3.,   2.,   2.,   2., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
                std::vector<double> expected{0., 1., 3., 100., 100., 100., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params = qac_parameter::create_repeating_constant_fill_parameters(
                    false, false, seconds(30), 1e-2, 100.
                );
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }

            SUBCASE("Shorter detection period span") {
                /* Test that a single repetition sequence is detected when the detection period is shorter than the repetition sequence.
                 */

                //                                         *     *     *
                std::vector<double> values  {0., 1., 3.,   2.,   2.,   2., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
                std::vector<double> expected{0., 1., 3., 100., 100., 100., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params = qac_parameter::create_repeating_constant_fill_parameters(
                    false, false, seconds(20), 1e-2, 100.
                );
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }

            SUBCASE("Longer detection period span") {
                /* Test that a single repetition sequence is not detected when the detection period is longer than the repetition sequence.
                 */

                //                                       *   *   *
                std::vector<double> values  {0., 1., 3., 2., 2., 2., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
                std::vector<double> expected{0., 1., 3., 2., 2., 2., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params = qac_parameter::create_repeating_constant_fill_parameters(
                    false, false, seconds(40), 1e-2, 100.
                );
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }
        
            SUBCASE("All repetition") {
                /* Test that in a series with only repetition the entire series is flagged.
                 */

                std::vector<double> values  {  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.};
                std::vector<double> expected{100., 100., 100., 100., 100., 100., 100., 100., 100., 100., 100., 100.};
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params = qac_parameter::create_repeating_constant_fill_parameters(
                    false, false, seconds(20), 1e-2, 100.
                );
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }
        }

        SUBCASE("Adjacent repetitions") {
            /* Test that in a series with adjacent repeting seqences values are replaced as expected.
             */

            //                                         *     *     *     *     *     *     *
            std::vector<double> values  {0., 1., 3.,   2.,   2.,   2.,   3.,   3.,   3.,   3., 4., 2., 5., 2., 1. };
            std::vector<double> expected{0., 1., 3., 100., 100., 100., 100., 100., 100., 100., 4., 2., 5., 2., 1. };
            generic_dt ta{ _t(0), seconds(10), values.size() };

            apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

            qac_parameter params = qac_parameter::create_repeating_constant_fill_parameters(
                false, false, seconds(30), 1e-2, 100.
            );
            auto ts = qac_ts(data_ts, params);
            auto result = ts.values();

            for ( std::size_t i = 0; i < values.size(); ++i ) {
                FAST_CHECK_EQ(result[i], expected[i]);
            }
        }

        SUBCASE("Leading repetition") {
            /* Test that a leading repetition sequence is detected.
             */

            //                             *     *     *     *     *     *
            std::vector<double> values  {  1.,   1.,   1.,   1.,   1.,   1., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
            std::vector<double> expected{100., 100., 100., 100., 100., 100., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
            generic_dt ta{ _t(0), seconds(10), values.size() };

            apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

            qac_parameter params = qac_parameter::create_repeating_constant_fill_parameters(
                false, false, seconds(30), 1e-2, 100.
            );
            auto ts = qac_ts(data_ts, params);
            auto result = ts.values();

            for ( std::size_t i = 0; i < values.size(); ++i ) {
                FAST_CHECK_EQ(result[i], expected[i]);
            }
        }

        SUBCASE("End repetition") {
            /* Test that a repetition sequence at the end of the series is detected.
             */

            //                                                                 *     *     *     *     *     *
            std::vector<double> values  {3., 4., 5., 7., 4., 2., 5., 2., 3.,   1.,   1.,   1.,   1.,   1.,   1. };
            std::vector<double> expected{3., 4., 5., 7., 4., 2., 5., 2., 3., 100., 100., 100., 100., 100., 100. };
            generic_dt ta{ _t(0)
, seconds(10), values.size() };

            apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

            qac_parameter params = qac_parameter::create_repeating_constant_fill_parameters(
                false, false, seconds(30), 1e-2, 100.
            );
            auto ts = qac_ts(data_ts, params);
            auto result = ts.values();

            for ( std::size_t i = 0; i < values.size(); ++i ) {
                FAST_CHECK_EQ(result[i], expected[i]);
            }
        }

        SUBCASE("Repetition tolerance") {

            SUBCASE("Repetition around zero") {
                /* Test that repetitions around zero with a detection tolerance is replaced correctly.
                 */

                //                                            *     *     *     *     *
                std::vector<double> values  {1.1, 1.4, 2.5,  0.1, -0.2, -0.1,  0.0, -0.1, 2.5, 1.9, 2.1, 2.3, 3.1};
                std::vector<double> expected{1.1, 1.4, 2.5, 100., 100., 100., 100., 100., 2.5, 1.9, 2.1, 2.3, 3.1};
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params = qac_parameter::create_repeating_constant_fill_parameters(
                    false, false, seconds(40), 0.5, 100.
                );
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }

            SUBCASE("Repetition around value") {
                /* Test that repetitions around a value with a detection tolerance is replaced correctly.
                 */

                //                                            *     *     *     *     *
                std::vector<double> values  {1.1, 1.4, 2.5,  5.1,  4.9,  4.8,  5.0,  4.9, 2.5, 1.9, 2.1, 2.3, 3.1};
                std::vector<double> expected{1.1, 1.4, 2.5, 100., 100., 100., 100., 100., 2.5, 1.9, 2.1, 2.3, 3.1};
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params = qac_parameter::create_repeating_constant_fill_parameters(
                    false, false, seconds(40), 0.5, 100.
                );
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }

            SUBCASE("Repetition just inside tolerance") {
                /* Test that repetitions exactly within the tolerance threshold is replaced.
                 */

                const double eps = std::numeric_limits<double>::epsilon();

                std::vector<double> values  {1.1, 1.4, 2.5, 4.0 + 10*eps, 4.0 - 10*eps,  4.0 + 10*eps,  4.0 - 10*eps,  4.0 + 10*eps, 2.5, 1.9, 2.1, 2.3, 3.1};
                std::vector<double> expected{1.1, 1.4, 2.5,         100.,         100.,          100.,          100.,          100., 2.5, 1.9, 2.1, 2.3, 3.1};
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params = qac_parameter::create_repeating_constant_fill_parameters(
                    false, false, seconds(40), 20*eps, 100.
                );
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }

            SUBCASE("Repetition just outside tolerance") {
                /* Test that repetitions exactly outside the tolerance threshold is not replaced.
                */

                const double eps = std::numeric_limits<double>::epsilon();

                std::vector<double> values  {1.1, 1.4, 2.5, 4.0 + 10*eps, 4.0 - 10*eps,  4.0 + 10*eps,  4.0 - 10*eps,  4.0 + 10*eps, 2.5, 1.9, 2.1, 2.3, 3.1};
                std::vector<double> expected{1.1, 1.4, 2.5, 4.0 + 10*eps, 4.0 - 10*eps,  4.0 + 10*eps,  4.0 - 10*eps,  4.0 + 10*eps, 2.5, 1.9, 2.1, 2.3, 3.1};
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params = qac_parameter::create_repeating_constant_fill_parameters(
                    false, false, seconds(40), 15*eps, 100.
                );
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }
        }
    }
}
