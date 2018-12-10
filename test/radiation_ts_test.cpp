#include "test_pch.h"
#include "core/utctime_utilities.h"
#include <armadillo>
#include <vector>
#include <chrono>
#include <boost/math/constants/constants.hpp>
#include "core/radiation.h"
#include <cmath>
#include <random>
#include <string>
#include <stdexcept>

#include "core/routing.h"
#include "core/mstack_param.h"

#include "mocks.h"
#include "core/time_series.h"
#include "core/utctime_utilities.h"
#include "core/cell_model.h"



//namespace shyft::core::radiation{
//    /** \tparam C is a cell, like shyft-cell, */
//    template <class C>
//    vector<arma::vec> surface_normal( const vector<C>& cells){
//        vector<arma::vec> r;
//        for(const auto&c:cells) {
//            double x=c.geo.mid_point().x;
//            r.push_back(arma::vec({1.0,1.0,1.0}));
//        }
//        return r;
//    }
//}
namespace shyft::core::radiation_model{
    using namespace shyft::core;
    using namespace std;



    struct radiation_ts  {



        struct parameter{
            typedef radiation::parameter rad_parameter_t;


            parameter(const rad_parameter_t& rad18
            )
                    : rad(rad) { /*Do nothing */}

            parameter()=default;
            parameter(const parameter&)=default;
            parameter(parameter&&)=default;
            parameter& operator=(const parameter &c)=default;
            parameter& operator=(parameter&&c)=default;

            rad_parameter_t rad;						// I followed pt_gs_k but pt_hs_k differ


            ///<calibration support, needs vector interface to params, size is the total count
            size_t size() const { return 22; }
            ///<calibration support, need to set values from ordered vector
            void set(const vector<double>& p) {
                if (p.size() != size())
                    throw runtime_error("HBV_Stack Parameter Accessor: .set size missmatch");
                int i = 0;
                rad.albedo = p[i++];
                rad.turbidity = p[i++];
            }

            ///< calibration support, get the value of i'th parameter
            double get(size_t i) const {
                switch (i) {
                    case 0: return rad.albedo;
                    case 1: return rad.turbidity;
                    default:
                        throw runtime_error("HBV_stack Parameter Accessor:.get(i) Out of range.");
                }
                return 0;
            }

            ///< calibration and python support, get the i'th parameter name
            string get_name(size_t i) const {
                static const char *names[] = {
                        "rad.albedo",
                        "rad.turbidity"
                };
                if (i >= size())
                    throw runtime_error("test_stack Parameter Accessor:.get_name(i) Out of range.");
                return names[i];
            }
        };

        struct response {
            // Model responses
            typedef radiation::response  rad_response_t;

            rad_response_t rad;
;
        };

        template<template <typename, typename> class A, class R, class T_TS, class RH_TS, class T,
                 class GCD, class P,class RC>
                        void run_radiation_model(const GCD& geo_cell_data,
                                const P& parameter,
                                                 const T& time_axis, int start_step, int  n_steps,
                                                 const T_TS& temp,
                                                 const RH_TS& rel_hum,
                                                 RC& response_collector) {
            // Access time series input data through accessors of template A (typically a direct accessor).
            using temp_accessor_t = A<T_TS, T>;
            using rel_hum_accessor_t = A<RH_TS, T>;
            auto temp_accessor = temp_accessor_t(temp, time_axis);
            auto rel_hum_accessor = rel_hum_accessor_t(rel_hum, time_axis);

            radiation::calculator < typename P::rad_parameter_t, typename R::rad_response_t> sw_radiation(parameter.rad);
            R response;

            size_t i_begin = n_steps > 0 ? start_step : 0;
            size_t i_end = n_steps > 0 ? start_step + n_steps : time_axis.size();
            for (size_t i = i_begin; i < i_end; ++i) {
                utcperiod period = time_axis.period(i);
                utctime t = time_axis.time(i);
                double temp = temp_accessor.value(i);
                double rel_hum = rel_hum_accessor.value(i);
                arma::vec surface_normal({0.0,0.0,1.0});

                sw_radiation.psw_radiation(response.rad, geo_cell_data.mid_point().x, t, surface_normal, temp, rel_hum, geo_cell_data.mid_point().z);//

                //std::cout<<response.rad.psw_radiation<<"  ===== "<<std::endl;


                response_collector.collect(i, response);///< \note collect the response valid for the i'th period (current state is now at the end of period)

            }
            response_collector.set_end_response(response);
        }

        typedef parameter parameter_t;
        typedef response response_t;
        typedef shared_ptr<parameter_t> parameter_t_;
        typedef shared_ptr<response_t>  response_t_;
        struct all_response_collector {
            double destination_area;///< in [m^2]
            // these are the one that we collects from the response, to better understand the model::
            pts_t rad_output;///< potential evap mm/h

            response_t end_reponse;///<< end_response, at the end of collected

            explicit all_response_collector() {}
            all_response_collector(const timeaxis_t& time_axis)
                    : rad_output(time_axis, 0.0) {}

            /**\brief called before run to allocate space for results */
            void initialize(const timeaxis_t& time_axis,int start_step,int n_steps) {
                ts_init(rad_output, time_axis, start_step, n_steps, ts_point_fx::POINT_AVERAGE_VALUE);
            }

            /**\brief Call for each time step, to collect needed information from R
            *
            * The R in this case is the Response type defined for the pt_g_s_k stack
            * and in principle, we can pick out all needed information from this.
            * The values are put into the plain point time-series at the i'th position
            * corresponding to the i'th simulation step, .. on the time-axis.. that
            * again gives us the concrete period in time that this value applies to.
            *
            */
            void collect(size_t idx, const response_t& response) {
                rad_output.set(idx, response.rad.psw_radiation);
            }
            void set_end_response(const response_t& r) { end_reponse = r; }
        };

    };


}


// Some typedefs for clarity
using namespace shyft::core;
using namespace shyft::time_series;
using namespace shyft::core::radiation_model;

using namespace shyfttest::mock;
using namespace shyfttest;

namespace rad = shyft::core::radiation;


typedef TSPointTarget<ta::point_dt> catchment_t;

namespace shyfttest {
    namespace mock {
        // need specialization for hbv_stack_response_t above
        template<> template<>
        void ResponseCollector<ta::fixed_dt>::collect<radiation_ts::response>(size_t idx, const radiation_ts::response& response) {

        }
        template <> template <>
        void DischargeCollector<ta::fixed_dt>::collect<radiation_ts::response>(size_t idx, const radiation_ts::response& response) {
            // hbv_outflow is given in mm, so compute the totals
        }
    };
}; // End namespace shyfttest

namespace shyft::test {

    class trapezoidal_average {
    private:
        double area = 0.0;
        double f_a = 0.0;; // Left hand side of next integration subinterval
        double t_start = 0.0; // Start of integration period
        double t_a = 0.0; // Left hand side time of next integration subinterval
    public:
        explicit trapezoidal_average() {}

        /** \brief initialize must be called to reset states before being used during ode integration.
         */
        void initialize(double f0, double t_start) {
            this->f_a = f0;
            this->t_start = t_start;
            t_a = t_start;
            area = 0.0;
        }

        /** \brief Add contribution to average using a simple trapezoidal rule
         *
         * See: http://en.wikipedia.org/wiki/Numerical_integration
         */
        void add(double f, double t) {
            area += 0.5*(f_a + f)*(t - t_a);
            f_a = f;
            t_a = t;
        }

        double result() const { return area/(t_a - t_start); }
    };

}

TEST_SUITE("radiation_ts_model") {
    //using shyft::core::radiation::parameter;
    //using shyft::core::radiation::response;
    using shyft::core::radiation::calculator;
    //    using shyft::core::radiation::surface_normal;
    using shyft::core::calendar;
    using shyft::core::utctime;
    using shyft::test::trapezoidal_average;
    using shyft::core::radiation_model::radiation_ts;
    // test basics: creation, etc

    TEST_CASE("test_call_stack") {
        xpts_t temp;
        xpts_t prec;
        xpts_t rel_hum;
        xpts_t wind_speed;
        xpts_t radiation;

        calendar cal;
        utctime t0 = cal.time(YMDhms(2014, 8, 1, 0, 0, 0));
        size_t n_ts_points = 3 * 24;
        utctimespan dt = deltahours(1);
        utctime t1 = t0 + n_ts_points*dt;
        shyfttest::create_time_series(temp, prec, rel_hum, wind_speed, radiation, t0, dt, n_ts_points);

        auto model_dt = deltahours(24);
        vector<utctime> times;
        for (utctime i = t0; i <= t1; i += model_dt)
        times.emplace_back(i);
        ta::fixed_dt time_axis(t0, dt, n_ts_points);
        ta::fixed_dt state_time_axis(t0, dt, n_ts_points + 1);
        // Initialize parameters
        rad::parameter rad_param;

        // should work automagically snow_state.distribute(snow_param);

        // Initialize response
        //response run_response;

        // Initialize collectors
        //shyft::core::radiation_model::all_response_collector response_collector(1000 * 1000, time_axis);

        //radiation::parameter parameter(rad_param);
        geo_cell_data geo_cell_data;
        radiation_ts::parameter parameter(rad_param);
        radiation_ts rm;
        radiation_ts::all_response_collector response_collector(time_axis);

        rm.run_radiation_model<direct_accessor, radiation_ts::response>(geo_cell_data, parameter, time_axis,0,0, temp,rel_hum,response_collector);

//        std::cout<<response.rad.psw_radiation<<"  ===== "<<std::endl;
        auto swrad = response_collector.rad_output;
        for (size_t i = 0; i < swrad.size(); ++i) {
            TS_ASSERT(std::isfinite(swrad.get(i).v) && swrad.get(i).v >= 0);
            std::cout<<swrad.get(i).v<<std::endl;
        }
    }


}
