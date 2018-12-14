/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "priestley_taylor.h"
#include "kirchner.h"
#include "hbv_physical_snow.h"
#include "glacier_melt.h"
#include "actual_evapotranspiration.h"
#include "precipitation_correction.h"
#include "unit_conversion.h"
#include "routing.h"
#include "mstack_param.h"

namespace shyft {
  namespace core {
    namespace pt_hps_k {
        using namespace std;

        struct parameter {
            typedef priestley_taylor::parameter pt_parameter_t;
            typedef hbv_physical_snow::parameter hps_parameter_t;
            typedef actual_evapotranspiration::parameter ae_parameter_t;
            typedef kirchner::parameter kirchner_parameter_t;
            typedef precipitation_correction::parameter precipitation_correction_parameter_t;
            typedef glacier_melt::parameter glacier_parameter_t;
            typedef routing::uhg_parameter routing_parameter_t;
            typedef mstack_parameter mstack_parameter_t;

            pt_parameter_t pt;
            hps_parameter_t hps;
            ae_parameter_t ae;
            kirchner_parameter_t  kirchner;
            precipitation_correction_parameter_t p_corr;
            glacier_parameter_t gm;
            routing_parameter_t routing;
            mstack_parameter_t msp;

            parameter(const pt_parameter_t& pt,
                        const hps_parameter_t& hps,
                        const ae_parameter_t& ae,
                        const kirchner_parameter_t& kirchner,
                        const precipitation_correction_parameter_t& p_corr,
                        glacier_parameter_t gm = glacier_parameter_t(),
                        routing_parameter_t routing=routing_parameter_t(),
                        mstack_parameter_t msp=mstack_parameter_t()
                      ) // for backwards compatibility pass default glacier parameter
             : pt(pt), hps(hps), ae(ae), kirchner(kirchner), p_corr(p_corr), gm(gm),routing(routing),msp{msp} { /* Do nothing */ }

			parameter()=default;
			parameter(const parameter&)=default;
			parameter(parameter&&)=default;
			parameter& operator=(const parameter &c)=default;
			parameter& operator=(parameter&&c)=default;
            ///< Calibration support, size is the total number of calibration parameters
            size_t size() const { return 24; }

            void set(const vector<double>& p) {
                if (p.size() != size())
                    throw runtime_error("pt_ss_k parameter accessor: .set size missmatch");
                int i = 0;
                kirchner.c1 = p[i++];
                kirchner.c2 = p[i++];
                kirchner.c3 = p[i++];
                ae.ae_scale_factor = p[i++];
                hps.lw = p[i++];
                hps.tx = p[i++];
                hps.cfr = p[i++];
                hps.wind_scale = p[i++];
                hps.wind_const = p[i++];
                hps.surface_magnitude = p[i++];
                hps.max_albedo = p[i++];
                hps.min_albedo = p[i++];
                hps.fast_albedo_decay_rate=p[i++];
                hps.slow_albedo_decay_rate=p[i++];
                hps.snowfall_reset_depth=p[i++];
                hps.calculate_iso_pot_energy = std::fabs(p[i++])<0.0001 ?false:true;
                gm.dtf = p[i++];
                p_corr.scale_factor = p[i++];
				pt.albedo = p[i++];
				pt.alpha = p[i++];
                routing.velocity = p[i++];
                routing.alpha = p[i++];
                routing.beta  = p[i++];
                msp.reservoir_direct_response_fraction = p[i++];
            }
            //
            ///< calibration support, get the value of i'th parameter
            double get(size_t i) const {
                switch (i) {
                    case  0:return kirchner.c1;
                    case  1:return kirchner.c2;
                    case  2:return kirchner.c3;
                    case  3:return ae.ae_scale_factor;
                    case  4:return hps.lw;
                    case  5:return hps.tx;
                    case  6:return hps.cfr;
                    case  7:return hps.wind_scale;
                    case  8:return hps.wind_const;
                    case  9:return hps.surface_magnitude;
                    case 10:return hps.max_albedo;
                    case 11:return hps.min_albedo;
                    case 12:return hps.fast_albedo_decay_rate;
                    case 13:return hps.slow_albedo_decay_rate;
                    case 14:return hps.snowfall_reset_depth;
                    case 15:return hps.calculate_iso_pot_energy?1.0:0.0;
                    case 16:return gm.dtf;
                    case 17:return p_corr.scale_factor;
					case 18:return pt.albedo;
					case 19:return pt.alpha;
                    case 20:return routing.velocity;
                    case 21:return routing.alpha;
                    case 22:return routing.beta;
                    case 23:return msp.reservoir_direct_response_fraction;

                default:
                    throw runtime_error("pt_hp_k parameter accessor:.get(i) Out of range.");
                }
                return 0.0;
            }

            ///< calibration and python support, get the i'th parameter name
            string get_name(size_t i) const {
                static const char *names[] = {
                    "kirchner.c1",
                    "kirchner.c2",
                    "kirchner.c3",
                    "ae.ae_scale_factor",
                    "hps.lw",
                    "hps.tx",
                    "hps.cfr",
                    "hps.wind_scale",
                    "hps.wind_const",
                    "hps.surface_magnitude",
                    "hps.max_albedo",
                    "hps.min_albedo",
                    "hps.fast_albedo_decay_rate",
                    "hps.slow_albedo_decay_rate",
                    "hps.snowfall_reset_depth",
                    "hps.calculate_iso_pot_energy",
                    "gm.dtf",
                    "p_corr.scale_factor",
					"pt.albedo",
					"pt.alpha",
                    "routing.velocity",
                    "routing.alpha",
                    "routing.beta",
                    "msp.reservoir_direct_response_fraction"
				};
                if (i >= size())
                    throw runtime_error("pt_hps_k parameter accessor:.get_name(i) Out of range.");
                return names[i];
            }

        };


        struct state {
            typedef hbv_physical_snow::state hps_state_t;
            typedef kirchner::state kirchner_state_t;

            hps_state_t hps;
            kirchner_state_t kirchner;
            state() {}
            state(const hps_state_t& hps, const kirchner_state_t& kirchner)
             : hps(hps), kirchner(kirchner) { /* Do nothing */ }
            state(const state& state) : hps(state.hps), kirchner(state.kirchner) {}
            bool operator==(const state& x) const {return hps==x.hps && kirchner==x.kirchner;}
            /**adjust kirchner q with the  specified scale-factor
            * to support the process of tuning output of a step
            * to a specified observed/wanted average
            */
            void adjust_q(double scale_factor) {kirchner.adjust_q(scale_factor);}
            state scale_snow(const double& snow_storage_area_fraction) const {
                state c{*this};
                c.hps.swe *= snow_storage_area_fraction;
                return c;
            }
            x_serialize_decl();
        };

        struct response {
            typedef priestley_taylor::response pt_response_t;
            typedef hbv_physical_snow::response hps_response_t;
            typedef actual_evapotranspiration::response ae_response_t;
            typedef kirchner::response kirchner_response_t;
            pt_response_t pt;
            hps_response_t hps;
            ae_response_t ae;
            kirchner_response_t kirchner;
            double gm_melt_m3s;
            // Stack response
            double total_discharge;
            double charge_m3s;
                        // scale snow parts relative snow_storage_area 
            response scale_snow(const double& snow_storage_area_fraction) const {
                response c{*this};
                c.hps.storage *= snow_storage_area_fraction;
                c.hps.outflow *=snow_storage_area_fraction;
                // are there others that we should also scale, sca, is a still meaningful, unscaled ?
                return c;
            }

        };



        template<template <typename, typename> class A, class R, class T_TS, class P_TS, class WS_TS, class RH_TS, class RAD_TS, class T,
        class S, class GEOCELLDATA, class P, class SC, class RC >
        void run(const GEOCELLDATA& geo_cell_data,
            const P& parameter,
            const T& time_axis, int start_step, int  n_steps,
            const T_TS& temp,
            const P_TS& prec,
            const WS_TS& wind_speed,
            const RH_TS& rel_hum,
            const RAD_TS& rad,
            S& state,
            SC& state_collector,
            RC& response_collector
            ) {
            // Access time series input data through accessors of template A (typically a direct accessor).
            using temp_accessor_t = A < T_TS, T > ;
            using prec_accessor_t = A < P_TS, T > ;
            using wind_speed_accessor_t = A<WS_TS, T>;
            using rel_hum_accessor_t = A < RH_TS, T > ;
            using rad_accessor_t = A < RAD_TS, T > ;

            auto temp_accessor = temp_accessor_t(temp, time_axis);
            auto prec_accessor = prec_accessor_t(prec, time_axis);
            auto wind_speed_accessor = wind_speed_accessor_t(wind_speed, time_axis);
            auto rel_hum_accessor = rel_hum_accessor_t(rel_hum, time_axis);
            auto rad_accessor = rad_accessor_t(rad, time_axis);


            // Initialize the method stack
            precipitation_correction::calculator p_corr(parameter.p_corr.scale_factor);
            priestley_taylor::calculator pt(parameter.pt.albedo, parameter.pt.alpha);
            const hbv_physical_snow::calculator<typename P::hps_parameter_t, typename S::hps_state_t, typename R::hps_response_t> hbv_physical_snow(parameter.hps);
            kirchner::calculator<kirchner::trapezoidal_average, typename P::kirchner_parameter_t> kirchner(parameter.kirchner);
            // extra check for state.snow. that need to match up with dimension of the parameter.hps, to avoid core-dumps if user pass an 'empty' state
            // into the system(it dimensions from the parameters to get right)
            state.hps.distribute(parameter.hps,false);// fixup state to match parameter dimensions for snow-bins, but only if they are empty

            R response;
            const double glacier_fraction = geo_cell_data.land_type_fractions_info().glacier();
            const double gm_direct = parameter.gm.direct_response; //glacier melt directly out of cell
            const double gm_routed = 1-gm_direct; // glacier melt routed through kirchner
            const double snow_storage_fraction = geo_cell_data.land_type_fractions_info().snow_storage();// on this part, snow builds up, and melts.-> season time-response
            const double kirchner_routed_prec =  geo_cell_data.land_type_fractions_info().reservoir()*(1.0-parameter.msp.reservoir_direct_response_fraction) + geo_cell_data.land_type_fractions_info().lake();
            const double direct_response_fraction = glacier_fraction*gm_direct + geo_cell_data.land_type_fractions_info().reservoir()*parameter.msp.reservoir_direct_response_fraction;// only direct response on reservoirs
            const double kirchner_fraction = 1 - direct_response_fraction;
            const double cell_area_m2 = geo_cell_data.area();
            const double glacier_area_m2 = geo_cell_data.area()*glacier_fraction;

            size_t i_begin = n_steps > 0 ? start_step : 0;
            size_t i_end = n_steps > 0 ? start_step + n_steps : time_axis.size();
            for (size_t i = i_begin; i < i_end; ++i) {
                utcperiod period = time_axis.period(i);
                double temp = temp_accessor.value(i);
                double rad = rad_accessor.value(i);
                double rel_hum = rel_hum_accessor.value(i);
                double prec = p_corr.calc(prec_accessor.value(i));
                double wind_speed = wind_speed_accessor.value(i);
                auto state_snow_scaled=state.scale_snow(snow_storage_fraction);
                state_collector.collect(i, state_snow_scaled);///< \note collect the state at the beginning of each period (the end state is saved anyway)

                hbv_physical_snow.step(state.hps, response.hps, period.start, period.timespan(), temp, rad, prec, wind_speed, rel_hum); // outputs mm/h, interpreted as over the entire area

                response.gm_melt_m3s = glacier_melt::step(parameter.gm.dtf,temp,geo_cell_data.area()*state.hps.sca,glacier_area_m2);// m3/s, that is, how much flow from the snow free glacier parts

                response.pt.pot_evapotranspiration = pt.potential_evapotranspiration(temp, rad, rel_hum)*to_seconds(calendar::HOUR);// mm/s -> mm/h, interpreted as over the entire area(!)
                response.ae.ae = actual_evapotranspiration::calculate_step(state.kirchner.q, response.pt.pot_evapotranspiration,
                                    parameter.ae.ae_scale_factor,std::max(state.hps.sca,glacier_fraction),  // a evap only on non-snow/non-glac area
                                    period.timespan());
                double gm_mmh= shyft::m3s_to_mmh(response.gm_melt_m3s, cell_area_m2);
                kirchner.step(period.start, period.end, state.kirchner.q, response.kirchner.q_avg, response.hps.outflow*snow_storage_fraction+prec*kirchner_routed_prec + gm_routed*gm_mmh, response.ae.ae); //all units mm/h over 'same' area

                response.total_discharge =
                      std::max(0.0, prec - response.ae.ae)*direct_response_fraction // when it rains, remove ae. from direct response
                    + gm_direct*gm_mmh  // glacier melt direct response                    
                    + response.kirchner.q_avg*kirchner_fraction;
                response.charge_m3s =
                    + shyft::mmh_to_m3s(prec, cell_area_m2)
                    - shyft::mmh_to_m3s(response.ae.ae, cell_area_m2)
                    + response.gm_melt_m3s
                    - shyft::mmh_to_m3s(response.total_discharge, cell_area_m2);
                response.hps.hps_state=state_snow_scaled.hps;//< note/sih: we need snow in the response due to calibration

                // Possibly save the calculated values using the collector callbacks.
                response_collector.collect(i, response.scale_snow(snow_storage_fraction));
                if(i+1==i_end)
                    state_collector.collect(i+1, state.scale_snow(snow_storage_fraction));///< \note last iteration,collect the  final state as well.

            }
            response_collector.set_end_response(response.scale_snow(snow_storage_fraction));
        }
    }
  } // core
} // shyft
  //-- serialization support shyft
x_serialize_export_key(shyft::core::pt_hps_k::state);
