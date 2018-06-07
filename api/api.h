/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#pragma once
#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <stdexcept>


/**
 * This file contains mostly typedefs and some few helper classes to
 * provide a python friendly header file to the api.python project.
 *
 *
 */

#include "core/utctime_utilities.h"
#include "core/time_axis.h"
#include "core/geo_point.h"
#include "core/geo_cell_data.h"
#include "core/time_series.h"
#include "core/region_model.h"
#include "core/model_calibration.h"
#include "core/bayesian_kriging.h"
#include "core/inverse_distance.h"

#include "core/pt_gs_k_cell_model.h"
#include "core/pt_hs_k_cell_model.h"
#include "core/pt_ss_k_cell_model.h"
#include "core/pt_hps_k_cell_model.h"
#include "core/hbv_stack_cell_model.h"

#include "core/time_series_dd.h"

namespace shyft {
  using namespace shyft::core;
  using namespace std;

  namespace api {
  using shyft::time_series::dd::apoint_ts;
  using shyft::time_series::dd::ats_vector;
  using shyft::core::stat_scope;
  
    /** \brief TsFactor provides time-series creation function using supplied primitives like vector of double, start, delta-t, n etc.
     */
    struct TsFactory {

        apoint_ts
        create_point_ts(int n, utctime tStart, utctimespan dt,
                        const std::vector<double>& values,
                        ts_point_fx interpretation=POINT_INSTANT_VALUE){
            return apoint_ts( time_axis::fixed_dt(tStart,dt, n), values, interpretation);
        }


        apoint_ts
        create_time_point_ts(utcperiod period, const std::vector<utctime>& times,
                             const std::vector<double>& values,
                             ts_point_fx interpretation=POINT_INSTANT_VALUE) {
            if (times.size() == values.size() + 1) {
                return apoint_ts( time_axis::point_dt(times), values, interpretation);
            } else if (times.size() == values.size()) {
                auto tx(times);
                tx.push_back(period.end > times.back()?period.end:times.back() + utctimespan(1));
                return apoint_ts( time_axis::point_dt(tx), values, interpretation);
            } else {
                throw std::runtime_error("create_time_point_ts times and values arrays must have corresponding count");
            }
        }
    };


    /** \brief GeoPointSource contains common properties, functions
     * for the geo-point located time-series.
     * Typically it contains a GeoPoint (3d position), plus a timeseries
     */
    struct GeoPointSource {
        GeoPointSource() =default;
        GeoPointSource(const geo_point& midpoint, const apoint_ts& ts)
          : mid_point_(midpoint), ts(ts) {}
        virtual ~GeoPointSource() {}
        typedef apoint_ts ts_t;
        typedef geo_point geo_point_t;

        geo_point mid_point_;
        apoint_ts ts;
		string uid;///< user-defined id, for external mapping,
        // boost python fixes for attributes and shared_ptr
        apoint_ts get_ts()  {return ts;}
        void set_ts(apoint_ts x) {ts=x;}
        bool is_equal(const GeoPointSource& x) const {return mid_point_==x.mid_point_ && ts == x.ts;}
        geo_point mid_point() const { return mid_point_; }
        bool operator==(const GeoPointSource&x) const {return is_equal(x);}
    };

    struct TemperatureSource : GeoPointSource {
        TemperatureSource()=default;
        TemperatureSource(const geo_point& p, const apoint_ts& ts)
         : GeoPointSource(p, ts) {}
        const apoint_ts& temperatures() const { return ts; }
		void set_temperature(size_t ix, double v) { ts.set(ix, v); } //btk dst compliant signature, used during btk-interpolation in gridpp exposure
		void set_value(size_t ix, double v) { ts.set(ix, v); }
        bool operator==(const TemperatureSource& x) const {return is_equal(x);}
    };

    struct PrecipitationSource : GeoPointSource {
        PrecipitationSource()=default;
        PrecipitationSource(const geo_point& p, const apoint_ts& ts)
         : GeoPointSource(p, ts) {}
        const apoint_ts& precipitations() const { return ts; }
		void set_value(size_t ix, double v) { ts.set(ix, v); }
		bool operator==(const PrecipitationSource& x) const {return is_equal(x);}
    };

    struct WindSpeedSource : GeoPointSource {
        WindSpeedSource()=default;
        WindSpeedSource(const geo_point& p, const apoint_ts& ts)
         : GeoPointSource(p, ts) {}
        bool operator==(const WindSpeedSource& x) const {return is_equal(x);}
    };

    struct RelHumSource : GeoPointSource {
        RelHumSource()=default;
        RelHumSource(const geo_point& p, const apoint_ts& ts)
         : GeoPointSource(p, ts) {}
        bool operator==(const RelHumSource& x) const {return is_equal(x);}
    };

    struct RadiationSource : GeoPointSource {
        RadiationSource()=default;
        RadiationSource(const geo_point& p, const apoint_ts& ts)
         : GeoPointSource(p, ts) {}
        bool operator==(const RadiationSource& x) const {return is_equal(x);}
    };

    struct a_region_environment {
            typedef PrecipitationSource precipitation_t;
            typedef TemperatureSource temperature_t;
            typedef RadiationSource radiation_t;
            typedef RelHumSource rel_hum_t;
            typedef WindSpeedSource wind_speed_t;
            /** make vectors non nullptr by  default */
            a_region_environment() {
                temperature=make_shared<vector<TemperatureSource>>();
                precipitation=make_shared<vector<PrecipitationSource>>();
                radiation=make_shared<vector<RadiationSource>>();
                rel_hum=make_shared<vector<RelHumSource>>();
                wind_speed=make_shared<vector<WindSpeedSource>>();
            }
            shared_ptr<vector<TemperatureSource>>   temperature;
            shared_ptr<vector<PrecipitationSource>> precipitation;
            shared_ptr<vector<RadiationSource>>     radiation;
            shared_ptr<vector<WindSpeedSource>>     wind_speed;
            shared_ptr<vector<RelHumSource>>        rel_hum;

            // our boost python needs these methods to get properties straight (most likely it can be fixed by other means but..)
            shared_ptr<vector<TemperatureSource>> get_temperature() {return temperature;}
            void set_temperature(shared_ptr<vector<TemperatureSource>> x) {temperature=x;}
            shared_ptr<vector<PrecipitationSource>> get_precipitation() {return precipitation;}
            void set_precipitation(shared_ptr<vector<PrecipitationSource>> x) {precipitation=x;}
            shared_ptr<vector<RadiationSource>> get_radiation() {return radiation;}
            void set_radiation(shared_ptr<vector<RadiationSource>> x) {radiation=x;}
            shared_ptr<vector<WindSpeedSource>> get_wind_speed() {return wind_speed;}
            void set_wind_speed(shared_ptr<vector<WindSpeedSource>> x) {wind_speed=x;}
            shared_ptr<vector<RelHumSource>> get_rel_hum() {return rel_hum;}
            void set_rel_hum(shared_ptr<vector<RelHumSource>> x) {rel_hum=x;}
    };

    typedef shyft::time_series::point_ts<time_axis::fixed_dt> result_ts_t;
    typedef std::shared_ptr<result_ts_t> result_ts_t_;
    using core::cell_statistics;
    /** \brief A class that facilitates fast state io, the yaml in Python is too slow
     *
     */

    template <typename cell>
    struct basic_cell_statistics {
        shared_ptr<vector<cell>> cells;
        explicit basic_cell_statistics( const shared_ptr<vector<cell>>& cells):cells(cells) {}

		double total_area(const vector<int>& indexes, stat_scope ix_type) const {
			double sum = 0.0;
			if (indexes.size() == 0) {
				for (const auto &c : *cells) sum += c.geo.area();
		    } else {
                cell_statistics::verify_cids_exist(*cells, indexes, ix_type);
				for (auto cid : indexes)
					for (size_t j=0;j<cells->size();++j) {
                        const auto&c = (*cells)[j];
						if (cell_statistics::is_match(c,ix_type,cid,j)) sum += c.geo.area();
                    }
			}
			return sum;
		}
		double forest_area(const vector<int>& indexes, stat_scope ix_type) const {
			double sum = 0.0;
			if (indexes.size() == 0) {
				for (const auto &c : *cells) sum += c.geo.area()*c.geo.land_type_fractions_info().forest();
			} else {
                cell_statistics::verify_cids_exist(*cells, indexes, ix_type);
				for (auto cid : indexes)
					for (const auto&c : *cells)
						if ( (int)c.geo.catchment_id() == cid) sum += c.geo.area()*c.geo.land_type_fractions_info().forest();
			}
			return sum;
		}
		double glacier_area(const vector<int>& indexes, stat_scope ix_type) const {
			double sum = 0.0;
			if (indexes.size() == 0) {
				for (const auto &c : *cells) sum += c.geo.area()*c.geo.land_type_fractions_info().glacier();
			} else {
                cell_statistics::verify_cids_exist(*cells, indexes, ix_type);
				for (auto cid : indexes)
					for (const auto&c : *cells)
						if ((int)c.geo.catchment_id() == cid) sum += c.geo.area()*c.geo.land_type_fractions_info().glacier();
			}
			return sum;
		}
		double lake_area(const vector<int>& indexes, stat_scope ix_type) const {
			double sum = 0.0;
			if (indexes.size() == 0) {
				for (const auto &c : *cells) sum += c.geo.area()*c.geo.land_type_fractions_info().lake();
			} else {
                cell_statistics::verify_cids_exist(*cells, indexes, ix_type);
				for (auto cid : indexes)
					for (const auto&c : *cells)
						if ((int)c.geo.catchment_id() == cid) sum += c.geo.area()*c.geo.land_type_fractions_info().lake();
			}
			return sum;
		}
		double reservoir_area(const vector<int>& indexes, stat_scope ix_type) const {
			double sum = 0.0;
			if (indexes.size() == 0) {
				for (const auto &c : *cells) sum += c.geo.area()*c.geo.land_type_fractions_info().reservoir();
			} else {
                cell_statistics::verify_cids_exist(*cells, indexes, ix_type);
				for (auto cid : indexes)
					for (const auto&c : *cells)
						if ((int)c.geo.catchment_id() == cid) sum += c.geo.area()*c.geo.land_type_fractions_info().reservoir();
			}
			return sum;
		}
		double unspecified_area(const vector<int>& indexes, stat_scope ix_type) const {
			double sum = 0.0;
			if (indexes.size() == 0) {
				for (const auto &c : *cells) sum += c.geo.area()*c.geo.land_type_fractions_info().unspecified();
			} else {
                cell_statistics::verify_cids_exist(*cells, indexes, ix_type);
				for (auto cid : indexes)
					for (const auto&c : *cells)
						if ((int)c.geo.catchment_id() == cid) sum += c.geo.area()*c.geo.land_type_fractions_info().unspecified();
			}
			return sum;
		}
        apoint_ts discharge(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                     sum_catchment_feature(*cells, indexes,
                            [](const cell& c){ return c.rc.avg_discharge; },ix_type));
        }
		vector<double> discharge(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.avg_discharge; }, ith_timestep,ix_type);
		}
		double discharge_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				sum_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.avg_discharge; }, ith_timestep,ix_type);
		}
        apoint_ts charge(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                    [](const cell& c) { return c.rc.charge_m3s; },ix_type));
        }
        vector<double> charge(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
            return cell_statistics::
                catchment_feature(*cells, indexes,
                    [](const cell& c) { return c.rc.charge_m3s; }, ith_timestep,ix_type);
        }
        double charge_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
            return cell_statistics::
                sum_catchment_feature_value(*cells, indexes,
                    [](const cell& c) { return c.rc.charge_m3s; }, ith_timestep,ix_type);
        }

		apoint_ts temperature(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                     average_catchment_feature(*cells, indexes,
                            [](const cell& c){ return c.env_ts.temperature; },ix_type));
        }
		vector<double> temperature(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.env_ts.temperature; }, ith_timestep,ix_type);
		}
		double temperature_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.env_ts.temperature; }, ith_timestep,ix_type);
		}

		apoint_ts precipitation(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                     average_catchment_feature(*cells, indexes,
                            [](const cell& c){ return c.env_ts.precipitation; },ix_type));
        }
		vector<double> precipitation(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.env_ts.precipitation; }, ith_timestep,ix_type);
		}
		double precipitation_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.env_ts.precipitation; }, ith_timestep,ix_type);
		}

		apoint_ts radiation(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                     average_catchment_feature(*cells, indexes,
                            [](const cell& c){ return c.env_ts.radiation; },ix_type));
        }
		vector<double> radiation(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.env_ts.radiation; }, ith_timestep,ix_type);
		}
		double radiation_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.env_ts.radiation; }, ith_timestep,ix_type);
		}

		apoint_ts wind_speed(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                     average_catchment_feature(*cells, indexes,
                            [](const cell& c){ return c.env_ts.wind_speed; },ix_type));
        }
		vector<double> wind_speed(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.env_ts.wind_speed; }, ith_timestep,ix_type);
		}
		double wind_speed_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.env_ts.wind_speed; }, ith_timestep,ix_type);
		}

		apoint_ts rel_hum(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                     average_catchment_feature(*cells, indexes,
                            [](const cell& c){ return c.env_ts.rel_hum; },ix_type));
        }
		vector<double> rel_hum(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.env_ts.rel_hum; }, ith_timestep,ix_type);
		}
		double rel_hum_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.env_ts.rel_hum; }, ith_timestep,ix_type);
		}
	};

    template <typename cell>
    struct kirchner_cell_state_statistics {
        shared_ptr<vector<cell>> cells;
        explicit kirchner_cell_state_statistics(const shared_ptr<vector<cell>>& cells) :cells(cells) {}

		apoint_ts discharge(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.kirchner_discharge; },ix_type));
        }
		vector<double> discharge(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.kirchner_discharge; }, ith_timestep,ix_type);
		}
		double discharge_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				sum_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.kirchner_discharge; }, ith_timestep,ix_type);
		}
	};

	template <typename cell>
	struct hbv_soil_cell_state_statistics {
		shared_ptr<vector<cell>> cells;
		explicit hbv_soil_cell_state_statistics(const shared_ptr<vector<cell>>& cells) :cells(cells) {}

		apoint_ts discharge(const vector<int>& indexes, stat_scope ix_type) const {
			return apoint_ts(*cell_statistics::
				sum_catchment_feature(*cells, indexes,
					[](const cell& c) { return c.sc.soil_moisture; },ix_type));
		}
		vector<double> discharge(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
					[](const cell& c) { return c.sc.soil_moisture; }, ith_timestep,ix_type);
		}
		double discharge_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				sum_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.soil_moisture; }, ith_timestep,ix_type);
		}
	};

	template <typename cell>											//To be checked & controlled
	struct hbv_tank_cell_state_statistics {
		shared_ptr<vector<cell>> cells;
		explicit hbv_tank_cell_state_statistics(const shared_ptr<vector<cell>>& cells) :cells(cells) {}

		apoint_ts discharge(const vector<int>& indexes, stat_scope ix_type) const {
			return apoint_ts(*cell_statistics::
				sum_catchment_feature(*cells, indexes,
					[](const cell& c) { return c.sc.tank_uz; },ix_type));			//to be modified
		}
		vector<double> discharge(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
					[](const cell& c) { return c.sc.tank_uz; }, ith_timestep,ix_type);		//to be modified
		}
		double discharge_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				sum_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.tank_uz; }, ith_timestep,ix_type);			//to be modified
		}
	};


    ///< cells with gamma_snow state collection gives access to time-series for state
    template <typename cell>
    struct gamma_snow_cell_state_statistics {
        shared_ptr<vector<cell>> cells;
        explicit gamma_snow_cell_state_statistics(const shared_ptr<vector<cell>>& cells) :cells(cells) {}

		apoint_ts albedo(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.gs_albedo; },ix_type));
        }
		vector<double> albedo(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.gs_albedo; }, ith_timestep,ix_type);
		}
		double albedo_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
				[](const cell& c) { return c.sc.gs_albedo; }, ith_timestep,ix_type);
		}

		apoint_ts lwc(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.gs_lwc; },ix_type));
        }
		vector<double> lwc(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.gs_lwc; }, ith_timestep,ix_type);
		}
		double lwc_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
				[](const cell& c) { return c.sc.gs_lwc; }, ith_timestep,ix_type);
		}

		apoint_ts surface_heat(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.gs_surface_heat; },ix_type));
        }
		vector<double> surface_heat(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.gs_surface_heat; }, ith_timestep,ix_type);
		}
		double surface_heat_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
				[](const cell& c) { return c.sc.gs_surface_heat; }, ith_timestep,ix_type);
		}

		apoint_ts alpha(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.gs_alpha; },ix_type));
        }
		vector<double> alpha(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.gs_alpha; }, ith_timestep,ix_type);
		}
		double alpha_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
				[](const cell& c) { return c.sc.gs_alpha; }, ith_timestep,ix_type);
		}

		apoint_ts sdc_melt_mean(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.gs_sdc_melt_mean; },ix_type));
        }
		vector<double> sdc_melt_mean(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.gs_sdc_melt_mean; }, ith_timestep,ix_type);
		}
		double sdc_melt_mean_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
				[](const cell& c) { return c.sc.gs_sdc_melt_mean; }, ith_timestep,ix_type);
		}

		apoint_ts acc_melt(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.gs_acc_melt; },ix_type));
        }
		vector<double> acc_melt(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.gs_acc_melt; }, ith_timestep,ix_type);
		}
		double acc_melt_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
				[](const cell& c) { return c.sc.gs_acc_melt; }, ith_timestep,ix_type);
		}

		apoint_ts iso_pot_energy(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.gs_iso_pot_energy; },ix_type));
        }
		vector<double> iso_pot_energy(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.gs_iso_pot_energy; }, ith_timestep,ix_type);
		}
		double iso_pot_energy_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
				[](const cell& c) { return c.sc.gs_iso_pot_energy; }, ith_timestep,ix_type);
		}

		apoint_ts temp_swe(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.gs_temp_swe; },ix_type));
        }
		vector<double> temp_swe(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.gs_temp_swe; }, ith_timestep,ix_type);
		}
		double temp_swe_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
				[](const cell& c) { return c.sc.gs_temp_swe; }, ith_timestep,ix_type);
		}
	};

    ///< cells with gamma_snow response give access to time-series for these
    template <typename cell>
    struct gamma_snow_cell_response_statistics {
        shared_ptr<vector<cell>> cells;
        explicit gamma_snow_cell_response_statistics(const shared_ptr<vector<cell>>& cells) :cells(cells) {}

		apoint_ts sca(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.snow_sca; },ix_type));
        }
		vector<double> sca(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.snow_sca; }, ith_timestep,ix_type);
		}
		double sca_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
				[](const cell& c) { return c.rc.snow_sca; }, ith_timestep,ix_type);
		}

		apoint_ts swe(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.snow_swe; },ix_type));
        }
		vector<double> swe(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.snow_swe; }, ith_timestep,ix_type);
		}
		double swe_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
				[](const cell& c) { return c.rc.snow_swe; }, ith_timestep,ix_type);
		}

		apoint_ts outflow(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.snow_outflow; },ix_type));
        }
		vector<double> outflow(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.snow_outflow; }, ith_timestep,ix_type);
		}
		double outflow_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				sum_catchment_feature_value(*cells, indexes,
				[](const cell& c) { return c.rc.snow_outflow; }, ith_timestep,ix_type);
		}
        apoint_ts glacier_melt(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                [](const cell& c) {return c.rc.glacier_melt;},ix_type)
            );
        }
        vector<double> glacier_melt(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
            return cell_statistics::
                catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.glacier_melt; }, ith_timestep,ix_type);
        }
        double glacier_melt_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
            return cell_statistics::
                sum_catchment_feature_value(*cells, indexes,
                [](const cell& c) { return c.rc.glacier_melt; }, ith_timestep,ix_type);
        }

	};

	///<-- universal snow section
    template <typename cell>
    struct universal_snow_cell_state_statistics {
        shared_ptr<vector<cell>> cells;
        universal_snow_cell_state_statistics(shared_ptr<vector<cell>> cells) :cells(cells) {}

		apoint_ts albedo(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.us_albedo; },ix_type));
        }
		vector<double> albedo(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.us_albedo; }, ith_timestep,ix_type);
		}
		double albedo_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.us_albedo; }, ith_timestep,ix_type);
		}

		apoint_ts lwc(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.us_lwc; },ix_type));
        }
		vector<double> lwc(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.us_lwc; }, ith_timestep,ix_type);
		}
		double lwc_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.us_lwc; }, ith_timestep,ix_type);
		}

		apoint_ts surface_heat(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.us_surface_heat; },ix_type));
        }
		vector<double> surface_heat(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.us_surface_heat; }, ith_timestep,ix_type);
		}
		double surface_heat_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.us_surface_heat; }, ith_timestep,ix_type);
		}

		apoint_ts alpha(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.us_alpha; },ix_type));
        }
		vector<double> alpha(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.us_alpha; }, ith_timestep,ix_type);
		}
		double alpha_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.us_alpha; }, ith_timestep,ix_type);
		}

		apoint_ts sdc_melt_mean(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.us_sdc_melt_mean; },ix_type));
        }
		vector<double> sdc_melt_mean(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.us_sdc_melt_mean; }, ith_timestep,ix_type);
		}
		double sdc_melt_mean_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.us_sdc_melt_mean; }, ith_timestep,ix_type);
		}

		apoint_ts acc_melt(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.us_acc_melt; },ix_type));
        }
		vector<double> acc_melt(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.us_acc_melt; }, ith_timestep,ix_type);
		}
		double acc_melt_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.us_acc_melt; }, ith_timestep,ix_type);
		}

		apoint_ts iso_pot_energy(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.us_iso_pot_energy; },ix_type));
        }
		vector<double> iso_pot_energy(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.us_iso_pot_energy; }, ith_timestep,ix_type);
		}
		double iso_pot_energy_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.us_iso_pot_energy; }, ith_timestep,ix_type);
		}

		apoint_ts temp_swe(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.us_temp_swe; },ix_type));
        }
		vector<double> temp_swe(const vector<int>& indexes,size_t ith_timestep,stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.us_temp_swe; }, ith_timestep,ix_type);
		}
		double temp_swe_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.us_temp_swe; }, ith_timestep,ix_type);
		}
	};

    ///< cells with universal_snow response give access to time-series for these
    template <typename cell>
    struct universal_snow_cell_response_statistics {
        shared_ptr<vector<cell>> cells;
        universal_snow_cell_response_statistics(shared_ptr<vector<cell>> cells) :cells(cells) {}

		apoint_ts sca(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.snow_sca; },ix_type));
        }
		vector<double> sca(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.snow_sca; }, ith_timestep,ix_type);
		}
		double sca_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.snow_sca; }, ith_timestep,ix_type);
		}

		apoint_ts swe(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.snow_swe; },ix_type));
        }
		vector<double> swe(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.snow_swe; }, ith_timestep,ix_type);
		}
		double swe_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.snow_swe; }, ith_timestep,ix_type);
		}

		apoint_ts outflow(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.snow_outflow; },ix_type));
        }
		vector<double> outflow(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.snow_outflow; }, ith_timestep,ix_type);
		}
		double outflow_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				sum_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.snow_outflow; }, ith_timestep,ix_type);
		}
        apoint_ts glacier_melt(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                    [](const cell& c) {
                return c.rc.glacier_melt;
            }
                )
            );
        }
        vector<double> glacier_melt(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
            return cell_statistics::
                catchment_feature(*cells, indexes,
                    [](const cell& c) { return c.rc.glacier_melt; }, ith_timestep,ix_type);
        }
        double glacier_melt_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
            return cell_statistics::
                sum_catchment_feature_value(*cells, indexes,
                    [](const cell& c) { return c.rc.glacier_melt; }, ith_timestep,ix_type);
        }

	};

	///<-- end universal snow section
    ///< access to skaugen's snow routine's state statistics
    template <typename cell>
    struct skaugen_cell_state_statistics {
        shared_ptr<vector<cell>> cells;
        explicit skaugen_cell_state_statistics(const shared_ptr<vector<cell>>& cells) :cells(cells) {}

		apoint_ts alpha(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.snow_alpha; },ix_type));
        }
		vector<double> alpha(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.snow_alpha; }, ith_timestep,ix_type);
		}
		double alpha_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.snow_alpha; }, ith_timestep,ix_type);
		}

		apoint_ts nu(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.snow_nu; },ix_type));
        }
		vector<double> nu(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.snow_nu; }, ith_timestep,ix_type);
		}
		double nu_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.snow_nu; }, ith_timestep,ix_type);
		}

		apoint_ts lwc(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.snow_lwc; },ix_type));
        }
		vector<double> lwc(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.snow_lwc; }, ith_timestep,ix_type);
		}
		double lwc_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.snow_lwc; }, ith_timestep,ix_type);
		}

		apoint_ts residual(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.snow_residual; },ix_type));
        }
		vector<double> residual(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.snow_residual; }, ith_timestep,ix_type);
		}
		double residual_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.snow_residual; }, ith_timestep,ix_type);
		}

		apoint_ts swe(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.snow_swe; },ix_type));
        }
		vector<double> swe(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.snow_swe; }, ith_timestep,ix_type);
		}
		double swe_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.snow_swe; }, ith_timestep,ix_type);
		}

		apoint_ts sca(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.snow_sca; },ix_type));
        }
		vector<double> sca(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.snow_sca; }, ith_timestep,ix_type);
		}
		double sca_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.snow_sca; }, ith_timestep,ix_type);
		}
	};

    ///< access to skaugen's snow routine response statistics
    template <typename cell>
    struct skaugen_cell_response_statistics {
        shared_ptr<vector<cell>> cells;
        explicit skaugen_cell_response_statistics(const shared_ptr<vector<cell>>& cells) : cells(cells) {}

		apoint_ts outflow(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.snow_outflow; },ix_type));
        }
		vector<double> outflow(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.snow_outflow; }, ith_timestep,ix_type);
		}
		double outflow_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				sum_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.snow_outflow; }, ith_timestep,ix_type);
		}

		apoint_ts total_stored_water(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.snow_total_stored_water; },ix_type));
        }
		vector<double> total_stored_water(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.snow_total_stored_water; }, ith_timestep,ix_type);
		}
		double total_stored_water_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.snow_total_stored_water; }, ith_timestep,ix_type);
		}
        apoint_ts glacier_melt(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                    [](const cell& c) {
                return c.rc.glacier_melt;
            }
                )
            );
        }
        vector<double> glacier_melt(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
            return cell_statistics::
                catchment_feature(*cells, indexes,
                    [](const cell& c) { return c.rc.glacier_melt; }, ith_timestep,ix_type);
        }
        double glacier_melt_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
            return cell_statistics::
                sum_catchment_feature_value(*cells, indexes,
                    [](const cell& c) { return c.rc.glacier_melt; }, ith_timestep,ix_type);
        }

	};

    ///< access to hbv snow routine's state statistics
    template <typename cell>
    struct hbv_snow_cell_state_statistics {
        shared_ptr<vector<cell>> cells;
        explicit hbv_snow_cell_state_statistics(const shared_ptr<vector<cell>>& cells) : cells(move(cells)) {}

		apoint_ts swe(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.snow_swe; },ix_type));
        }
		vector<double> swe(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.snow_swe; }, ith_timestep,ix_type);
		}
		double swe_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.snow_swe; }, ith_timestep,ix_type);
		}

		apoint_ts sca(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.sc.snow_sca; },ix_type));
        }
		vector<double> sca(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.sc.snow_sca; }, ith_timestep,ix_type);
		}
		double sca_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.snow_sca; }, ith_timestep,ix_type);
		}
		size_t find_size(const vector<int>& indexes, stat_scope ix_type) const {
			if (indexes.size() == 0)
				return (*cells)[0].sc.sp.size();
			auto ci = indexes[0];
			for (const auto& c : *cells) {
				if (c.geo.catchment_id() == ci)
					return c.sc.sp.size();
			}
			return 0;
		}
		ats_vector sp(const vector<int>& indexes, stat_scope ix_type) const {
			ats_vector r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(apoint_ts(*cell_statistics::
					average_catchment_feature(*cells, indexes,
						[&i](const cell& c) { return c.sc.sp[i]; }, ix_type)));
			}
			return r;
		}
		vector<vector<double>> sp(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			vector<vector<double>> r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(cell_statistics::
					catchment_feature(*cells, indexes,
						[&i](const cell& c) { return c.sc.sp[i]; }, ith_timestep,ix_type));
			}
			return r;
		}
		vector<double> sp_value(const vector<int>& indexes,
			size_t ith_timestep) const {
			vector<double> r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(cell_statistics::
					average_catchment_feature_value(*cells, indexes,
						[&i](const cell& c) {return c.sc.sp[i]; },
						ith_timestep));
			}
			return r;
		}
		ats_vector sw(const vector<int>& indexes, stat_scope ix_type) const {
			ats_vector r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(apoint_ts(*cell_statistics::
					average_catchment_feature(*cells, indexes,
						[&i](const cell& c) { return c.sc.sw[i]; }, ix_type)));
			}
			return r;
		}
		vector<vector<double>> sw(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			vector<vector<double>> r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(cell_statistics::
					catchment_feature(*cells, indexes,
						[&i](const cell& c) { return c.sc.sw[i]; }, ith_timestep,ix_type));
			}
			return r;
		}
		vector<double> sw_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			vector<double> r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(cell_statistics::
					average_catchment_feature_value(*cells, indexes,
						[&i](const cell& c) {return c.sc.sw[i]; },
						ith_timestep, ix_type));
			}
			return r;
		}

	};

    ///< access to hbv snow routine response statistics
    template <typename cell>
    struct hbv_snow_cell_response_statistics {
        shared_ptr<vector<cell>> cells;
        explicit hbv_snow_cell_response_statistics(const shared_ptr<vector<cell>>& cells) : cells(cells) {}

		apoint_ts outflow(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.snow_outflow; },ix_type));
        }
		vector<double> outflow(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.snow_outflow; }, ith_timestep,ix_type);
		}
		double outflow_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				sum_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.snow_outflow; }, ith_timestep,ix_type);
		}

		apoint_ts glacier_melt(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                    [](const cell& c) {
                        return c.rc.glacier_melt;
                    }
                )
            );
        }
        vector<double> glacier_melt(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.glacier_melt; }, ith_timestep,ix_type);
		}
		double glacier_melt_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				sum_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.glacier_melt; }, ith_timestep,ix_type);
		}

	};


    ///< access to hbv physical snow routine's state statistics
    template <typename cell>
    struct hbv_physical_snow_cell_state_statistics {
        shared_ptr<vector<cell>> cells;
        hbv_physical_snow_cell_state_statistics(shared_ptr<vector<cell>> cells) : cells(move(cells)) {}

		apoint_ts swe(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                    [](const cell& c) { return c.sc.hps_swe; },ix_type));
        }
		vector<double> swe(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
                    [](const cell& c) { return c.sc.hps_swe; }, ith_timestep,ix_type);
		}
		double swe_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.hps_swe; }, ith_timestep,ix_type);
		}

		apoint_ts sca(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                    [](const cell& c) { return c.sc.hps_sca; },ix_type));
        }
		vector<double> sca(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
                    [](const cell& c) { return c.sc.hps_sca; }, ith_timestep,ix_type);
		}
		double sca_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.sc.hps_sca; }, ith_timestep,ix_type);
		}

        apoint_ts surface_heat(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [] (const cell& c) { return c.sc.hps_surface_heat; },ix_type));
        }

        vector<double> surface_heat(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
            return cell_statistics::
                catchment_feature(*cells, indexes,
                    [] (const cell& c) { return c.sc.hps_surface_heat; },
                    ith_timestep,ix_type);
        }

        double surface_heat_value(const vector<int>& indexes,
                                  size_t ith_timestep, stat_scope ix_type) const {
            return cell_statistics::
                average_catchment_feature_value(*cells, indexes,
                        [] (const cell& c) {return c.sc.hps_surface_heat; },
                        ith_timestep,ix_type);
        }
		size_t find_size(const vector<int>& indexes, stat_scope ix_type) const {
			if (indexes.size() == 0)
				return (*cells)[0].sc.sp.size();
			auto ci = indexes[0];
			for (size_t j=0;j<cells->size();++j) {
                const auto &c =*cells[j];
				if (cell_statistics::is_match(c,ix_type,ci,j))
					return c.sc.sp.size();
			}
			return 0;
		}
		ats_vector albedo(const vector<int>& indexes, stat_scope ix_type) const {
			ats_vector r;
			size_t size = find_size(indexes,ix_type);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(apoint_ts(*cell_statistics::
					average_catchment_feature(*cells, indexes,
						[&i](const cell& c) { return c.sc.albedo[i]; },ix_type)));
			}
			return r;
		}
		vector<vector<double>> albedo(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			vector<vector<double>> r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(cell_statistics::
					catchment_feature(*cells, indexes,
						[&i](const cell& c) { return c.sc.albedo[i]; }, ith_timestep,ix_type));
			}
			return r;
		}
		vector<double> albedo_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			vector<double> r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(cell_statistics::
					average_catchment_feature_value(*cells, indexes,
						[&i](const cell& c) {return c.sc.albedo[i]; },
						ith_timestep, ix_type));
			}
			return r;
		}
		ats_vector iso_pot_energy(const vector<int>& indexes, stat_scope ix_type) const {
			ats_vector r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(apoint_ts(*cell_statistics::
					average_catchment_feature(*cells, indexes,
						[&](const cell& c) { return c.sc.iso_pot_energy[i]; },ix_type)));
			}
			return r;
		}
		vector<vector<double>> iso_pot_energy(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			vector<vector<double>> r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(cell_statistics::
					catchment_feature(*cells, indexes,
						[&i](const cell& c) { return c.sc.iso_pot_energy[i]; }, ith_timestep,ix_type));
			}
			return r;
		}

		vector<double> iso_pot_energy_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			vector<double> r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(cell_statistics::
					average_catchment_feature_value(*cells, indexes,
						[&i](const cell& c) {return c.sc.iso_pot_energy[i]; },
						ith_timestep,ix_type));
			}
			return r;
		}
		ats_vector sp(const vector<int>& indexes, stat_scope ix_type) const {
			ats_vector r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(apoint_ts(*cell_statistics::
					average_catchment_feature(*cells, indexes,
						[&i](const cell& c) { return c.sc.sp[i]; },ix_type)));
			}
			return r;
		}
		vector<vector<double>> sp(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			vector<vector<double>> r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(cell_statistics::
					catchment_feature(*cells, indexes,
						[&i](const cell& c) { return c.sc.sp[i]; }, ith_timestep,ix_type));
			}
			return r;
		}
		vector<double> sp_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			vector<double> r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(cell_statistics::
					average_catchment_feature_value(*cells, indexes,
						[&i](const cell& c) {return c.sc.sp[i]; },
						ith_timestep, ix_type));
			}
			return r;
		}
		ats_vector sw(const vector<int>& indexes, stat_scope ix_type) const {
			ats_vector r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(apoint_ts(*cell_statistics::
					average_catchment_feature(*cells, indexes,
						[&i](const cell& c) { return c.sc.sw[i]; },ix_type)));
			}
			return r;
		}
		vector<vector<double>> sw(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			vector<vector<double>> r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(cell_statistics::
					catchment_feature(*cells, indexes,
						[&i](const cell& c) { return c.sc.sw[i]; }, ith_timestep,ix_type));
			}
			return r;
		}
		vector<double> sw_value(const vector<int>& indexes,	size_t ith_timestep, stat_scope ix_type) const {
			vector<double> r;
			size_t size = find_size(indexes);
			for (size_t i = 0; i < size; ++i) {
				r.push_back(cell_statistics::
					average_catchment_feature_value(*cells, indexes,
						[&i](const cell& c) {return c.sc.sw[i]; },
						ith_timestep,ix_type));
			}
			return r;
		}
	};


    ///< access to hbv physical snow routine response statistics
    template <typename cell>
    struct hbv_physical_snow_cell_response_statistics {
        shared_ptr<vector<cell>> cells;
        hbv_physical_snow_cell_response_statistics(shared_ptr<vector<cell>> cells) : cells(cells) {}

		apoint_ts outflow(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.hps_outflow; },ix_type));
        }
		vector<double> outflow(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.hps_outflow; }, ith_timestep,ix_type);
		}
		double outflow_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				sum_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.hps_outflow; }, ith_timestep,ix_type);
		}

		apoint_ts glacier_melt(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                sum_catchment_feature(*cells, indexes,
                    [](const cell& c) {
                        return c.rc.glacier_melt;
                    }
                )
            );
        }
        vector<double> glacier_melt(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.glacier_melt; }, ith_timestep,ix_type);
		}
		double glacier_melt_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				sum_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.glacier_melt; }, ith_timestep,ix_type);
		}

	};


    template <typename cell>
    struct priestley_taylor_cell_response_statistics {
        shared_ptr<vector<cell>> cells;
        explicit priestley_taylor_cell_response_statistics(const shared_ptr<vector<cell>>& cells) :cells(cells) {}

		apoint_ts output(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.pe_output; },ix_type));
        }
		vector<double> output(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.pe_output; }, ith_timestep,ix_type);
		}
		double output_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.pe_output; }, ith_timestep,ix_type);
		}
	};


	template <typename cell>
	struct hbv_soil_cell_response_statistics {
		shared_ptr<vector<cell>> cells;
		explicit hbv_soil_cell_response_statistics(const shared_ptr<vector<cell>>& cells) :cells(cells) {}

		apoint_ts output(const vector<int>& indexes, stat_scope ix_type) const {
			return apoint_ts(*cell_statistics::
				average_catchment_feature(*cells, indexes,
					[](const cell& c) { return c.rc.soil_outflow; },ix_type));
		}
		vector<double> output(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
					[](const cell& c) { return c.rc.soil_outflow; }, ith_timestep,ix_type);
		}
		double output_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.soil_outflow; }, ith_timestep,ix_type);
		}
	};

    /** helper class for pot_ratio method of ae cell_response */
    struct pot_ratio_single_value_ts {
        double v;
        double value(size_t) const {return v;};
    };
    template <typename cell>
    struct actual_evapotranspiration_cell_response_statistics {
        shared_ptr<vector<cell>> cells;
        explicit actual_evapotranspiration_cell_response_statistics(const shared_ptr<vector<cell>>& cells) :cells(cells) {}

		apoint_ts output(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) { return c.rc.ae_output; },ix_type));
        }
		vector<double> output(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[](const cell& c) { return c.rc.ae_output; }, ith_timestep,ix_type);
		}
		double output_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.ae_output; }, ith_timestep,ix_type);
		}

		/**  pot_ratio is needed to study the dry-level of kirchner and pot/act evapo ratio.
         *
         *  pot_ratio(t) = (1.0 - std::exp(-water_level*3.0/ scale_factor))  [ratio]
         *
         * where
         *   water_level = kirchner.state.q(t) [mm/h]
         *   scale_factor = cell.parameter.ae.scale_factor
         */
		apoint_ts pot_ratio(const vector<int>& indexes, stat_scope ix_type) const {
            return apoint_ts(*cell_statistics::
                average_catchment_feature(*cells, indexes,
                [](const cell& c) {
                    vector<double>  v; v.reserve(c.sc.kirchner_discharge.size());
                    auto scale_factor = c.parameter->ae.ae_scale_factor;
                    auto area_m2 = c.geo.area();
                    for(size_t i=0;i<c.sc.kirchner_discharge.size();++i) {
                        v.emplace_back(
                            shyft::core::actual_evapotranspiration::calc_pot_ratio(m3s_to_mmh(c.sc.kirchner_discharge.value(i),area_m2),scale_factor)
                        );
                    }
                    return shyft::core::pts_t(c.sc.kirchner_discharge.time_axis(),move(v),c.sc.kirchner_discharge.point_interpretation());
                },ix_type)
            );
        }
		vector<double> pot_ratio(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
				[ith_timestep](const cell& c) {
                    auto scale_factor = c.parameter->ae.ae_scale_factor;
                    auto area_m2 = c.geo.area();
                    auto water_level = m3s_to_mmh(c.sc.kirchner_discharge.value(ith_timestep),area_m2);
                    double pot_ratio_value=shyft::core::actual_evapotranspiration::calc_pot_ratio(water_level,scale_factor);
                    return pot_ratio_single_value_ts{pot_ratio_value};
                }, ith_timestep,ix_type);
		}
		double pot_ratio_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[ith_timestep](const cell& c) {
                        auto scale_factor = c.parameter->ae.ae_scale_factor;
                        auto area_m2 = c.geo.area();
                        auto water_level = m3s_to_mmh(c.sc.kirchner_discharge.value(ith_timestep),area_m2);
                        double pot_ratio_value=shyft::core::actual_evapotranspiration::calc_pot_ratio(water_level,scale_factor);
                        return pot_ratio_single_value_ts{pot_ratio_value};
                    }, ith_timestep,ix_type);
		}

    };

	template <typename cell>
	struct hbv_actual_evapotranspiration_cell_response_statistics {
		shared_ptr<vector<cell>> cells;
		explicit hbv_actual_evapotranspiration_cell_response_statistics(const shared_ptr<vector<cell>>& cells) :cells(cells) {}

		apoint_ts output(const vector<int>& indexes, stat_scope ix_type) const {
			return apoint_ts(*cell_statistics::
				average_catchment_feature(*cells, indexes,
					[](const cell& c) { return c.rc.ae_output; },ix_type));
		}
		vector<double> output(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				catchment_feature(*cells, indexes,
					[](const cell& c) { return c.rc.ae_output; }, ith_timestep,ix_type);
		}
		double output_value(const vector<int>& indexes, size_t ith_timestep, stat_scope ix_type) const {
			return cell_statistics::
				average_catchment_feature_value(*cells, indexes,
					[](const cell& c) { return c.rc.ae_output; }, ith_timestep,ix_type);
		}
	};

	/**\brief geo_cell_data_io provide fast conversion to-from string
	 *
	 * In the context of orchestration/repository, we found that it could be useful to
     * cache information as retrieved from a GIS system
     *
     * This only makes sense to keep as long as we need it for performance reasons
     *
	 */
	struct geo_cell_data_io {
        static size_t size() {return 11;}// number of doubles to store a gcd
        static void push_to_vector(vector<double>&v,const geo_cell_data& gcd ){
            v.push_back(gcd.mid_point().x);
            v.push_back(gcd.mid_point().y);
            v.push_back(gcd.mid_point().z);
            v.push_back(gcd.area());
            v.push_back(int(gcd.catchment_id()));
            v.push_back(gcd.radiation_slope_factor());
            v.push_back(gcd.land_type_fractions_info().glacier());
            v.push_back(gcd.land_type_fractions_info().lake());
            v.push_back(gcd.land_type_fractions_info().reservoir());
            v.push_back(gcd.land_type_fractions_info().forest());
            v.push_back(gcd.land_type_fractions_info().unspecified());// not really needed, since it can be computed from 1- theothers
        }
        static vector<double> to_vector(const geo_cell_data& gcd) {
            vector<double> v;v.reserve(11);
            push_to_vector(v,gcd);
            return v;
        }
        static geo_cell_data from_raw_vector(const double *v) {
            land_type_fractions ltf; ltf.set_fractions(v[6],v[7],v[8],v[9]);
            //                               x    y    z    a    cid       rsl
            return geo_cell_data(geo_point(v[0],v[1],v[2]),v[3],int(v[4]),v[5],ltf);
        }
        static geo_cell_data from_vector(const vector<double>&v) {
            if(v.size()!= size())
                throw invalid_argument("geo_cell_data_io::from_vector: size of vector must be equal to geo_cell_data_io::size()");
            return from_raw_vector(v.data());

        }


    };
  } // Namespace api
}// Namespace shyft
