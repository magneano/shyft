/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
///	Copyright 2012 Statkraft Energi A/S
///
///	This file is part of Shyft.
///
///	Shyft is free software: you can redistribute it and/or modify it under the terms of
/// the GNU Lesser General Public License as published by the Free Software Foundation,
/// either version 3 of the License, or (at your option) any later version.
///
///	Shyft is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
/// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
/// PURPOSE. See the GNU Lesser General Public License for more details.
///
///	You should have received a copy of the GNU Lesser General Public License along with
/// Shyft, usually located under the Shyft root directory in two files named COPYING.txt
/// and COPYING_LESSER.txt.	If not, see <http://www.gnu.org/licenses/>.
///
/// Adapted from early enki method programmed by Kolbjørn Engeland and Sjur Kolberg
///
#pragma once

#include "core/utctime_utilities.h"
#include "core/geo_cell_data.h"
#include <armadillo>
#include <vector>
#include <chrono>
#include <boost/math/constants/constants.hpp>
//#include <proj.h>

namespace shyft::core {
	using std::vector;
	using std::cos;
	using std::sin;
	using std::pow;
	using std::exp;
	using namespace std;
	const double pi = boost::math::constants::pi<double>();

	namespace radiation_dingman {

		struct parameter {
			double albedo = 0.1; // average albedo of the surrounding ground surface:0.15-0.25 -- grass, 0.10-0.15 -- coniferous forest, 0.15 - 0.25 -- deciduous forest, 0.04-0.08 -- open water, 0.15-0.35 -- bare soil
			double turbidity = 1.0; // 1.0 - clean, 0.5 -- dusty
			parameter(double albedo = 0.1, double turbidity = 1.0):albedo(albedo),turbidity(turbidity){}
		};

		//struct state {}; // No state variables for this method

		struct response {
			double rcs_radiation = 0.0;
			double slope_factor = 1.0;
			double db_radiation = 0.0;
			double dif_radiation = 0.0;
			double bs_radiation = 0.0;
		};

		template <class P, class R>
		struct calculator {
			P param;
			R response;
			double rcs_radiation_ = 0.0;
			double slope_factor_ = 0.0;
			explicit calculator(const P& p, R& r):param(p),response(r) {}
			void compute_ra_radiation(R& response, double latitude, double rhumidity, double temperature, utctime t,  double slope = 0.0, double aspect = 0.0, double elevation = 0.0){
				doy_ = utc.day_of_year(t);
				// std::cout << "doy: " << doy_ << std::endl;
				//lt_ = utc.calendar_units(t).hour + utc.calendar_units(t).minute/60.0;
				delta_ = compute_earth_declination(doy_);
				//std::cout << "declination: " << delta_*180/pi << std::endl;
				phi_ = latitude*pi/180;

				Ts_ = compute_rise_set(phi_,delta_);
				//std::cout << "time_sunset: " << Ts_ << std::endl;
				double KET = compute_total_extraterrestrial_radiation(phi_,delta_,Ts_,doy_);
//				std::cout << "ket: " << KET << std::endl;

				double phi_eq = equivalent_latitude(phi_,slope,aspect);
//				std::cout << "phi_eq: " << phi_eq*180/pi << std::endl;
				double diff = long_difference(phi_eq,slope,aspect);
//				std::cout << "diff: " << diff*180/pi << std::endl;
				double ts_sloped = compute_rise_set(phi_eq,delta_, diff);
//				std::cout << "time_sunset_slope: " << ts_sloped << std::endl;

				double KETs = compute_total_extraterrestrial_radiation(phi_eq,delta_,ts_sloped,doy_);
//				std::cout << "kets: " << KETs << std::endl;

				compute_transmissivities(rhumidity, temperature,Ts_, phi_, delta_, elevation);
				double Kdirh = direct_beam_radiation(phi_,delta_,Ts_,doy_);
//				std::cout << "kdirh: " << Kdirh << std::endl;

				double Kdirhs = direct_beam_radiation(phi_eq,delta_,ts_sloped,doy_);
//				std::cout << "kdirhs: " << Kdirhs << std::endl;

				double Kdif = diffuse_radiation(phi_,delta_, Ts_,doy_);
//				std::cout << "kdif: " << Kdif << std::endl;

				double Kbs = backscattered_radiation(phi_,delta_, Ts_, doy_);
//				std::cout << "kbs: " << Kbs << std::endl;
//                std::cout << "===============================" << std::endl;

				double rcs_horizontal = Kdirh+Kdif+Kbs;
				double rcs_slope = Kdirhs+Kdif+Kbs;
				response.rcs_radiation = rcs_slope;
				response.slope_factor = rcs_slope/rcs_horizontal;
				response.db_radiation = Kdirhs;
				response.dif_radiation = Kdif;
				response.bs_radiation = Kbs;
				return;

			}
		private:
			double phi_ = 0.0;
			double delta_ = 0.0;
			double tau_wa_ = 1.0; //(1-alpha_w), where alpha_w is a fraction of solar radiation absorbed by water vapor
			double tau_da_ = 1.0; //(1-alpha_d), where alpha_d is the fraction of solar radiation absorbed by dust and other solid aerosols
			double tau_ws_ = 1.0; //(1-rho_w), where rho_w is the fraction of solar radiation scattered by water vapor
			double tau_rs_ = 1.0; //(1-rho_rs), where rho_rs is the fraction of solar radiation scattered byair molecules (Rayleigh scattering, which is responsible for the blur of the sky)
			double tau_ds_ = 1.0; //(1-rho_d), where rho_d is the fraction of solar radiation scattered by dust and soil
			double total_tau_direct_ = 1.0;
			double total_tau_diffuse_ = 1.0;

			//const double gsc = 1367; // W/m2 -- solar constant
			const double gsc = 117.8; // MJ/m2*day -- solar constant
			//const double gsc = 4.910; // MJ/m2*hr -- solar constant
			const double omega = 0.2618;// angular velocity of the earth [rad/hr]
			const double deg2rad=pi/180; // degrees to radians
			const double rad2deg=180/pi; // rad to deg
			calendar utc;
			double doy_; // day of the yearI
			double Ts_;// time sunset

			double hour_angle(double lt){
				return 15*(lt-12)*deg2rad;
			}
			/**\brief computes earth declination angle
             * // ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, eq.(D.5)
             * \param doy -- day of the year
             * \return declination, [rad] */
			double  compute_earth_declination(double doy){
				double G = day_angle(doy);
				return  0.006918 - 0.399912*cos(G)+0.070257*sin(G)- 0.006758*cos(2*G)+0.000907*sin(2*G) - 0.002697*cos(3*G)+0.00148*sin(3*G);
			}
			/**\brief computes orbital eccentricity
             * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, eq.(D.4)*/
			double orbital_eccentricity(double doy){
				double G = day_angle(doy);
				return 1.000110 + 0.034221*cos(G) + 0.001280*sin(G)+0.000719*cos(2*G)+0.000077*sin(2*G);
			}
			/**\brief computes day angle
             * \param doy -- day of the year
             * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.572, eq.(D.3)*/
			double day_angle(double doy){
				return 2*pi/365.0*(doy-1);
			}
			/**\brief compute zenith angle
             * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, eq.(D.6)
             * \param latitude, [rad]
             * \param declination, [rad]
             * \param hour_angle, [rad]
             * \return zenoth angle, [rad]*/
			double zenith_angle(double latitude, double declination, double hour_angle){
				return acos(sin(latitude)*sin(declination)+cos(latitude)*cos(declination)*cos(hour_angle));
			}
			/**\brief computes extraterrestrial solar radiation
             * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.575, eq.(D.9)
             * \param latitude,[rad]
             * \param declination,[rad]
             * \param Ts,[hour] sunset time
             * return KET, [W/m^2]*/
			double compute_total_extraterrestrial_radiation(double latitude, double declination, double Ts, double doy){
				const double omega = 0.2618;// angular velocity of the earth [rad/hr]
				return 2*gsc*orbital_eccentricity(doy)*(sin(latitude)*sin(declination)*Ts+cos(latitude)*cos(declination)*sin(omega*Ts)/omega);
			}
			/**\brief compute time of sunset
             * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.575, eq.(D.7)
             * \param latitude
             * ]param declination
             * \return time_sunset*/
			double compute_rise_set(double latitude, double declination, double difference = 0.0){
				double time_sunset;
				if (abs(latitude - declination) >= pi / 2) { time_sunset = 0; }
				if (abs(latitude - declination) < pi / 2 and (abs(latitude + declination) >= pi/2)) { time_sunset = 12; }
				if (abs(latitude - declination) < pi/2 and (abs(latitude + declination) < pi/2)){time_sunset = (acos(-tan(declination)*tan(latitude))-difference)/omega;} // for horizontal surface
				return time_sunset;
			}
			/**\brief compute dayly average optical air mass
             * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.575, eq.(D.12)
             * \param elevation, [m]*/
			double optical_air_mass(double time_sunset, double latitude, double declination, double elevation){
				double A = 0.008307+sin(latitude)*sin(declination);
				double B = cos(latitude)*cos(declination);
				double C = -0.01259;
				double M = 0.0;
				if (time_sunset == 0.0){M = 39.7;}
				else if (B<A){
					M = 1/time_sunset*(1/omega*(1.021/pow((A*A-B*B),0.5))*acos((B+A*cos(omega*time_sunset))/(A+B*cos(omega*time_sunset)))+C*time_sunset);
				}
				else if (B>A){
					double log1 = (pow((B+A)*(1+cos(omega*time_sunset)),0.5));
					double log2 = (pow((B-A)*(1-cos(omega*time_sunset)),0.5));
					M =  1/time_sunset*(1/omega*(1.021/pow(-A*A+B*B,0.5))*log((log1+log2)/(log1-log2))+C*time_sunset);
				}//here is a mistake of B>A sqrt(A*A-B*B) will give complex value
				else {M = 1/time_sunset*(1/omega*1.021/A*tan(omega*time_sunset*0.5)+C);}
				return M*exp(-elevation/7000);
			}
			/**\brief computes transmissivities
             * \param rhumidity,[percent] relative humidity
             * \param temperature, [K]
             * \param time_sunset, [hour]
             * \param latitude, [rad]
             * \param declination, [rad]
             * \param elevation, [m] */
			void compute_transmissivities(double rhumidity, double temperature, double time_sunset, double latitude, double declination, double elevation = 0.0){
				double M = optical_air_mass(time_sunset,latitude,declination,elevation);
				std::cout << "M: " << M << std::endl;
				double W = 0.00493*(rhumidity/temperature)*exp(26.23-5416/temperature);
//				std::cout << "W: " << W << std::endl;
				tau_wa_  = 1 - 0.077 *pow(M*W,0.3);
//				std::cout << "tau_wa: " << tau_wa_ << std::endl;
				tau_da_ = pow(0.965,M);
//				std::cout << "tau_da_: " << tau_da_ << std::endl;
				tau_ws_ = 1 - 0.0225*M*W;
//				std::cout << "tau_ws_: " << tau_ws_ << std::endl;
				tau_rs_ = 0.972 - 0.08262*M+0.00933*M*M - 0.00095*M*M*M + 0.0000437*pow(M,4);
//				std::cout << "tau_rs: " << tau_rs_ << std::endl;
				tau_ds_ = pow(0.965,M);
				total_tau_direct_ = tau_wa_*tau_da_*tau_ws_*tau_rs_*tau_ds_;
				total_tau_diffuse_ = 0.5*tau_wa_*tau_da_*(1-tau_ws_*tau_rs_*tau_ds_);
			}

			double direct_beam_radiation(double latitude, double declination, double time_sunset, double doy){
				return compute_total_extraterrestrial_radiation(latitude, declination, time_sunset, doy)*total_tau_direct_;
			}
			double diffuse_radiation(double latitude, double declination, double time_sunset, double doy){
				return compute_total_extraterrestrial_radiation(latitude, declination, time_sunset, doy)*total_tau_diffuse_;
			}
			double backscattered_radiation(double latitude, double declination, double time_sunset, double doy){
				double kdirh = direct_beam_radiation(latitude,declination,time_sunset,doy);
				double kdif = diffuse_radiation(latitude, declination,time_sunset,doy);
				return param.albedo*(kdirh+kdif)*total_tau_diffuse_;
			}
			double equivalent_latitude(double latitude, double slope=0.0, double aspect=0.0){
				return asin(sin(slope)*cos(aspect)*cos(latitude)+cos(slope)*sin(latitude));
			}
			double long_difference(double latitude, double slope, double aspect){
				return atan(sin(slope)*sin(aspect)/cos(latitude)*cos(slope)-sin(latitude)*sin(slope)*cos(aspect));
			}
		};

/** \tparam C is a cell, like shyft-cell, */
		template <class C>
		vector<arma::vec> surface_normal( const vector<C>& cells){
			vector<arma::vec> r;
			for(const auto&c:cells) {
				double x=c.geo.mid_point().x;
				r.push_back(arma::vec({1.0,1.0,1.0}));
			}
			return r;
		}
	}

}