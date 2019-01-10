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
#include <boost/math/constants/constants.hpp>
//#include <proj.h>

namespace shyft::core {


	namespace hydro_functions {
		using namespace std;
		const double pi = boost::math::constants::pi<double>();
		const double sigma = 4.90*0.000000001; //Stefan-Boltzmann constant  MJ/m^2 K^4
        const double gsc = 1367; // W/m2 -- solar constant
		const double Pa2kPa = 0.001;
        const double deg2rad = pi / 180; // degrees to radians
        const double rad2deg = 180 / pi; // rad to deg
        const double MJm2d2Wm2 = 0.086400;// MJ/m^2/day to W/m^2
        const double MJm2h2Wm2 = 0.0036; // MJ/m^2/hour to W/m^2
        const double lambda = 2.45; // latent heat of vaporization, [MJ/kg]
        const double inverse_lambda = 0.408; // inverse lambda, [kg/MJ]
        const double rhow = 1000; // water density, [kg/m^3]

		/** \brief computes standard atmospheric pressure
         * \param height, [m] -- elevation of the point
         * \return p, [kPa] -- atmospheric pressure */
		inline double atm_pressure(double height, double T0 = 288.15){ // height < 11000
			const double p0 = 101325.0; //[Pa]standard sea level pressure
			const double L = 0.0065; //[K/m] temperature lapse rate
			const double g = 9.80665; //[m/s2] earth-surface gravitational acceleration
			const double R0 = 8.31447;//[J/mol/K] Universal gas constant
			const double M = 0.0289644; //[kg/mol] molar mass of dry air
			//const double T0 = 288.15; // [K] sea level standard temperature
			return p0*pow((1 - L*height/T0),(g*M/R0/L))*Pa2kPa;
		}
		/**\brief computes actual vapor pressure from dewpoint temperature
         * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.113
         * \param temperature, [degC]
         * \param rhumidity, [percent] -- relative humidity
         * \return e, [kPa] -- actual vapor pressure*/
		inline double actual_vp(double temperature, double rhumidity){
			double es = (temperature>=0.0) ? (6108 * exp(17.27*temperature/(temperature + 237.3))) : 6108*exp(21.87*temperature/(temperature+265.5)) ; // saturation vapor pressure,[kPa], eq.(3.9)
			return rhumidity/100.0*es*Pa2kPa;//[kPa], based on eq.(3.12)
		}

		// the whole section is based on supplementary materials of
		// McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and McVicar, T. R.:
		// Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis,
		// Hydrol. Earth Syst. Sci., 17, 1331-1363, https://doi.org/10.5194/hess-17-1331-2013, 2013.
		// TODO: check units
		/**\brief wet-bulb temperature  McJannet et al., 2008b, Equation 2*/
		inline double wetbulb_temperature(double air_temp, double dew_temp, double actual_vp) {
			double tt = (dew_temp+237.3)*(dew_temp+237.3);
			double a = 0.00066*100*air_temp + 4098*actual_vp*dew_temp/tt;
			double b = 0.00066*100 + 4098*actual_vp/tt;
			return a/b;
		}
		/**\brief dew point temperature    McJannet et al., 2008b, Equation 26*/
		inline double dew_temperature(double act_vp) {
			return 116.9+237.3*log(act_vp)/(16.78 - log(act_vp));
		}
		/**\briet slope of the saturation vapor pressure curve  Allen et al., 1998, Equation 13 */
		inline double svp_slope(double air_temp) {
			return 4098*(0.6108*exp(17.27*air_temp/(237.3+air_temp)))/(air_temp+237.3)/(air_temp+237.3);
		}
		/**\brief saturation vapor pressure at temperature  Allen et al., 1998, Equation 11*/
		inline double svp(double temp){
			return 0.6108*exp(17.27*temp/(temp+237.3));
		}
		/**\brief daily saturation vapor pressure  Allen et al., 1998, Equation 12*/
		inline double svp_daily(double tmax, double tmin) {
			return 0.5*(svp(tmax)+svp(tmin));
		}
		/**\brief mean daily actual vapor pressure  Allen et al., 1998, Equation 17*/
		inline double avp_daily_mean(double tmax, double tmin, double rhmax, double rhmin) {
			return 0.5*(svp(tmax)*rhmax/100+svp(tmin)*rhmin/100);
		}
		/**\brief mean daily actual vapor pressure using dew point temperature*/
		inline double avp_daily_mean(double dew_temp) {
			return svp(dew_temp);
		}
		/**\brief psychrometric constant Allen et al., 1998, Equation 8*/
		inline double gamma(double pressure, double lambda = 2.45){
		    return 0.00163*pressure/lambda;
		}

	}


}