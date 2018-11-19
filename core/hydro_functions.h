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
	using std::vector;
	using std::cos;
	using std::sin;
	using std::pow;
	using std::exp;
	using namespace std;
	const double pi = boost::math::constants::pi<double>();

	/** \brief computes standard atmospheric pressure
     * \param height, [m] -- elevation of the point
     * \return p, [kPa] -- atmospheric pressure */
	double atm_pressure(double height){ // height < 11000
		const double p0 = 101325.0; //[Pa]standard sea level pressure
		const double L = 0.0065; //[K/m] temperature lapse rate
		const double g = 9.80665; //[m/s2] earth-surface gravitational acceleration
		const double R0 = 8.31447;//[J/mol/K] Universal gas constant
		const double M = 0.0289644; //[kg/mol] molar mass of dry air
		const double T0 = 288.15; // [K] sea level standard temperature
		return p0*pow((1 - L*height/T0),(g*M/R0/L))*Pa2kPa;
	}
	/**\brief computes actual vapor pressure from dewpoint temperature
     * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.113
     * \param temperature, [degC]
     * \param rhumidity, [percent] -- relative humidity
     * \return e, [kPa] -- actual vapor pressure*/
	double actual_vap_pressure(double temperature, double rhumidity){
		double es = (temperature>=0.0) ? (611 * exp(17.27*temperature/(temperature + 237.3))) : 611*exp(21.87*temperature/(temperature+265.5)) ; // saturation vapor pressure,[kPa], eq.(3.9)
		return rhumidity/100.0*es*Pa2kPa;//[kPa], based on eq.(3.12)
	}

}