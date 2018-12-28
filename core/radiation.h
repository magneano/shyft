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
#include "core/hydro_functions.h"
//#include <proj.h>
namespace shyft {
//    namespace utility {
//        namespace geometry{
//            template <class Point>
//            struct triangle{
//                Point p1, p2, p3;
//                explicit triangle(const Point& p1, const Point& p2, const Point& p3):p1(p1),p2(p2), p3(p3) {
//                    make_triangle();
//                }
//
//                Point normal = std::make_tuple(0.0,0.0,0.0);
//                Point mid_point = std::make_tuple(0.0,0.0,0.0);
//                double slope{0.0};
//                double aspect{0.0};
//                double elevation{0.0};
//
//            private:
//                void make_triangle(){
//                    // find mid-point:
//                    double x1 = std::get<0>(p1);
//                    double x2 = std::get<0>(p2);
//                    double x3 = std::get<0>(p3);
//                    double y1 = std::get<1>(p1);
//                    double y2 = std::get<1>(p2);
//                    double y3 = std::get<1>(p3);
//                    double z1 = std::get<2>(p1);
//                    double z2 = std::get<2>(p2);
//                    double z3 = std::get<2>(p3);
//                    double x = (x1+x2+x3)/3.0;
//                    double y = (y1+y2+y3)/3.0;
//                    double z = (z1+z2+z3)/3.0;
//                    mid_point = std::make_tuple(x,y,z);
//
//                    //find normal vector
//                    Point p1p2 = std::make_tuple(x2-x1,y2-y1,z2-z1);
//                    Point p1p3 = std::make_tuple(x3-x1,y3-y1,z3-z1);
//                    double normal_x = std::get<1>(p1p2)*std::get<2>(p1p3) - std::get<2>(p1p2)*std::get<1>(p1p3);
//                    double normal_y = -std::get<0>(p1p2)*std::get<2>(p1p3) - std::get<2>(p1p2)*std::get<0>(p1p3);
//                    double normal_z = std::get<0>(p1p2)*std::get<1>(p1p3) - std::get<1>(p1p2)*std::get<0>(p1p3);
//                    normal = std::make_tuple(normal_x,normal_y,normal_z);
//
//                    // slope, aspect, elevation
//                    slope = atan2(pow(pow(normal_x,2)+pow(normal_y,2),0.5),normal_z);
//                    aspect = atan2(normal_y, normal_x);
//                    elevation = std::get<2>(mid_point);
//                }
//
//
//            };
//        }
//    }

    namespace core {

        namespace radiation {
            using namespace std;
            const double pi = boost::math::constants::pi<double>();
            using namespace hydro_functions;

            struct parameter {
                double albedo = 0.1; // average albedo of the surrounding ground surface:0.15-0.25 -- grass, 0.10-0.15 -- coniferous forest, 0.15 - 0.25 -- deciduous forest, 0.04-0.08 -- open water, 0.15-0.35 -- bare soil
                double turbidity = 1.0; // 1.0 - clean, 0.5 -- dusty
                parameter(double albedo = 0.1, double turbidity = 1.0) : albedo(albedo), turbidity(turbidity) {}
            };

            //struct state {}; // No state variables for this method

            struct response {
                double dir_radiation = 0.0; // direct beam
                double dif_radiation = 0.0; // diffuse
                double ref_radiation = 0.0; // reflected
                double psw_radiation = 0.0; // predicted clear sky solar radiation for inclined surface [W/m2]
                double tsw_radiation = 0.0; // translated  solar radiation on a sloping surface based on measured horizontal radiation [W/m^2]
                double lw_radiation = 0.0; // long-wave radiation ///TODO
                double net_radiation = 0.0; // net radiation /// TODO
            };

            template<class P, class R>
            struct calculator {
                P param;

                explicit calculator(const P &p) : param(p) {}

                double slope() const {
                    return slope_ * rad2deg;
                } // horizontal slope, [deg], should be calculated from normal vector
                double aspect() const {
                    return aspect_ * rad2deg;
                } // surface aspect angle, [deg], 0 -- due south, -pi/2 -- due east, pi/2 -- due west, +-pi -- due north

                double latitude() const {
                    return phi_ * rad2deg;
                }// latitude, [deg] should be available from cell?/// TODO: add PROJ4 for conversion from cartesian to wgs84

                double
                ra_radiation() const { return ra_; } // extraterrestrial solar radiation for inclined surface[W/m2]
                double ra_radiation_hor() { return rahor_; } // extraterrestrial solar radiation for horizontal surfaces

                double sun_rise() const { return omega1_24_ * rad2deg; }

                double sun_set() const { return omega2_24_ * rad2deg; }

                /** \brief computes instantaneous predicted short-wave clear-sky radiation (direct, diffuse, reflected) for inclined surfaces
                 * ref.: Allen, R. G.; Trezza, R. & Tasumi, M. Analytical integrated functions for daily solar radiation on slopes Agricultural and Forest Meteorology, 2006, 139, 55-73
                 * \param latitude, [deg]
                 * \param utctime,
                 * \param surface_normal
                 * \param temperature, [degC]
                 * \param rhumidity, [percent]
                 * \param elevation
                 * \return */
                template<class V>
                void psw_radiation(R &response, double latitude, utctime t, const V &surface_normal,
                                   double temperature = 0.0, double rhumidity = 40.0, double elevation = 0.0) {
                    doy_ = utc.day_of_year(t);
                    lt_ = utc.calendar_units(t).hour + utc.calendar_units(t).minute / 60.0;
                    delta_ = compute_earth_declination(doy_);
                    omega_ = hour_angle(lt_); // earth hour angle

                    /// TODO: consider implementation of algorithm from Dozier, J. & Frew, J. Rapid calculation of terrain parameters for radiation modeling from digital elevation data IEEE Transactions on Geoscience and Remote Sensing, 1990
                    arma::vec normal0{0.0, 0.0, 1.0}; // normal to horizontal plane
                    slope_ = acos(arma::norm_dot(surface_normal, normal0));
                    //slope_ = atan2(pow(pow(surface_normal(0),2)+pow(surface_normal(1),2),0.5),surface_normal(2));
                    //aspect_ = pi - atan2(surface_normal(1),surface_normal(0));
                    aspect_ = atan2(surface_normal(1), surface_normal(0));


                    phi_ = latitude * pi / 180;
                    compute_abc(delta_, phi_, slope_, aspect_);
                    costt_ = costt(omega_); // eq.(14)
                    compute_abc(delta_, phi_, 0.0, 0.0);
                    costthor_ = costt(omega_);

                    compute_sun_rise_set(delta_, phi_, slope_, aspect_);
                    if (omega_ > omega1_24_ and omega_ < omega2_24_) {
                        rahor_ = max(0.0, compute_ra(costthor_, doy_)); // eq.(1) with cos(theta)hor
                        //ra_ = min(rahor_,max(0.0,compute_ra(costt_,doy_))); // eq.(1)
                        ra_ = max(0.0, compute_ra(costt_, doy_)); // eq.(1)
                    } else {
                        ra_ = 0.0;
                        rahor_ = 0.0;
                    };

                    double W; //equivalent depth of precipitable water in the atmosphere[mm]
                    eatm_ = atm_pressure(
                            elevation); // [kPa] atmospheric pressure as a function of elevation ///TODO: get elevation from cell.midpoint().z
                    ea_ = actual_vp(temperature, rhumidity); //[kPa] actual vapor pressure
                    W = 0.14 * ea_ * eatm_ + 2.1; // eq.(18)

                    double Kbo;
                    double sin_beta, sin_betahor;
                    sin_betahor = abs(
                            costthor_); // eq.(20) equal to (4), cos_tthor = sin_betahor /// TODO: check if cos_tt = sin_beta is true for inclined surface
                    sin_beta = abs(costt_);
                    // clearness index for direct beam radiation

                    Kbo = min(1.0, max(-0.4, 0.98 * exp(-0.00146 * eatm_ / param.turbidity / sin_beta -
                                                        0.075 * pow((W / sin_beta), 0.4)))); // eq.(17)
                    double Kbohor = min(1.0, max(-0.4, 0.98 * exp(-0.00146 * eatm_ / param.turbidity / sin_betahor -
                                                                  0.075 * pow((W / sin_betahor), 0.4)))); // eq.(17)

                    double Kdo; // transmissivity of diffuse radiation, eq.(19)a,b,c
                    if (Kbo >= 0.15) { Kdo = 0.35 - 0.36 * Kbo; }
                    else if (Kbo < 0.15 and Kbo > 0.065) { Kdo = 0.18 + 0.82 * Kbo; }
                    else { Kdo = 0.10 + 2.08 * Kbo; }

                    double Kdohor;
                    if (Kbohor >= 0.15) { Kdohor = 0.35 - 0.36 * Kbohor; }
                    else if (Kbohor < 0.15 and Kbohor > 0.065) { Kdohor = 0.18 + 0.82 * Kbohor; }
                    else { Kdohor = 0.10 + 2.08 * Kbohor; }

                    double fi_ = fi();//eq.(32)

                    // fb_ = min(5.0,Kbo/Kbohor*ra_/(rahor_>0.0?rahor_:max(0.00001,ra_)));//eq.(34)
                    fb_ = ra_ / (rahor_ > 0.0 ? rahor_ : max(0.00001, ra_));//eq.(34)

                    double fia_ = fia(Kbohor, Kdohor); //eq.(33)

                    if (omega1_24_ > omega2_24_) {
                        omega1_24_ = omega2_24_;
                        ra_ = 0.0;
                    }//slope is always shaded
                    //if ((omega_ > omega2_24b_) and (omega_<omega1_24b_) and (omega1_24_<omega2_24b_) and (omega1_24b_<omega2_24_)){ra_ = 0.0;}
                    // rso_ = max(0.0,Kbo*ra_ + (fia_*Kdo + param.albedo*(1-fi_)*(Kbo+Kdo))*rahor_); // eq.(37)     direct beam + diffuse + reflected, only positive values accepted
                    response.dir_radiation = Kbo * ra_;
                    response.dif_radiation = fia_ * Kdo * rahor_;
                    response.ref_radiation = param.albedo * (1 - fi_) * (Kbo + Kdo) * rahor_;
                    response.psw_radiation = response.dir_radiation + response.dif_radiation +
                                             response.ref_radiation; // clear sky solar radiation for inclined surface [W/m2]
                    response.tsw_radiation = response.psw_radiation; // if no measurements -- than we use predicted value
                    return;
                }

                /**\brief translates measured solar radiation from horizontal surfaces to slopes
                 * ref.: Allen, R. G.; Trezza, R. & Tasumi, M. Analytical integrated functions for daily solar radiation on slopes Agricultural and Forest Meteorology, 2006, 139, 55-73
                 * \param latitude, [deg]
                 * \param utctime,
                 * \param surface_normal
                 * \param temperature, [degC]
                 * \param rhumidity, [percent]
                 * \param elevation, [m]
                 * \param rsm,[W/m^2] -- measured solar radiation
                 * \return */
                template<class V>
                void tsw_radiation(R &response, double latitude, utctime t, const V &surface_normal,
                                   double temperature = 0.0, double rhumidity = 40.0, double elevation = 0.0,
                                   double rsm = 0.0) {
                    // first calculate all predicted values
                    psw_radiation(response, latitude, t, surface_normal, temperature, rhumidity, elevation);
                    double tauswhor = rsm > 0.0 ? rsm / (rahor_ > 0.0 ? rahor_ : rsm)
                                                : 1.0; //? not sure if we use a theoretical rahor here
                    double KBhor;
                    if (tauswhor >= 0.42) { KBhor = 1.56 * tauswhor - 0.55; }
                    else if (tauswhor > 0.175 and tauswhor < 0.42) {
                        KBhor = 0.022 - 0.280 * tauswhor + 0.828 * tauswhor * tauswhor + 0.765 * pow(tauswhor, 3);
                    }
                    else { KBhor = 0.016 * tauswhor; }
                    double KDhor = tauswhor - KBhor;
                    response.tsw_radiation = rsm * (fb_ * KBhor / tauswhor + fia(KBhor, KDhor) * KDhor / tauswhor +
                                                    param.albedo * (1 - fi()));
                    return;
                }

            private:
                double delta_;
                double omega_;
                double phi_;
                double slope_;
                double aspect_;
                double ea_;
                double eatm_;

                double ra_ = 0.0; // extraterrestrial solar radiation for inclined surface[W/m2]
                double rahor_ = 0.0; // extraterrestrial solar radiation for horizontal surfaces

                const double gsc = 1367; // W/m2 -- solar constant
                const double Pa2kPa = 0.001; // Pa to kPa
                const double deg2rad = pi / 180; // degrees to radians
                const double rad2deg = 180 / pi; // rad to deg
                calendar utc;
                double doy_; // day of the yearI
                double lt_; // local time
                double costt_, costthor_;
                double a_, b_, c_;
                double omega1_24_, omega2_24_, omega1_24b_, omega2_24b_; //integration limits, actual periods of sun
                double fb_;

                /** \brief computes necessary geometric parameters
                 * \param omega, [rad] -- hour angle
                 * \return cos(theta) -- theta -- angle of incidence             * */
                double costt(double omega) {
                    return -a_ + b_ * cos(omega) + c_ * sin(omega); // eq.14
                }

                /** \brief computes necessary geometric parameters
                 * \param phi, [rad] -- latitude
                 * \param s, [rad] -- slope angle
                 * \param gamma, [rad] -- aspect angle
                 * \return cos(theta) -- theta -- angle of incidence             * */
                void compute_abc(double delta, double phi, double s = 0.0, double gamma = 0.0) {
                    a_ = sin(delta) * cos(phi) * sin(s) * cos(gamma) - sin(delta) * sin(phi) * cos(s); // eq.11a
                    b_ = cos(delta) * cos(phi) * cos(s) + cos(delta) * sin(phi) * sin(s) * cos(gamma);//eq 11b
                    c_ = cos(delta) * sin(s) * sin(gamma);
                }

                /**\brief compute sun rise and sun set values
                 * \param delta,[rad] - earrh declination angle
                 * \param phi,[rad] -- latitude
                 * \param slope,[rad]
                 * \param aspect,[rad]
                 * calculates  local variables omega1_24_, omega2_24_, omega1_24b_, omega2_24b_*/
                void compute_sun_rise_set(double delta, double phi, double slope, double aspect) {
                    ///TODO get info about hemisphere from data, don't see anything in geopoint, but it is inside netcdf file -- add geopoint_with_crs to interface
                    double omega_s; // omega_s -- time of potential horizontal sunset, -omega_s --time of horizontal sunrize
                    // this solar noon and sunrise taken from ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.575
                    if (abs(phi - delta) >= pi / 2) { omega_s = pi; }
                    if (abs(phi - delta) < pi / 2 and (abs(phi + delta) >= pi / 2)) { omega_s = 0; }
                    if (abs(phi - delta) < pi / 2 and (abs(phi + delta) < pi / 2)) {
                        omega_s = acos(-tan(delta) * tan(phi));
                    } // for horizontal surface

                    compute_abc(delta, phi, slope, aspect);
                    double costt_sunset = costt(omega_s);
                    double costt_sunrise = costt(-omega_s);

                    /// TODO: verify this weird sun rise and sunset procedure

                    double bbcc = b_ * b_ + c_ * c_ > 0.0 ? b_ * b_ + c_ * c_ : b_ * b_ + c_ * c_ + 0.0001;
                    double sqrt_bca = bbcc - a_ * a_ > 0.0 ? bbcc - a_ * a_ : 0.0; // the authors suggest here 0.0001???
                    double sin_omega1 = min(1.0, max(-1.0, (a_ * c_ - b_ * pow(sqrt_bca, 0.5)) / bbcc));//eq.(13a)
                    double omega1 = asin(sin_omega1);
                    omega1_24_ = omega1;
                    double costt_omega1 = costt(omega1);
                    if ((costt_sunrise <= costt_omega1) and (costt_omega1 < 0.001)) { omega1_24_ = omega1; }
                    else {
                        omega1 = -pi - omega1;
                        if (costt(omega1) > 0.001) { omega1_24_ = -omega_s; }
                        else {
                            if (omega1 <= -omega_s) { omega1_24_ = -omega_s; }
                            else { omega1_24_ = -omega1; }
                        }
                    }
                    if (omega1_24_ < -omega_s) { omega1_24_ = -omega_s; }
                    double sin_omega2 = min(1.0, max(-1.0, (a_ * c_ + b_ * pow(sqrt_bca, 0.5)) / bbcc));//eq.(13b)
                    double omega2 = asin(sin_omega2);
                    omega2_24_ = omega2;
                    double costt_omega2 = costt(omega2);
                    if (costt_sunset <= costt_omega2 and costt_omega2 < 0.001) { omega2_24_ = omega2; }
                    else {
                        omega2 = pi - omega2;
                        if (costt(omega2) > 0.001) { omega2_24_ = omega_s; }
                        else {
                            if (omega2 >= omega_s) { omega2_24_ = omega_s; }
                            else { omega2_24_ = omega2; }
                        }
                    }
                    if (omega2_24_ > omega_s) { omega2_24_ = omega_s; }

//
                    // two periods of direct beam radiation (eq.7)
                    if (sin(slope) > sin(phi) * cos(delta) + cos(phi) * sin(delta)) {
                        double sinA = min(1.0, max(-1.0, (a_ * c_ + b_ * pow(sqrt_bca, 0.5)) / bbcc));
                        double A = asin(sinA);
                        double sinB = min(1.0, max(-1.0, (a_ * c_ + b_ * pow(sqrt_bca, 0.5)) / bbcc));
                        double B = asin(sinB);
                        omega2_24b_ = min(A, B);
                        omega1_24b_ = max(A, B);
                        compute_abc(delta_, phi_, slope_, aspect_);
                        double costt_omega2_24b = costt(omega2_24b_);
                        if (costt_omega2_24b < -0.001 or costt_omega2_24b > 0.001) { omega2_24b_ = -pi - omega2_24b_; }
                        double costt_omega1_24b = costt(omega1_24b_);
                        if ((costt_omega1_24b < -0.001 or costt_omega1_24b > 0.001)) { omega1_24b_ = pi - omega1_24b_; }
                        if ((omega2_24b_ > omega1_24_) or (omega1_24b_ < omega2_24_)) {
                            omega2_24b_ = omega1_24_;
                            omega1_24b_ = omega1_24_;
                        } // single period of sun
                    }
                }

//                /** \brief computes standard atmospheric pressure
//                 * \param height, [m] -- elevation of the point
//                 * \return p, [kPa] -- atmospheric pressure */
//                double atm_pressure(double height) { // height < 11000
//                    const double p0 = 101325.0; //[Pa]standard sea level pressure
//                    const double L = 0.0065; //[K/m] temperature lapse rate
//                    const double g = 9.80665; //[m/s2] earth-surface gravitational acceleration
//                    const double R0 = 8.31447;//[J/mol/K] Universal gas constant
//                    const double M = 0.0289644; //[kg/mol] molar mass of dry air
//                    const double T0 = 288.15; // [K] sea level standard temperature
//                    return p0 * pow((1 - L * height / T0), (g * M / R0 / L)) * Pa2kPa;
//                }
//
//                /**\brief computes actual vapor pressure from dewpoint temperature
//                 * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.113
//                 * \param temperature, [degC]
//                 * \param rhumidity, [percent] -- relative humidity
//                 * \return e, [kPa] -- actual vapor pressure*/
//                double actual_vp(double temperature, double rhumidity) {
//                    double es = (temperature >= 0.0) ? (611 * exp(17.27 * temperature / (temperature + 237.3))) : 611 *
//                                                                                                                  exp(21.87 *
//                                                                                                                      temperature /
//                                                                                                                      (temperature +
//                                                                                                                       265.5)); // saturation vapor pressure,[sPa], eq.(3.9)
//                    return rhumidity / 100.0 * es * Pa2kPa;//[kPa], based on eq.(3.12)
//                }

                /**\brief computes solar hour angle from local time
                 * ref.: https://en.wikipedia.org/wiki/Equation_of_time
                 * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, no EOT correction provided
                 * \param lt -- local time (LT) [h]
                 * \param longitute = (0 for UTC time)
                 * \param tz = 0 -- time zone [h], difference of LT from UTC
                 * we use UTC, so longitude = 0.0, tz = 0.0
                 * \return HRA, [rad] -- Earth hour angle*/
                double hour_angle(double lt, double longitude = 0.0, double tz = 0.0) {
                    double LSTM = 15 * tz; // local standard time meridian
//                double B = 360/365*(doy)*(-81);
//                double EoT = 9.87*sin(2*B) - 7.3*cos(B) - 1.5*sin(B);// equation of time
                    double M = 2 * pi / 365.2596358 * doy_ * pi / 180;
                    double EoT = -7.659 * sin(M) +
                                 9.863 * sin(2 * M + 3.5932);//https://en.wikipedia.org/wiki/Equation_of_time*/
                    double TC = 4 * (longitude - LSTM) +
                                EoT; //time correction factor including EoT correction for Earth eccentricity
                    double LST = lt - TC / 60; // local solar time
                    return 15 * (lt - 12) *
                           deg2rad; ///TODO: find right EOT and data for validation, so it should return: 15*(LST-12)
                }

                /**\brief computes earth declination angle
                 * // ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, eq.(D.5)
                 * \param doy -- day of the year
                 * \return declination, [rad] */
                double compute_earth_declination(double doy) {
                    double G = 2 * pi / 365.0 * (doy - 1);
                    return 0.006918 - 0.399912 * cos(G) + 0.070257 * sin(G) - 0.006758 * cos(2 * G) +
                           0.000907 * sin(2 * G) - 0.002697 * cos(3 * G) + 0.00148 * sin(3 * G);
                }

                /**\brief computes extraterrestrial solar radiation
                 * \param cos(theta)
                 * \param doy -- day of the year
                 * return ra, [W/m^2]*/
                double compute_ra(double cos_theta, double doy) {
                    return gsc * cos_theta * (1 + 0.0033 * cos(doy * 2 * pi / 365)); // eq.(1)
                }

                double fi() { return 0.75 + 0.25 * cos(slope_) - 0.5 / pi * slope_;/*eq.(32)*/}

                double fia(double kb, double kd) {
                    return (1 - kb) *
                           (1 + pow(kb / ((kb + kd) > 0.0 ? (kb + kd) : 0.0001), 0.5) * pow(sin(slope_ / 2), 3)) *
                           fi() + fb_ * kb; /*eq.(33)*/
                }

                /**\brief compute 24h parameter*/
                double costt24(double omega1, double omega2) {
                    //return -a_*(omega2-omega1) + b_ * (sin(omega2)-sin(omega1))+c_*(cos(omega2)-cos(omega1));
                    return 2 * (-a_ * (omega2) + b_ * (sin(omega2)) + c_ * (cos(omega2) - 1));
                }
            };


        }

    }
}