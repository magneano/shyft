#include "test_pch.h"
#include "core/utctime_utilities.h"
#include "core/geo_cell_data.h"
#include <armadillo>
#include <vector>
#include <chrono>
#include <boost/math/constants/constants.hpp>

namespace shyft::core {
    using std::vector;
    using std::cos;
    using std::sin;
    using std::pow;
    using std::exp;
    using namespace std;
    const double pi = boost::math::constants::pi<double>();

    namespace radiation {

        struct calculator {

            double doy; // day of the year
            double delta = 0.1; // declination of the earth, positive for northern hemisphere summer
            double omega; // earth hour angle

            double s = 0.0; // horizontal slope, should be calculated from normal vector
            double gamma = 0.0; // surface aspect angle, 0 -- due south, -pi/2 -- due east, pi/2 -- due west, +-pi -- due north

            double phi = 0.0 ;// latitude, should be available from cell?

            double ra = 0.0; // extraterrestrial solar radiation for inclined surface[W/m2]
            double rahor = 0.0; // extraterrestrial solar radiation for horizontal surfaces
            double rso = 0.0; // clear sky solar radiation for inclined surface [W/m2]

            double alpha = 0.1; // albedo
            string local_time_str;
            double local_time;

            /** \brief computes instantaneous extraterrestrial solar radiation
             * \param latitude (phi)
             * \param time
             * \param surface normal */

            template <class V>
            double ra_radiation(double latitude, double longitude,  utctime t, const V & surface_normal){
                compute_doy(t);
                auto c = utc.calendar_units(t);
                auto us = c.micro_second;
                local_time_str = utc.to_string(t);
                local_time = utc.calendar_units(t).hour + utc.calendar_units(t).minute/60.0;
                // https://en.wikipedia.org/wiki/Position_of_the_Sun
//                delta = -0.39779*cos(((0.98565*(doy + 10))+1.914*sin(0.98565*(doy-2)*pi/180.0))*pi/180.0);
                // ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574
                double G = 2*pi/365.0*(doy-1);
                delta = 0.006918 - 0.399912*cos(G)+0.070257*sin(G)- 0.006758*cos(2*G)+0.000907*sin(2*G) - 0.002697*cos(3*G)+0.00148*sin(3*G);
                // delta = 0.39795*cos(0.98563*(doy-173)*pi/180.0); //http://mypages.iit.edu/~maslanka/SolarGeo.pdf -> this one gives -22.8 <= delta < 22.8
                omega = hour_angle(utc.calendar_units(t).hour+utc.calendar_units(t).minute/60.0); // earth hour angle
                compute_costt(omega, latitude);
                ra = gsc*cos_tt*(1+0.0033*cos(doy*2*pi/365)); // eq.(1)
                rahor = gsc*cos_tthor*(1+0.0033*cos(doy*2*pi/365)); // eq.(1) with cos(theta)hor
                // std::cout<<ra<<std::endl;
                return ra;
            }
            /** \brief computes instantaneous clear-sky radiation  for inclined surfaces*/
            template <class C>
            double rso_cs_radiation(utctime t, const C& cell){
                omega = hour_angle(utc.calendar_units(t).hour); // hour angle
                double W; //equivalent depth of precipitable water in the atmosphere[mm]
                double P = 101.325;// [kPa] atmospheric pressure /// TODO calculate as a function of elevation (should be available for cell) cell.midpoint().z?
                double ea = 101.325; //[kPa] actual vapor pressure ///TODO: calculate mean saturation pressure from temperature, next actual vapor pressure from relative humidity data
                W = 0.14*ea*P + 2.1; // eq.(18)
                double Kt = 1.0; // turbidity, 1.0 - clean, 0.5 -- dusty
                double Kbo;
                double sin_beta, sin_betahor;
                sin_betahor = cos_tthor; // eq.(20) equal to (4), cos_tthor = sin_betahor /// TODO: check if cos_tt = sin_beta is true for inclined surface
                sin_beta = cos_tt;
                // clearness index for direct beam radiation
                Kbo = 0.98 * exp(-0.00146*P/Kt/sin_beta - 0.075*pow(W/sin_beta,0.4)); // eq.(17)
                double Kbohor = 0.98 * exp(-0.00146*P/Kt/sin_betahor - 0.075*pow(W/sin_betahor,0.4)); // eq.(17)
                double Kdo; // transmissivity of diffuse radiation, eq.(19)a,b,c
                if (Kbo >= 0.15){Kdo = 0.35 - 0.36*Kbo;}
                else if (Kbo < 0.15 and Kbo > 0.065){Kdo = 0.18 + 0.82*Kbo;}
                else {Kdo = 0.10 + 2.08 * Kbo;}
                double Kdohor;
                if (Kbohor >= 0.15){Kdohor = 0.35 - 0.36*Kbohor;}
                else if (Kbohor < 0.15 and Kbohor > 0.065){Kdohor = 0.18 + 0.82*Kbohor;}
                else {Kdohor = 0.10 + 2.08 * Kbohor;}

                double fi = 0.75 + 0.25*cos(s) - 0.5/pi*s;//eq.(32)

                double fb = Kbo/Kbohor*ra/rahor;//eq.(34)
                double fia = (1 - Kbohor)*(1+pow(Kbohor/(Kbohor+Kdohor),0.5)*pow(sin(s/2),3))*fi + fb*Kbohor; //eq.(33)
                rso = Kbo*ra + (fia*Kdo + alpha*(1-fi)*(Kbo+Kdo))*rahor; // eq.(37)     direct beam + diffuse + reflected
                return rso; // eq.(37)
            }
        private:
            const double gsc = 1367; // W/m2 -- solar constant    
            calendar utc;
            double a, b, c, cos_tt, cos_tthor;
            /** \brief computes necessary geometric parameters
             * \param omega -- hour angle
             * */
            void compute_costt(double omega, double phi){
                a = sin(delta)*cos(phi)*sin(s)*cos(gamma) - sin(delta)*sin(phi)*cos(s); // eq.11a
                b = cos(delta)*cos(phi)*cos(s) + cos(delta)*sin(phi)*sin(s)*cos(gamma);//eq 11b
                c = cos(delta)*sin(s)*sin(gamma);
                cos_tt = -a + b*cos(omega) + c*sin(omega); // eq.14
                cos_tthor = sin(delta)*sin(phi) + cos(delta)*cos(phi)*cos(omega);
            }
            /** \brief computes standard atmospheric pressure
             * \param height -- [m] elevation of the point */
            double atm_pressure(double height){ // heoght < 11000
                const double p0 = 101325.0; //standard sea level pressure
                const double L = 0.0065; //[K/m] temperature lapse rate
                const double g = 9.80665; //[m/s2] earth-surface gravitational acceleration
                const double R0 = 8.31447;//[J/mol/K] Universal gas constant
                const double M = 0.0289644; //[kg/mol] molar mass of dry air
                const double T0 = 288.15; // sea level standard temperature
                return p0*pow((1 - L*height/T0),(g*M/R0/L));
            }
            /**\brief computes actual vapor pressure from dewpoint temperature
             * \param temperature
             * \param humidity
             * ///TODO: compute using more accurate formula*/
            double actual_vap_pressure(double temperature, double humidity){
                double es = 0.6108 * exp(17.27*temperature/(temperature + 237.3));
                return humidity/100*es;
            }
            /**\brief computes solar hour angle from local time
             * ref.: https://en.wikipedia.org/wiki/Equation_of_time
             * \param lt -- local time (LT) [h]
             * \param longitute = (0 for UTC time)
             * \param tz = 0 -- time zone [h], difference of LT from UTC
             * we use UTC, so longitude = 0.0, tz = 0.0 */
            double hour_angle(double lt, double longitude=0.0, double tz=0.0){
                double LSTM = 15*tz; // local standard time meridian
//                double B = 360/365*(doy)*(-81);
//                double EoT = 9.87*sin(2*B) - 7.3*cos(B) - 1.5*sin(B);// equation of time
                double M = 2*pi/365.2596358*doy*pi/180;
                double EoT = -7.659*sin(M)+9.863*sin(2*M+3.5932);//https://en.wikipedia.org/wiki/Equation_of_time*/

                double TC = 4*(longitude - LSTM) + EoT; //time correction factor including EoT correction for Earth eccentricity
                double LST = lt - TC/60; // local solar time

                return 15*(lt-12); ///TODO: find right EOT and data for validation
            }
            void compute_doy(utctime t){
                doy = utc.day_of_year(t);
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

namespace shyft::test {
    using shyft::core::geo_cell_data;
    struct cell {
        geo_cell_data geo;
    };
}

TEST_SUITE("radiation") {
    using shyft::core::radiation::calculator;
    using shyft::core::radiation::surface_normal;
    using shyft::core::calendar;
    using shyft::test::cell;
    using std::vector;
    using shyft::core::utctime;
    // test basics: creation, etc
    // /TODO: prepare set of tests: 1. check vapor pressure result, 2. check atm_pressure vs elevation

    TEST_CASE("check_doy"){
        calculator r;
        calendar utc_cal;
        double lat = 56.0;
        double lon = 66.31;
        utctime t1, t2, t3;
        t1 = utc_cal.time(1970, 1, 1, 10, 30, 0, 0);
        t2 = utc_cal.time(1970, 6, 1, 12, 0, 0, 0);
        t3 = utc_cal.time(1970, 9, 1, 23, 30, 0, 0);
        r.ra_radiation(lat, lon, t1, 1);
        FAST_CHECK_EQ(r.doy, doctest::Approx(1.0));// doy 1 January
        FAST_CHECK_EQ(r.local_time, doctest::Approx(10.50));// hour 1 January
        r.ra_radiation(lat, lon, t2, 1);
        FAST_CHECK_EQ(r.doy, doctest::Approx(152.0));// doy 1 June
        r.ra_radiation(lat, lon, t3, 1);
        FAST_CHECK_EQ(r.doy, doctest::Approx(244.0));// doy 1 September

    }
    TEST_CASE("check_hour_angle"){
        /** compare to http://www.jgiesen.de/astro/suncalc/index.htm,  if withput EOT the hour angle is calculated as 15*(LT-12),
         * so for 10.30 it is -22.5, see:https://en.wikipedia.org/wiki/Hour_angle*/
        calculator r;
        calendar utc_cal;
        double lat = 56.0;
        double lon = 66.31;
        utctime t1, t2, t3;
        t1 = utc_cal.time(1970, 1, 1, 10, 30, 0, 0);
        t2 = utc_cal.time(1970, 6, 1, 12, 0, 0, 0);
        t3 = utc_cal.time(1970, 9, 1, 23, 30, 0, 0);
        // using simple formulae for checking, needto provide test for EOT
        r.ra_radiation(lat, lon, t1, 1);
                FAST_CHECK_EQ(r.omega, doctest::Approx(-22.5));// earth hour angle
        r.ra_radiation(lat, lon, t2, 1);
                FAST_CHECK_EQ(r.omega, doctest::Approx(0.0));// earth hour angle
        r.ra_radiation(lat, lon, t3, 1);
                FAST_CHECK_EQ(r.omega, doctest::Approx(172.5));// earth hour angle

    }
    TEST_CASE("check_declination_angle"){
        /** The equation from referenced book gives the values closest to 23.5 */
        calculator r;
        calendar utc_cal;
        double lat = 56.0;
        double lon = 66.31;
        utctime t1, t2, t3, t4;
        t1 = utc_cal.time(1970, 12, 21, 12, 30, 0, 0);
        t2 = utc_cal.time(1970, 6, 21, 12, 30, 0, 0);
        t3 = utc_cal.time(1970, 3, 21, 12, 30, 0, 0);
        t4 = utc_cal.time(1970, 9, 24, 12, 30, 0, 0);
        r.ra_radiation(lat, lon, t1, 1);
                FAST_CHECK_EQ(r.delta*180/shyft::core::pi, doctest::Approx(-23.5).epsilon(0.1));// earth declination angle 21/22 December (min)
        r.ra_radiation(lat, lon, t2, 1);
                FAST_CHECK_EQ(r.delta*180/shyft::core::pi, doctest::Approx(23.5).epsilon(0.1));// earth declination angle 21/22 June (max)
        r.ra_radiation(lat, lon, t3, 1);
                FAST_CHECK_EQ(r.delta*180/shyft::core::pi, doctest::Approx(0.0).epsilon(0.1));// earth declination angle 20/21 March 0
        r.ra_radiation(lat, lon, t4, 1);
                //FAST_CHECK_EQ(r.delta*180/shyft::core::pi, doctest::Approx(0.0).epsilon(0.1));// earth declination angle 22/23 September 0 this gives a bit higher error
    }

    TEST_CASE("surface_normal_from_cells") {
        vector<cell> cells;
        auto r= surface_normal(cells);
    }
}
