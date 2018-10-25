#include "test_pch.h"
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

    namespace radiation {

        struct parameter {
            double albedo = 0.1; // average albedo of the surrounding ground surface:0.15-0.25 -- grass, 0.10-0.15 -- coniferous forest, 0.15 - 0.25 -- deciduous forest, 0.04-0.08 -- open water, 0.15-0.35 -- bare soil
            parameter(double albedo = 0.1):albedo(albedo){}
        };

        //struct state {}; // No state variables for this method

        struct response {
            double cs_radiation = 0.0;
        };

        template <class P>
        struct calculator {
            P param;
            explicit calculator(const P& p):param(p) {}

            double declination() const {return delta_*180/pi;}  // declination of the earth, [deg], positive for northern hemisphere summer
            double hra() const {return omega_;} // earth hour angle, [deg]

            double slope() const {return slope_;} // horizontal slope, should be calculated from normal vector
            double aspect() const {return aspect_;} // surface aspect angle, 0 -- due south, -pi/2 -- due east, pi/2 -- due west, +-pi -- due north

            double latitude() const {return phi_;}// latitude, should be available from cell?/// TODO: add PROJ4 for conversion from cartesian to wgs84

            double ra = 0.0; // extraterrestrial solar radiation for inclined surface[W/m2]
            double rahor = 0.0; // extraterrestrial solar radiation for horizontal surfaces
            double rso = 0.0; // clear sky solar radiation for inclined surface [W/m2]

            double patm() const {return eatm_;}
            double ea() const {return ea_;}



            /** \brief computes instantaneous extraterrestrial solar radiation
             * \param latitude (phi)
             * \param time
             * \param surface normal
             * \return ra [W/m2]*/

            template <class V>
            double ra_radiation(double latitude, utctime t, const V & surface_normal, double temperature=0.0, double rhumidity=40.0, double elevation = 0.0 ){
                compute_day_and_time(t);
                // ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, eq.(D.5)
                double G = 2*pi/365.0*(doy_-1);
                delta_ = 0.006918 - 0.399912*cos(G)+0.070257*sin(G)- 0.006758*cos(2*G)+0.000907*sin(2*G) - 0.002697*cos(3*G)+0.00148*sin(3*G);
                omega_ = hour_angle(lt_); // earth hour angle
                costt_ = costt(omega_, delta_, latitude, slope_, aspect_);
                costthor_ = costt(omega_, delta_, latitude);
                ra = gsc*costt_*(1+0.0033*cos(doy_*2*pi/365)); // eq.(1)
                rahor = gsc*costthor_*(1+0.0033*cos(doy_*2*pi/365)); // eq.(1) with cos(theta)hor

                ///TODO: clean up
                eatm_ = atm_pressure(elevation); // [kPa] atmospheric pressure as a function of elevation
                ea_ = actual_vap_pressure(temperature, rhumidity); //[kPa] actual vapor pressure

                return ra;
            }
            /** \brief computes instantaneous clear-sky radiation  for inclined surfaces
             * \param temperature, [degC]
             * \param rhumidity, [percent]
             * \param time -- utctime
             * \param cell
             * \return rso, [W/m^2] -- clear sky radiation*/
            template <class C>
            double rso_cs_radiation(double temperature, double rhumidity, utctime t, const C& cell){
                omega_ = hour_angle(utc.calendar_units(t).hour); // hour angle
                double W; //equivalent depth of precipitable water in the atmosphere[mm]
                eatm_ = atm_pressure(cell.midpoint().z); // [kPa] atmospheric pressure as a function of elevation
                ea_ = actual_vap_pressure(temperature, rhumidity); //[kPa] actual vapor pressure
                W = 0.14*ea_*eatm_ + 2.1; // eq.(18)
                double Kt = 1.0; // turbidity, 1.0 - clean, 0.5 -- dusty
                double Kbo;
                double sin_beta, sin_betahor;
                sin_betahor = costthor_; // eq.(20) equal to (4), cos_tthor = sin_betahor /// TODO: check if cos_tt = sin_beta is true for inclined surface
                sin_beta = costt_;
                // clearness index for direct beam radiation
                Kbo = 0.98 * exp(-0.00146*eatm_/Kt/sin_beta - 0.075*pow(W/sin_beta,0.4)); // eq.(17)
                double Kbohor = 0.98 * exp(-0.00146*eatm_/Kt/sin_betahor - 0.075*pow(W/sin_betahor,0.4)); // eq.(17)
                double Kdo; // transmissivity of diffuse radiation, eq.(19)a,b,c
                if (Kbo >= 0.15){Kdo = 0.35 - 0.36*Kbo;}
                else if (Kbo < 0.15 and Kbo > 0.065){Kdo = 0.18 + 0.82*Kbo;}
                else {Kdo = 0.10 + 2.08 * Kbo;}
                double Kdohor;
                if (Kbohor >= 0.15){Kdohor = 0.35 - 0.36*Kbohor;}
                else if (Kbohor < 0.15 and Kbohor > 0.065){Kdohor = 0.18 + 0.82*Kbohor;}
                else {Kdohor = 0.10 + 2.08 * Kbohor;}

                double fi = 0.75 + 0.25*cos(slope_) - 0.5/pi*slope_;//eq.(32)

                double fb = Kbo/Kbohor*ra/rahor;//eq.(34)
                double fia = (1 - Kbohor)*(1+pow(Kbohor/(Kbohor+Kdohor),0.5)*pow(sin(slope_/2),3))*fi + fb*Kbohor; //eq.(33)
                rso = Kbo*ra + (fia*Kdo + param.albedo*(1-fi)*(Kbo+Kdo))*rahor; // eq.(37)     direct beam + diffuse + reflected
                return rso; // eq.(37)
            }
        private:
            double delta_;
            double omega_;
            double phi_;
            double slope_;
            double aspect_;
            double ea_;
            double eatm_;

            const double gsc = 1367; // W/m2 -- solar constant
            const double PatokPa = 0.001; // Pa to kPa
            calendar utc;
            double doy_; // day of the year
            double lt_; // local time
            double costt_, costthor_;
            double a,b,c;
            /** \brief computes necessary geometric parameters
             * \param omega, [deg] -- hour angle
             * \param phi, [deg] -- latitude
             * \param s, [deg] -- slope angle
             * \param gamma, [deg] -- aspect angle
             * \return cos(theta) -- theta -- angle of incidence             * */
            double costt(double omega,  double delta, double phi, double s = 0.0, double gamma = 0.0){
                a = sin(delta)*cos(phi)*sin(s)*cos(gamma) - sin(delta)*sin(phi)*cos(s); // eq.11a
                b = cos(delta)*cos(phi)*cos(s) + cos(delta)*sin(phi)*sin(s)*cos(gamma);//eq 11b
                c = cos(delta)*sin(s)*sin(gamma);
                return  -a + b*cos(omega) + c*sin(omega); // eq.14
            }
            double sinomega(){
                return 0.0;
            }
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
                return p0*pow((1 - L*height/T0),(g*M/R0/L))*PatokPa;
            }
            /**\brief computes actual vapor pressure from dewpoint temperature
             * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.113
             * \param temperature, [degC]
             * \param rhumidity, [percent] -- relative humidity
             * \return e, [kPa] -- actual vapor pressure*/
            double actual_vap_pressure(double temperature, double rhumidity){
                double es = (temperature>=0.0) ? (611 * exp(17.27*temperature/(temperature + 237.3))) : 611*exp(21.87*temperature/(temperature+265.5)) ; // saturation vapor pressure,[sPa], eq.(3.9)
                return rhumidity/100.0*es*PatokPa;//[kPa], based on eq.(3.12)
            }
            /**\brief computes solar hour angle from local time
             * ref.: https://en.wikipedia.org/wiki/Equation_of_time
             * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, no EOT correction provided
             * \param lt -- local time (LT) [h]
             * \param longitute = (0 for UTC time)
             * \param tz = 0 -- time zone [h], difference of LT from UTC
             * we use UTC, so longitude = 0.0, tz = 0.0
             * \return HRA, [deg] -- Earth hour angle*/
            double hour_angle(double lt, double longitude=0.0, double tz=0.0){
                double LSTM = 15*tz; // local standard time meridian
//                double B = 360/365*(doy)*(-81);
//                double EoT = 9.87*sin(2*B) - 7.3*cos(B) - 1.5*sin(B);// equation of time
                double M = 2*pi/365.2596358*doy_*pi/180;
                double EoT = -7.659*sin(M)+9.863*sin(2*M+3.5932);//https://en.wikipedia.org/wiki/Equation_of_time*/
                double TC = 4*(longitude - LSTM) + EoT; //time correction factor including EoT correction for Earth eccentricity
                double LST = lt - TC/60; // local solar time
                return 15*(lt-12); ///TODO: find right EOT and data for validation, so it should return: 15*(LST-12)
            }
            void compute_day_and_time(utctime t){
                doy_ = utc.day_of_year(t);
                lt_ = utc.calendar_units(t).hour + utc.calendar_units(t).minute/60.0;
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
    using shyft::core::radiation::parameter;
    using shyft::core::radiation::calculator;
    using shyft::core::radiation::surface_normal;
    using shyft::core::calendar;
    using shyft::test::cell;
    using std::vector;
    using shyft::core::utctime;
    // test basics: creation, etc
    // /TODO: prepare set of tests: 1. check vapor pressure result, 2. check atm_pressure vs elevation

   /* TEST_CASE("check_doy"){
        parameter p;
        calculator r(p);
        calendar utc_cal;
        double lat = 56.0;
        utctime t1, t2, t3;
        t1 = utc_cal.time(1970, 1, 1, 10, 30, 0, 0);
        t2 = utc_cal.time(1970, 6, 1, 12, 0, 0, 0);
        t3 = utc_cal.time(1970, 9, 1, 23, 30, 0, 0);
        r.ra_radiation(lat, t1, 1);
        FAST_CHECK_EQ(r.doy_, doctest::Approx(1.0));// doy 1 January
        FAST_CHECK_EQ(r.lt_, doctest::Approx(10.50));// hour 1 January
        r.ra_radiation(lat, t2, 1);
        FAST_CHECK_EQ(r.doy_, doctest::Approx(152.0));// doy 1 June
        r.ra_radiation(lat, t3, 1);
        FAST_CHECK_EQ(r.doy_, doctest::Approx(244.0));// doy 1 September

    }*/ // no need for checking doy anymore
    TEST_CASE("check_hour_angle"){
        /** \brief Check Earth Hour Angle
         * based on ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, no EOT correction provided*/
        parameter p;
        calculator r(p);
        calendar utc_cal;
        double lat = 56.0;
        utctime t1, t2, t3;
        t1 = utc_cal.time(1970, 1, 1, 10, 30, 0, 0);
        t2 = utc_cal.time(1970, 6, 1, 12, 0, 0, 0);
        t3 = utc_cal.time(1970, 9, 1, 23, 30, 0, 0);
        // using simple formulae for checking, need to provide test for EOT
        r.ra_radiation(lat, t1, 1);
                FAST_CHECK_EQ(r.hra(), doctest::Approx(-22.5));// earth hour angle
        r.ra_radiation(lat, t2, 1);
                FAST_CHECK_EQ(r.hra(), doctest::Approx(0.0));// earth hour angle
        r.ra_radiation(lat, t3, 1);
                FAST_CHECK_EQ(r.hra(), doctest::Approx(172.5));// earth hour angle

    }
    TEST_CASE("check_declination_angle"){
        /** \brief Check Earth declination angle
         * based on ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, eq.(D.5) */
        parameter p;
        calculator r(p);
        calendar utc_cal;
        double lat = 56.0;
        utctime t1, t2, t3, t4;
        t1 = utc_cal.time(1970, 12, 21, 12, 30, 0, 0);
        t2 = utc_cal.time(1970, 6, 21, 12, 30, 0, 0);
        t3 = utc_cal.time(1970, 3, 21, 12, 30, 0, 0);
        t4 = utc_cal.time(1970, 9, 24, 12, 30, 0, 0);
        r.ra_radiation(lat, t1, 1);
                FAST_CHECK_EQ(r.declination(), doctest::Approx(-23.5).epsilon(0.1));// earth declination angle 21/22 December (min)
        r.ra_radiation(lat, t2, 1);
                FAST_CHECK_EQ(r.declination(), doctest::Approx(23.5).epsilon(0.1));// earth declination angle 21/22 June (max)
        r.ra_radiation(lat, t3, 1);
                FAST_CHECK_EQ(r.declination(), doctest::Approx(0.0).epsilon(0.1));// earth declination angle 20/21 March 0
        r.ra_radiation(lat, t4, 1);
                FAST_CHECK_EQ(r.declination(), doctest::Approx(0.0).epsilon(0.5));// earth declination angle 22/23 September 0 this gives a bit higher error
    }
    TEST_CASE("check_pressures"){
        parameter p;
        calculator r(p);
        calendar utc_cal;
        double lat = 56.0;
        utctime t;
        t = utc_cal.time(1970, 12, 21, 12, 30, 0, 0);
        r.ra_radiation(lat, t, 1, 0.0, 0.0, 0.0);
                FAST_CHECK_EQ(r.patm(), doctest::Approx(101.325).epsilon(0.1));
        r.ra_radiation(lat, t, 1, 0.0, 0.0, 1000.0);
                FAST_CHECK_EQ(r.patm(), doctest::Approx(89.88).epsilon(0.1));
        r.ra_radiation(lat, t, 1, 0.0, 0.0,  -1000.0);
                FAST_CHECK_EQ(r.patm(), doctest::Approx(113.93).epsilon(0.1));

    }

    TEST_CASE("surface_normal_from_cells") {
        vector<cell> cells;
        auto r= surface_normal(cells);
    }
}
