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

            /// TODO: remove accessing functions after all tests
            double declination() const {return delta_*rad2deg;}  // declination of the earth, [deg], positive for northern hemisphere summer
            double hra() const {return omega_*rad2deg;} // earth hour angle, [deg]

            double slope() const {return slope_*rad2deg;} // horizontal slope, [deg], should be calculated from normal vector
            double aspect() const {return aspect_*rad2deg;} // surface aspect angle, [deg], 0 -- due south, -pi/2 -- due east, pi/2 -- due west, +-pi -- due north

            double latitude() const {return phi_*rad2deg;}// latitude, [deg] should be available from cell?/// TODO: add PROJ4 for conversion from cartesian to wgs84

            double extrater_radiation() const {return ra_ ;} // extraterrestrial solar radiation for inclined surface[W/m2]
            double extrater_radiation_hor() {return rahor_;} // extraterrestrial solar radiation for horizontal surfaces
            double cs_radiation() const {return rso_ ;} // clear sky solar radiation for inclined surface [W/m2]

            double patm() const {return eatm_;}
            double ea() const {return ea_;}

            double sun_rise() const {return omega1_24_*rad2deg;}
            double sun_set() const {return omega2_24_*rad2deg;}


            /** \brief computes instantaneous clear-sky radiation  for inclined surfaces
             * \param temperature, [degC]
             * \param rhumidity, [percent]
             * \param time -- utctime
             * \param cell
             * \return rso, [W/m^2] -- clear sky radiation*/
            template <class V>
            double rso_cs_radiation(double latitude, utctime t, const V & surface_normal, double temperature=0.0, double rhumidity=40.0, double elevation = 0.0 ){
                doy_ = utc.day_of_year(t);
                lt_ = utc.calendar_units(t).hour + utc.calendar_units(t).minute/60.0;
                delta_ = compute_earth_declination(doy_);
                omega_ = hour_angle(lt_); // earth hour angle

                arma::vec normal0{0.0,0.0,1.0}; // normal to horizontal plane
                slope_ = acos(arma::norm_dot(surface_normal, normal0));
                //slope_ = atan2(pow(pow(surface_normal(0),2)+pow(surface_normal(1),2),0.5),surface_normal(2));
                aspect_ = pi - atan2(surface_normal(1),surface_normal(0));

                compute_abc(delta_,latitude,slope_,aspect_);
                costt_ = costt(omega_); // eq.(14)
                compute_abc(delta_,latitude);
                costthor_ = costt(omega_);

                ra_ = compute_ra(costt_,doy_); // eq.(1)
                rahor_ = compute_ra(costthor_,doy_); // eq.(1) with cos(theta)hor

                double W; //equivalent depth of precipitable water in the atmosphere[mm]
                eatm_ = atm_pressure(elevation); // [kPa] atmospheric pressure as a function of elevation ///TODO: get elevation from cell.midpoint().z
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

                double fb = Kbo/Kbohor*ra_/rahor_;//eq.(34)
                double fia = (1 - Kbohor)*(1+pow(Kbohor/(Kbohor+Kdohor),0.5)*pow(sin(slope_/2),3))*fi + fb*Kbohor; //eq.(33)
                rso_ = Kbo*ra_ + (fia*Kdo + param.albedo*(1-fi)*(Kbo+Kdo))*rahor_; // eq.(37)     direct beam + diffuse + reflected
                compute_sun_rise_set(delta_,latitude,slope_,aspect_);
                if (omega_ < omega1_24_ or omega_>omega2_24_){rso_=0.0;}
                return rso_; // eq.(37)
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
            double rso_ = 0.0; // clear sky solar radiation for inclined surface [W/m2]

            const double gsc = 1367; // W/m2 -- solar constant
            const double Pa2kPa = 0.001; // Pa to kPa
            const double deg2rad=pi/180; // degrees to radians
            const double rad2deg=180/pi; // rad to deg
            calendar utc;
            double doy_; // day of the yearI
            double lt_; // local time
            double costt_, costthor_;
            double a_,b_,c_;
            double omega1_24_,omega2_24_,omega1_24b_,omega2_24b_; //integration limits, actual periods of sun

            /** \brief computes necessary geometric parameters
             * \param omega, [rad] -- hour angle
             * \return cos(theta) -- theta -- angle of incidence             * */
            double costt(double omega){
                return  -a_ + b_*cos(omega) + c_*sin(omega); // eq.14
            }
            /** \brief computes necessary geometric parameters
             * \param phi, [rad] -- latitude
             * \param s, [rad] -- slope angle
             * \param gamma, [rad] -- aspect angle
             * \return cos(theta) -- theta -- angle of incidence             * */
            void compute_abc(double delta, double phi, double s = 0.0, double gamma = 0.0){
                a_ = sin(delta)*cos(phi)*sin(s)*cos(gamma) - sin(delta)*sin(phi)*cos(s); // eq.11a
                b_ = cos(delta)*cos(phi)*cos(s) + cos(delta)*sin(phi)*sin(s)*cos(gamma);//eq 11b
                c_ = cos(delta)*sin(s)*sin(gamma);
            }
            /**\brief compute sun rise and sun set values
             * \param delta,[rad] - earrh declination angle
             * \param phi,[rad] -- latitude
             * \param slope,[rad]
             * \param aspect,[rad]
             * \param hemisphere_str -- "N" or "S"
             * calculates  local variables omega1_24_, omega2_24_, omega1_24b_, omega2_24b_*/
            void compute_sun_rise_set(double delta, double phi, double slope, double aspect,std::string hemisphere="N"){
                ///TODO get info about hemisphere from data, don't see anything in geopoint, but it is inside netcdf file
                double omega_s; // omega_s -- time of potential horizontal sunset, -omega_s --time of horizontal sunrize
                if (hemisphere.compare("N")) {
                    if ((delta + phi) > pi / 2) { omega_s = pi; } // for northern hemisphere
                    if ((delta - phi) > pi / 2) { omega_s = 0; }
                }
                else {
                    if ((delta + phi) < -pi / 2) { omega_s = pi; } // for southern hemisphere
                    if ((delta - phi) < -pi / 2) { omega_s = 0; }
                }
                compute_abc(delta,phi,slope,aspect);
                double costt_sunset = costt(omega_s);
                double costt_sunrise = costt(-omega_s);

                /// TODO: verify this weird sun rise and sunset procedure
                // commented stuff is supposed to be calculated, but it gives zeroes

                double bbcc = b_*b_ + c_*c_ > 0.0 ? b_*b_+c_*c_ : b_*b_+c_*c_+0.0001;
                double sqrt_bca = bbcc-a_*a_>0.0 ? bbcc-a_*a_ : 0.0; // the authors suggest here 0.0001???
                double sin_omega1 = min(1.0, max(-1.0,(a_*c_ - b_*pow(sqrt_bca,0.5))/bbcc));//eq.(13a)
                double omega1 = asin(sin_omega1);
                omega1_24_=omega1;
//                double costt_omega1 = costt(omega1);
//                if ((costt_sunrise <= costt_omega1) and (costt_omega1<0.001)){omega1_24_=omega1;}
//                else {
//                    omega1 = -pi-omega1;
//                    if (costt(omega1) > 0.001){omega1_24_ = -omega_s;}
//                    else {
//                        if (omega1<=-omega_s){omega1_24_ = -omega_s;}
//                        else {omega1_24_ = -omega1;}
//                    }
//                }
//                if (omega1_24_<-omega_s){omega1_24_ = -omega_s;}
                double sin_omega2 = min(1.0,max(-1.0,(a_*c_ + b_*pow(sqrt_bca,0.5))/bbcc));//eq.(13b)
                double omega2 = asin(sin_omega2);
                omega2_24_=omega2;
//                double costt_omega2 = costt(omega2);
//                if (costt_sunset<=costt_omega2 and costt_omega2<0.001){omega2_24_ = omega2;}
//                else {
//                    omega2 = pi-omega2;
//                    if (costt(omega2)>0.001){omega2_24_ = omega_s;}
//                    else{
//                        if (omega2>=omega_s){omega2_24_ = omega_s;}
//                        else{omega2_24_ = omega2;}
//                    }
//                }
//                if (omega2_24_ > omega_s){omega2_24_ = omega_s;}
//                if (omega1_24_ < omega2_24_) {omega1_24_ = omega2_24_;}//slope is always shaded
//
//                // two periods of direct beam radiation (eq.7)
//                if (sin(slope)>sin(phi)*cos(delta)+cos(phi)*sin(delta)){
//                    double sinA = min(1.0,max(-1.0,(a_*c_ + b_*pow(sqrt_bca,0.5))/bbcc));
//                    double A = asin(sinA);
//                    double sinB = min(1.0,max(-1.0,(a_*c_ + b_*pow(sqrt_bca,0.5))/bbcc));
//                    double B = asin(sinB);
//                    omega2_24b_ = min(A,B);
//                    omega1_24b_ = max(A,B);
//                    double costt_omega2_24b = costt(omega2_24b_);
//                    if (costt_omega2_24b<-0.001 or costt_omega2_24b>0.001){omega2_24b_ = -pi-omega2_24b_;}
//                    double costt_omega1_24b=costt(omega1_24b_);
//                    if ((costt_omega1_24b<-0.001 or costt_omega1_24b>0.001)){omega1_24b_ = pi-omega1_24b_;}
//                    if ((omega2_24b_<omega1_24_) or (omega1_24b_<omega2_24_)){
//                        omega2_24b_ = omega1_24_;
//                        omega1_24b_ = omega1_24_;} // single period of sun
//                }
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
                return p0*pow((1 - L*height/T0),(g*M/R0/L))*Pa2kPa;
            }
            /**\brief computes actual vapor pressure from dewpoint temperature
             * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.113
             * \param temperature, [degC]
             * \param rhumidity, [percent] -- relative humidity
             * \return e, [kPa] -- actual vapor pressure*/
            double actual_vap_pressure(double temperature, double rhumidity){
                double es = (temperature>=0.0) ? (611 * exp(17.27*temperature/(temperature + 237.3))) : 611*exp(21.87*temperature/(temperature+265.5)) ; // saturation vapor pressure,[sPa], eq.(3.9)
                return rhumidity/100.0*es*Pa2kPa;//[kPa], based on eq.(3.12)
            }
            /**\brief computes solar hour angle from local time
             * ref.: https://en.wikipedia.org/wiki/Equation_of_time
             * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, no EOT correction provided
             * \param lt -- local time (LT) [h]
             * \param longitute = (0 for UTC time)
             * \param tz = 0 -- time zone [h], difference of LT from UTC
             * we use UTC, so longitude = 0.0, tz = 0.0
             * \return HRA, [rad] -- Earth hour angle*/
            double hour_angle(double lt, double longitude=0.0, double tz=0.0){
                double LSTM = 15*tz; // local standard time meridian
//                double B = 360/365*(doy)*(-81);
//                double EoT = 9.87*sin(2*B) - 7.3*cos(B) - 1.5*sin(B);// equation of time
                double M = 2*pi/365.2596358*doy_*pi/180;
                double EoT = -7.659*sin(M)+9.863*sin(2*M+3.5932);//https://en.wikipedia.org/wiki/Equation_of_time*/
                double TC = 4*(longitude - LSTM) + EoT; //time correction factor including EoT correction for Earth eccentricity
                double LST = lt - TC/60; // local solar time
                return 15*(lt-12)*deg2rad; ///TODO: find right EOT and data for validation, so it should return: 15*(LST-12)
            }
            /**\brief computes earth declination angle
             * // ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, eq.(D.5)
             * \param doy -- day of the year
             * \return declination, [rad] */
            double  compute_earth_declination(double doy){
                double G = 2*pi/365.0*(doy-1);
                return  0.006918 - 0.399912*cos(G)+0.070257*sin(G)- 0.006758*cos(2*G)+0.000907*sin(2*G) - 0.002697*cos(3*G)+0.00148*sin(3*G);
            }
            /**\brief computes extraterrestrial solar radiation
             * \param cos(theta)
             * ]param doy -- day of the year
             * return ra, [W/m^2]*/
            double compute_ra(double cos_theta, double doy ){
                return gsc*cos_theta*(1+0.0033*cos(doy*2*pi/365)); // eq.(1)
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

    TEST_CASE("check_hour_angle"){
        /** \brief Check Earth Hour Angle
         * based on ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, no EOT correction provided*/
        parameter p;
        calculator r(p);
        calendar utc_cal;
        double lat = 56.0;
        utctime t1, t2, t3;
        arma::vec surface_normal({0.0,0.0,1.0});
        t1 = utc_cal.time(1970, 1, 1, 10, 30, 0, 0);
        t2 = utc_cal.time(1970, 6, 1, 12, 0, 0, 0);
        t3 = utc_cal.time(1970, 9, 1, 23, 30, 0, 0);
        // using simple formulae for checking, need to provide test for EOT
        r.rso_cs_radiation(lat, t1, surface_normal);
                FAST_CHECK_EQ(r.hra(), doctest::Approx(-22.5));// earth hour angle
        r.rso_cs_radiation(lat, t2, surface_normal);
                FAST_CHECK_EQ(r.hra(), doctest::Approx(0.0));// earth hour angle
        r.rso_cs_radiation(lat, t3, surface_normal);
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
        arma::vec surface_normal({0.0,0.0,1.0});
        t1 = utc_cal.time(1970, 12, 21, 12, 30, 0, 0);
        t2 = utc_cal.time(1970, 6, 21, 12, 30, 0, 0);
        t3 = utc_cal.time(1970, 3, 21, 12, 30, 0, 0);
        t4 = utc_cal.time(1970, 9, 24, 12, 30, 0, 0);
        r.rso_cs_radiation(lat, t1, surface_normal);
                FAST_CHECK_EQ(r.declination(), doctest::Approx(-23.5).epsilon(0.1));// earth declination angle 21/22 December (min)
        r.rso_cs_radiation(lat, t2, surface_normal);
                FAST_CHECK_EQ(r.declination(), doctest::Approx(23.5).epsilon(0.1));// earth declination angle 21/22 June (max)
        r.rso_cs_radiation(lat, t3, surface_normal);
                FAST_CHECK_EQ(r.declination(), doctest::Approx(0.0).epsilon(0.1));// earth declination angle 20/21 March 0
        r.rso_cs_radiation(lat, t4, surface_normal);
                FAST_CHECK_EQ(r.declination(), doctest::Approx(0.0).epsilon(0.5));// earth declination angle 22/23 September 0 this gives a bit higher error
    }
    TEST_CASE("check_atm_pressures"){
        parameter p;
        calculator r(p);
        calendar utc_cal;
        double lat = 56.0;
        utctime t;
        arma::vec surface_normal({0.0,0.0,1.0});
        t = utc_cal.time(1970, 12, 21, 12, 30, 0, 0);
        r.rso_cs_radiation(lat, t, surface_normal, 0.0, 0.0, 0.0);
                FAST_CHECK_EQ(r.patm(), doctest::Approx(101.325).epsilon(0.1));
        r.rso_cs_radiation(lat, t, surface_normal, 0.0, 0.0, 1000.0);
                FAST_CHECK_EQ(r.patm(), doctest::Approx(89.88).epsilon(0.1)); // https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118568101.app2
        r.rso_cs_radiation(lat, t, surface_normal, 0.0, 0.0,  -1000.0);
                FAST_CHECK_EQ(r.patm(), doctest::Approx(113.93).epsilon(0.1));

    }
    TEST_CASE("check_vap_pressure"){
        /**\brief check saturation pressure (rhumidity=100)
         * based ref.:  Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.113, fig.3.1*/
        parameter p;
        calculator r(p);
        calendar utc_cal;
        double lat = 56.0;
        utctime t;
        arma::vec surface_normal({0.0,0.0,1.0});
        t = utc_cal.time(1970, 12, 21, 12, 30, 0, 0);
        r.rso_cs_radiation(lat, t, surface_normal, 21.0, 100.0, 0.0);
                FAST_CHECK_EQ(r.ea(), doctest::Approx(2.5).epsilon(0.01));
        r.rso_cs_radiation(lat, t, surface_normal, 21.0, 50.0, 0.0);
                FAST_CHECK_EQ(r.ea(), doctest::Approx(1.25).epsilon(0.01));
        r.rso_cs_radiation(lat, t, surface_normal, 0.0, 100.0, 0.0);
                FAST_CHECK_EQ(r.ea(), doctest::Approx(0.6).epsilon(0.01));
        r.rso_cs_radiation(lat, t, surface_normal, 27.0, 100.0, 0.0);
                FAST_CHECK_EQ(r.ea(), doctest::Approx(3.6).epsilon(0.01));
        r.rso_cs_radiation(lat, t, surface_normal, -31.0, 100.0, 0.0);
                FAST_CHECK_EQ(r.ea(), doctest::Approx(0.03).epsilon(0.01));
    }
    TEST_CASE("check_slope_aspect"){
        /**\brief check slope and aspect*/
        parameter p;
        calculator r(p);
        calendar utc_cal;
        double lat = 56.0;
        utctime t;
        arma::vec surface_normal({0.0,0.0,1.0});
        t = utc_cal.time(1970, 12, 21, 12, 30, 0, 0);
        r.rso_cs_radiation(lat, t, surface_normal, -31.0, 100.0, 0.0);
                FAST_CHECK_EQ(r.slope(), doctest::Approx(0.0).epsilon(0.01));
                FAST_CHECK_EQ(r.aspect(), doctest::Approx(180.0).epsilon(0.01));
    }

    TEST_CASE("check_sun_rise_and_set"){

        parameter p;
        p.albedo=0.2;
        calculator r(p);
        calendar utc_cal;
        double lat = 44.0;
        utctime t;
        arma::vec surface_normal({0.5,0.5,1.0});
        // January
        t = utc_cal.time(1970, 01, 1, 12, 30, 0, 0);
        r.rso_cs_radiation(lat, t, surface_normal, -31.0, 100.0, 150.0);
                FAST_CHECK_EQ(r.sun_rise(), doctest::Approx(10.0).epsilon(0.01));
                FAST_CHECK_EQ(r.sun_set(), doctest::Approx(23.0).epsilon(0.01));
        // March
        t = utc_cal.time(1970, 03, 1, 12, 30, 0, 0);
        r.rso_cs_radiation(lat, t, surface_normal, -31.0, 100.0, 150.0);
                FAST_CHECK_EQ(r.sun_rise(), doctest::Approx(10.0).epsilon(0.01));
                FAST_CHECK_EQ(r.sun_set(), doctest::Approx(23.0).epsilon(0.01));
        // July
        t = utc_cal.time(1970, 07, 16, 12, 30, 0, 0);
        r.rso_cs_radiation(lat, t, surface_normal, -31.0, 100.0, 150.0);
                FAST_CHECK_EQ(r.sun_rise(), doctest::Approx(10.0).epsilon(0.01));
                FAST_CHECK_EQ(r.sun_set(), doctest::Approx(23.0).epsilon(0.01));
                // December
        t = utc_cal.time(1970, 11, 21, 12, 30, 0, 0);
        r.rso_cs_radiation(lat, t, surface_normal, -31.0, 100.0, 150.0);
                FAST_CHECK_EQ(r.sun_rise(), doctest::Approx(10.0).epsilon(0.01));
                FAST_CHECK_EQ(r.sun_set(), doctest::Approx(23.0).epsilon(0.01));


    }
    TEST_CASE("check_solar_radiation"){
        parameter p;
        p.albedo=0.2;
        calculator r(p);
        calendar utc_cal;
        double lat = 44.0;
        utctime t;
        arma::vec surface_normal({0.5,0.5,1.0});
        t = utc_cal.time(1970, 06, 21, 12, 30, 0, 0); // June 12.30
        r.rso_cs_radiation(lat, t, surface_normal, -31.0, 100.0, 150.0);
                FAST_CHECK_EQ(r.cs_radiation(), doctest::Approx(500.0).epsilon(0.01));
        t = utc_cal.time(1970, 06, 21, 23, 30, 0, 0); // June 23.30 -- no radiation
        r.rso_cs_radiation(lat, t, surface_normal, -31.0, 100.0, 150.0);
                FAST_CHECK_EQ(r.cs_radiation(), doctest::Approx(0.0).epsilon(0.01));
        t = utc_cal.time(1970, 06, 21, 10, 30, 0, 0); // June 10.30
        r.rso_cs_radiation(lat, t, surface_normal, -31.0, 100.0, 150.0);
                FAST_CHECK_EQ(r.cs_radiation(), doctest::Approx(300.0).epsilon(0.01));
        t = utc_cal.time(1970, 11, 21, 12, 30, 0, 0); // December 12.30
        r.rso_cs_radiation(lat, t, surface_normal, -31.0, 100.0, 150.0);
                FAST_CHECK_EQ(r.cs_radiation(), doctest::Approx(150.0).epsilon(0.01));
        t = utc_cal.time(1970, 11, 21, 23, 30, 0, 0); // December 23.30
        r.rso_cs_radiation(lat, t, surface_normal, -31.0, 100.0, 150.0);
                FAST_CHECK_EQ(r.cs_radiation(), doctest::Approx(0.0).epsilon(0.01));
        t = utc_cal.time(1970, 11, 21, 10, 30, 0, 0); // December 10.30
        r.rso_cs_radiation(lat, t, surface_normal, -31.0, 100.0, 150.0);
                FAST_CHECK_EQ(r.cs_radiation(), doctest::Approx(50.0).epsilon(0.01));

    }
    TEST_CASE("surface_normal_from_cells") {
        vector<cell> cells;
        auto r= surface_normal(cells);
    }
}
