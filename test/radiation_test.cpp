#include "test_pch.h"
#include "core/utctime_utilities.h"
#include "core/geo_cell_data.h"
#include <armadillo>
#include <vector>
#include <boost/math/constants/constants.hpp>

namespace shyft::core {
    using std::vector;
    using std::cos;
    using std::sin;
    using std::pow;
    using std::exp;
    const double pi = boost::math::constants::pi<double>();

    namespace radiation {

        struct calculator {
            calendar utc;
            double doy; // day of the year
            double hour;
            double omega; // earth hour angle

            double s = 0.0; // horizontal slope, should be calculated from normal vector
            double gamma = 0.0; // surface aspect angle, 0 -- due south, -pi/2 -- due east, pi/2 -- due west, +-pi -- due north

            double delta = 0.1; // declination of the earth, positive for northern hemisphere summer

            double phi = 0.0 ;// latitude, should be available from cell?

            double ra = 0.0; // extraterrestrial solar radiation for inclined surface[W/m2]
            double rahor = 0.0; // extraterrestrial solar radiation for horizontal surfaces
            double rso = 0.0; // clear sky solar radiation for inclined surface [W/m2]

            double alpha = 0.1; // albedo

            /** \tparam V is typically arma::vec */
            template <class V>
            double real_radiation(double latitude, utctime t, const V & surface_normal) {
                doy=utc.day_of_year(t);
                hour = utc.calendar_units(t).hour;
                return hour;
            }
            const double gsc = 1367; // W/m2 -- solar constant

            /** \brief computes instantaneous extraterrestrial solar radiation
             * \param latitude (phi)
             * \param time
             * \param surface normal */

            template <class V>
            double ra_radiation(double latitude, utctime t, const V & surface_normal){
                phi = latitude;
                delta = 23.45*pi/180*sin(2*pi*(284+doy)/36.25); // earth declination angle based on eq.(3.1) https://www.itacanet.org/the-sun-as-a-source-of-energy/part-3-calculating-solar-angles/
                omega = utc.calendar_units(t).hour; // earth  hour angle /// TODO convert time to hour angle (HRA)
                compute_costt(omega);
                ra = gsc*cos_tt*(1+0.0033*cos(doy*2*pi/365)); // eq.(1)
                rahor = gsc*cos_tthor*(1+0.0033*cos(doy*2*pi/365)); // eq.(1) with cos(theta)hor
                // std::cout<<ra<<std::endl;
                return ra;
            }
            /** \brief computes instantaneous clear-sky radiation  for inclined surfaces*/
            template <class C>
            double rso_cs_radiation(utctime t, const C& cell){
                double longitude = 0.0; // longitude should come from cell? ///TODO get from cell, ask if we have it somewhere, or need to calculate from coordinates, where is wgs84?
                omega = hour_angle(longitude,utc.calendar_units(t).hour); // hour angle
                double tau_swo; // broadband atmospheric transmissivity for shortwave radiation for cloud-free conditions
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
                rso = Kbo*ra + (fia*Kdo + alpha*(1-fi)*(Kbo+Kdo))*rahor; // eq.(37)
                return rso; // eq.(37)
            }
        private:
            double a, b, c, cos_tt, cos_tthor;
            /** \brief computes necessary geometric parameters
             * \param omega -- hour angle
             * */
            void compute_costt(double omega){
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
             * ref.: https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-time
             * \param longitute
             * \param lt - local time [h]*/
            double hour_angle(double longitude, double lt){
                double DTutc = 0.0; // [hour], difference of LT from UTC///TODO get timezone from utctime utilities
                double LSTM = 15*DTutc; // local standard time meridian
                double B = 360/365*(doy-81);
                double EoT = 9.87*sin(2*B) - 7.3*cos(B) - 1.5*sin(B);// equation of time
                double TC = 4*(longitude - LSTM) + EoT; //time correction factor
                double LST = lt + TC/60; // loacl solar time
                return 15*(LST-12);
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
    // test basics: creation, etc
    
    TEST_CASE("compute_real_radiation") {
        calculator r;
        FAST_CHECK_EQ(r.real_radiation(1.0,calendar::HOUR*25,1),doctest::Approx(1));
        FAST_CHECK_EQ(r.ra_radiation(1.0,calendar::HOUR*25,1),doctest::Approx(1)); // should fail
    }
    
    TEST_CASE("surface_normal_from_cells") {
        vector<cell> cells;
        auto r= surface_normal(cells);
    }
}
