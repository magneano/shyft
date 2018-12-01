#include "test_pch.h"
#include "core/utctime_utilities.h"
#include <armadillo>
#include <vector>
#include <chrono>
#include <boost/math/constants/constants.hpp>
#include "core/radiation.h"
#include <cmath>
#include <random>

//namespace shyft::core::radiation{
//    /** \tparam C is a cell, like shyft-cell, */
//    template <class C>
//    vector<arma::vec> surface_normal( const vector<C>& cells){
//        vector<arma::vec> r;
//        for(const auto&c:cells) {
//            double x=c.geo.mid_point().x;
//            r.push_back(arma::vec({1.0,1.0,1.0}));
//        }
//        return r;
//    }
//}

namespace shyft::test {

    class trapezoidal_average {
    private:
        double area = 0.0;
        double f_a = 0.0;; // Left hand side of next integration subinterval
        double t_start = 0.0; // Start of integration period
        double t_a = 0.0; // Left hand side time of next integration subinterval
    public:
        explicit trapezoidal_average() {}

        /** \brief initialize must be called to reset states before being used during ode integration.
         */
        void initialize(double f0, double t_start) {
            this->f_a = f0;
            this->t_start = t_start;
            t_a = t_start;
            area = 0.0;
        }

        /** \brief Add contribution to average using a simple trapezoidal rule
         *
         * See: http://en.wikipedia.org/wiki/Numerical_integration
         */
        void add(double f, double t) {
            area += 0.5*(f_a + f)*(t - t_a);
            f_a = f;
            t_a = t;
        }

        double result() const { return area/(t_a - t_start); }
    };

}

TEST_SUITE("radiation") {
    using shyft::core::radiation::parameter;
    using shyft::core::radiation::response;
    using shyft::core::radiation::calculator;
//    using shyft::core::radiation::surface_normal;
    using shyft::core::calendar;
    using shyft::core::utctime;
    using shyft::test::trapezoidal_average;
    // test basics: creation, etc


//    TEST_CASE("check_slope_aspect"){
//        /**\brief check slope and aspect*/
//        parameter p;
//        calculator r(p);
//        calendar utc_cal;
//        double lat = 56.0;
//        utctime t;
//        arma::vec surface_normal({0.0,0.0,1.0});
//        t = utc_cal.time(1970, 12, 21, 12, 30, 0, 0);
//        r.psw_radiation(lat, t, surface_normal, -31.0, 100.0, 0.0);
//                FAST_CHECK_EQ(r.slope(), doctest::Approx(0.0).epsilon(0.01));
//                FAST_CHECK_EQ(r.aspect(), doctest::Approx(0.0).epsilon(0.01));
//    }

    TEST_CASE("check_solar_radiation_horizontal"){
        parameter p;
        response r;
        p.albedo = 0.2;
        p.turbidity = 1.0;
        calculator<parameter,response> rad(p);
        calendar utc_cal;
        double lat = 44.0;
        utctime t;
        // checking for horizontal surface Eugene, OR, p.64, fig.1b
        arma::vec surface_normal({0.0,0.0,1.0});
        utctime ta;
        trapezoidal_average av_rahor;
        trapezoidal_average av_ra;
        trapezoidal_average av_rso;
        trapezoidal_average av_rs;
        std::uniform_real_distribution<double> ur(150.0, 390.0);
        std::default_random_engine gen;
        std::cout << "========= Horizontal =======" << std::endl;
        SUBCASE("June") {
            std::cout << "========= June ========" << std::endl;
            ta = utc_cal.time(1970, 06, 21, 00, 00, 0, 0);
            //rad.psw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.psw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 06, 21, h, 00, 0, 0); // June
                //rad.psw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0);
                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.psw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }

            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(500.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rahor.result(), doctest::Approx(av_ra.result()).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(370.0).epsilon(0.05));

        }
        SUBCASE("January") {
            std::cout << "========= January =======" << std::endl;
            ta = utc_cal.time(1970, 01, 1, 00, 00, 0, 0);
            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.psw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 01, 1, h, 00, 0, 0); // January
                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.psw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }

            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(130.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rahor.result(), doctest::Approx(av_ra.result()).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(80.0).epsilon(0.05));
        }
        SUBCASE("December") {
            std::cout << "========= December =======" << std::endl;
            ta = utc_cal.time(1970, 12, 21, 00, 00, 0, 0);
            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.tsw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 12, 21, h, 00, 0, 0); // January
                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.psw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }

            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(130.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rahor.result(), doctest::Approx(av_ra.result()).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(80.0).epsilon(0.05));
        }

    }
   TEST_CASE("check_solar_radiation_slope_45s"){
        parameter p;
        response r;
        p.albedo = 0.2;
        p.turbidity = 1.0;
        calculator<parameter,response> rad(p);
        calendar utc_cal;
        double lat = 44.0;
        utctime t;
        // checking for horizontal surface Eugene, OR, p.64, fig.1d
        // 24h  average radiation
        double slope = 45*shyft::core::radiation::pi/180; // 45 S
        double proj = sin(slope);
        double aspect = 0.0*shyft::core::radiation::pi/180;// facing south
        arma::vec surface_normal({proj*cos(aspect),proj*sin(aspect),cos(slope)});
        utctime ta;
        trapezoidal_average av_rahor;
        trapezoidal_average av_ra;
        trapezoidal_average av_rso;
        trapezoidal_average av_rs;
        std::uniform_real_distribution<double> ur(100.0, 390.0);
        std::default_random_engine gen;
        std::cout << "========= Slope 45S =======" << std::endl;
        SUBCASE("June") {
            std::cout << "========= June ========" << std::endl;
            ta = utc_cal.time(1970, 06, 21, 00, 00, 0, 0);
            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.psw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 06, 21, h, 00, 0, 0); // June
                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.psw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }

            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
            std::cout << "sun_rise: " << rad.sun_rise() << std::endl;
            std::cout << "sun_set: " << rad.sun_set() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(390.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(310.0).epsilon(0.05));
        }
        SUBCASE("January") {
            std::cout << "========= January ========" << std::endl;
            ta = utc_cal.time(1970, 01, 1, 00, 00, 0, 0);
            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.psw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 01, 1, h, 00, 0, 0); // June
                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.psw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
            std::cout << "sun_rise: " << rad.sun_rise() << std::endl;
            std::cout << "sun_set: " << rad.sun_set() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(390.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(200.0).epsilon(0.05));
        }
        SUBCASE("December") {
            std::cout << "========= December ========" << std::endl;
            ta = utc_cal.time(1970, 12, 12, 00, 00, 0, 0);
            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.psw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 12, 12, h, 00, 0, 0); // June
                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.psw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(390.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(200.0).epsilon(0.05));
        }

    }
    TEST_CASE("check_solar_radiation_slope_90S"){
        parameter p;
        response r;
        p.albedo = 0.05;
        p.turbidity = 1.0;
        calculator<parameter,response> rad(p);
        calendar utc_cal;
        double lat = 44.0;
        utctime t;
        // checking for horizontal surface Eugene, OR, p.64, fig.1e
        double slope = 90*shyft::core::radiation::pi/180; // 90 S
        double proj = sin(slope);
        double aspect = 0.0;// facing south
        arma::vec surface_normal({proj*cos(aspect),proj*sin(aspect),cos(slope)});

        utctime ta;
        trapezoidal_average av_rahor;
        trapezoidal_average av_ra;
        trapezoidal_average av_rso;
        trapezoidal_average av_rs;
        std::uniform_real_distribution<double> ur(0.0, 250.0);
        std::default_random_engine gen;
        std::cout << "========= Slope 90S =======" << std::endl;
        SUBCASE("June") {
            std::cout << "========= June ========" << std::endl;
            ta = utc_cal.time(1970, 06, 21, 00, 00, 0, 0);
            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.psw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 06, 21, h, 00, 0, 0); // June
                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.psw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(110.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(85.0).epsilon(0.05));
        }
        SUBCASE("January") {
            std::cout << "========= January ========" << std::endl;
            ta = utc_cal.time(1970, 01, 1, 00, 00, 0, 0);
            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.psw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 01, 1, h, 00, 0, 0); // June
                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.psw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(410.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(200.0).epsilon(0.05));
        }
        SUBCASE("December") {
            std::cout << "========= December ========" << std::endl;
            ta = utc_cal.time(1970, 12, 12, 00, 00, 0, 0);
            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.psw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 12, 12, h, 00, 0, 0); // June
                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.psw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(410.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(200.0).epsilon(0.05));
        }

    }
    TEST_CASE("check_solar_radiation_slope_90N"){
        parameter p;
        response r;
        p.albedo = 0.05;
        p.turbidity = 1.0;
        calculator<parameter,response> rad(p);
        calendar utc_cal;
        double lat = 44.0;
        utctime t;
        // checking for horizontal surface Eugene, OR, p.64, fig.1f
        double slope = 90*shyft::core::radiation::pi/180; // 90 N
        double proj = sin(slope);
        double aspect = 180*shyft::core::radiation::pi/180;// facing north
        arma::vec surface_normal({proj*cos(aspect),proj*sin(aspect),cos(slope)});
        double ra_sum = 0.0;
        double rso_sum = 0.0;
        double rahor_sum = 0.0;
        utctime ta;
        trapezoidal_average av_rahor;
        trapezoidal_average av_ra;
        trapezoidal_average av_rso;
        trapezoidal_average av_rs;
        std::uniform_real_distribution<double> ur(10.0, 100.0);
        std::default_random_engine gen;
        std::cout << "========= Slope 90N =======" << std::endl;
        SUBCASE("June") {
            std::cout << "========= June ========" << std::endl;
            ta = utc_cal.time(1970, 06, 21, 00, 00, 0, 0);
            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.psw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 06, 21, h, 00, 0, 0); // June
                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.psw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(105.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(50.0).epsilon(0.05));
        }
        SUBCASE("January") {
            std::cout << "========= January ========" << std::endl;
            ta = utc_cal.time(1970, 01, 1, 00, 00, 0, 0);
            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.psw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 01, 1, h, 00, 0, 0); // June
                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.psw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(0.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(9.0).epsilon(0.05));
        }
        SUBCASE("December") {
            std::cout << "========= December ========" << std::endl;
            ta = utc_cal.time(1970, 12, 12, 00, 00, 0, 0);
            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.psw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 12, 12, h, 00, 0, 0); // June
                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.psw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(0.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(9.0).epsilon(0.05));
        }

    }
    /*TEST_CASE("check_solar_radiation_slope_90N_golden"){
        parameter p;
        response res;
        p.albedo = 0.05;
        p.turbidity = 1.0;
        calculator<parameter,response> r(p);
        calendar utc_cal;
        double lat = 39.74;
        double elevation = 1829;
        double temperature = 0.0;
        double rhumidity = 50.0;
        utctime t;
        // checking for horizontal surface Eugene, OR, p.64, fig.1f
        double slope = 90*shyft::core::radiation::pi/180; // 90 N
        double proj = sin(slope);
        double aspect = 180*shyft::core::radiation::pi/180;// facing north
        arma::vec surface_normal({proj*cos(aspect),proj*sin(aspect),cos(slope)});
        double ra_sum = 0.0;
        double rso_sum = 0.0;
        double rahor_sum = 0.0;
        utctime ta;
        trapezoidal_average av_rahor;
        trapezoidal_average av_ra;
        trapezoidal_average av_rso;
        trapezoidal_average av_rs;
        std::uniform_real_distribution<double> ur(10.0, 100.0);
        std::default_random_engine gen;
        std::cout << "========= Slope 90N golden =======" << std::endl;
                SUBCASE("June") {
            std::cout << "========= June ========" << std::endl;
            ta = utc_cal.time(1970, 06, 21, 00, 00, 0, 0);
            rad.psw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, elevation, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.tsw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 06, 21, h, 00, 0, 0); // June
                rad.psw_radiation(r, lat, t, surface_normal, 20.0, 50.0, elevation, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.tsw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(100.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(90.0).epsilon(0.05));
        }
                SUBCASE("January") {
            std::cout << "========= January ========" << std::endl;
            ta = utc_cal.time(1970, 01, 1, 00, 00, 0, 0);
            rad.psw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, elevation, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.tsw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 01, 1, h, 00, 0, 0); // June
                rad.psw_radiation(r, lat, t, surface_normal, 20.0, 50.0, elevation, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.tsw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(0.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(19.0).epsilon(0.05));
        }
                SUBCASE("December") {
            std::cout << "========= December ========" << std::endl;
            ta = utc_cal.time(1970, 12, 12, 00, 00, 0, 0);
            rad.psw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, elevation, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rso.initialize(r.tsw_radiation, 0.0);
            av_rs.initialize(r.tsw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(1970, 12, 12, h, 00, 0, 0); // June
                rad.psw_radiation(r, lat, t, surface_normal, 20.0, 50.0, elevation, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rso.add(r.tsw_radiation, h);
                av_rs.add(r.tsw_radiation, h);
            }
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rso: " << av_rso.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(0.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(19.0).epsilon(0.05));
        }

    }*/
//    TEST_CASE("surface_normal_from_cells") {
//        vector<cell> cells;
//        auto r= surface_normal(cells);
//    }
}
