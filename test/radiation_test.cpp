#include "test_pch.h"
#include "core/utctime_utilities.h"
#include <armadillo>
#include <vector>
#include <chrono>
#include <boost/math/constants/constants.hpp>
#include "core/radiation.h"
#include <cmath>
#include <random>
#include <tuple>

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
namespace rasputin {
    using namespace std;
    using Point = std::array<double, 3>;
    using PointList = std::vector<Point>;
    using VectorList = PointList;
    using Vector = Point;
    using Face = std::array<int, 3>;
    using FaceList = std::vector<Face>;

    VectorList surface_normals(const PointList &pts, const FaceList &faces) {
        VectorList result;
        result.reserve(faces.size());
        for (const auto face: faces) {
            const auto p0 = pts[face[0]];
            const auto p1 = pts[face[1]];
            const auto p2 = pts[face[2]];
            const arma::vec::fixed<3> v0 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
            const arma::vec::fixed<3> v1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
            const arma::vec::fixed<3> n = arma::cross(v0/arma::norm(v0), v1/arma::norm(v1));
            result.emplace_back(n.at(2) >= 0.0 ? Vector{n.at(0), n.at(1), n.at(2)} : Vector{-n.at(0), -n.at(1), -n.at(2)});
        }
        return result;
    }
    std::vector<double> slopes(const PointList &pts, const FaceList &faces){
        VectorList normals = surface_normals(pts,faces);
        std::vector<double> result;
        result.reserve(normals.size());
        for (const auto &normal: normals) {
            result.push_back(atan2(pow(pow(normal.at(0),2)+pow(normal.at(1),2),0.5),normal.at(2)));
        }
        return result;
    }
    std::vector<double> aspects(const PointList &pts, const FaceList &faces){
        VectorList normals = surface_normals(pts,faces);
        std::vector<double> result;
        result.reserve(normals.size());
        for (const auto &normal: normals) {
            result.push_back(atan2(normal.at(1),normal.at(0)));
        }
        return result;
    }

    double elevation{0.0};
}

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



    TEST_CASE("geometry"){

        std::array<double,3> point1{{0.0,0.0,0.0}};
        std::array<double,3> point2{{1.0,0.0,0.0}};
        std::array<double,3> point3{{0.0,1.0,0.0}};
        std::array<int,3> face{{0,1,2}};
        std::vector<std::array<double,3>> points;
        points.push_back(point1);
        points.push_back(point2);
        points.push_back(point3);

        std::vector<std::array<int,3>> faces;
        faces.push_back(face);

        std::vector<std::array<double,3>> normals;
        normals = rasputin::surface_normals(points,faces);
        for (const auto n:normals)
            std::cout<<n.at(0)<<"; "<<n.at(1) << "; "<<n.at(2)<<std::endl;

        std::vector<double> slopes;

        slopes = rasputin::slopes(points,faces);
        for (const auto s:slopes)
            std::cout<<s<<std::endl;
        //FAST_CHECK_EQ(tr.slope, doctest::Approx(0.0).epsilon(0.01));
//        FAST_CHECK_EQ(tr.aspect, doctest::Approx(0.0).epsilon(0.01));




    }
//    TEST_CASE("check_slope_aspect"){
//        /**\brief check slope and aspect*/
//        parameter p;
//        response r;
//        p.albedo = 0.2;
//        p.turbidity = 1.0;
//        calculator<parameter,response> rad(p);
//        calendar utc_cal;
//        double lat = 44.0;
//        utctime t;
//        // checking for horizontal surface Eugene, OR, p.64, fig.1b
//        arma::vec surface_normal({0.0,0.0,1.0});
//        t = utc_cal.time(1970, 12, 21, 12, 30, 0, 0);
//        rad.psw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0);
//        FAST_CHECK_EQ(rad.slope(), doctest::Approx(0.0).epsilon(0.01));
//        FAST_CHECK_EQ(rad.aspect(), doctest::Approx(0.0).epsilon(0.01));
//    }

//    TEST_CASE("check_solar_radiation_horizontal"){
//        parameter p;
//        response r;
//        p.albedo = 0.2;
//        p.turbidity = 1.0;
//        calculator<parameter,response> rad(p);
//        calendar utc_cal;
//        double lat = 44.0;
//        utctime t;
//        // checking for horizontal surface Eugene, OR, p.64, fig.1b
//        arma::vec surface_normal({0.0,0.0,1.0});
//        utctime ta;
//        trapezoidal_average av_rahor;
//        trapezoidal_average av_ra;
//        trapezoidal_average av_rso;
//        trapezoidal_average av_rs;
//        std::uniform_real_distribution<double> ur(150.0, 390.0);
//        std::default_random_engine gen;
//        std::cout << "========= Horizontal =======" << std::endl;
//        SUBCASE("June") {
//            std::cout << "========= June ========" << std::endl;
//            ta = utc_cal.time(1970, 06, 21, 00, 00, 0, 0);
//            //rad.psw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
//            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
//            av_ra.initialize(rad.ra_radiation(), 0.0);
//            av_rso.initialize(r.psw_radiation, 0.0);
//            av_rs.initialize(r.tsw_radiation, 0.0);
//            for (int h = 1; h < 24; ++h) {
//                t = utc_cal.time(1970, 06, 21, h, 00, 0, 0); // June
//                //rad.psw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0);
//                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//                av_rahor.add(rad.ra_radiation_hor(), h);
//                av_ra.add(rad.ra_radiation(), h);
//                av_rso.add(r.psw_radiation, h);
//                av_rs.add(r.tsw_radiation, h);
//            }
//
//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rso: " << av_rso.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
//                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(500.0).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rahor.result(), doctest::Approx(av_ra.result()).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(370.0).epsilon(0.05));
//
//        }
//        SUBCASE("January") {
//            std::cout << "========= January =======" << std::endl;
//            ta = utc_cal.time(1970, 01, 1, 00, 00, 0, 0);
//            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
//            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
//            av_ra.initialize(rad.ra_radiation(), 0.0);
//            av_rso.initialize(r.psw_radiation, 0.0);
//            av_rs.initialize(r.tsw_radiation, 0.0);
//            for (int h = 1; h < 24; ++h) {
//                t = utc_cal.time(1970, 01, 1, h, 00, 0, 0); // January
//                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0);
//                av_rahor.add(rad.ra_radiation_hor(), h);
//                av_ra.add(rad.ra_radiation(), h);
//                av_rso.add(r.psw_radiation, h);
//                av_rs.add(r.tsw_radiation, h);
//            }
//
//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rso: " << av_rso.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
//                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(130.0).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rahor.result(), doctest::Approx(av_ra.result()).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(80.0).epsilon(0.05));
//        }
//        SUBCASE("December") {
//            std::cout << "========= December =======" << std::endl;
//            ta = utc_cal.time(1970, 12, 21, 00, 00, 0, 0);
//            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
//            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
//            av_ra.initialize(rad.ra_radiation(), 0.0);
//            av_rso.initialize(r.tsw_radiation, 0.0);
//            av_rs.initialize(r.tsw_radiation, 0.0);
//            for (int h = 1; h < 24; ++h) {
//                t = utc_cal.time(1970, 12, 21, h, 00, 0, 0); // January
//                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0);
//                av_rahor.add(rad.ra_radiation_hor(), h);
//                av_ra.add(rad.ra_radiation(), h);
//                av_rso.add(r.psw_radiation, h);
//                av_rs.add(r.tsw_radiation, h);
//            }
//
//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rso: " << av_rso.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
//                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(130.0).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rahor.result(), doctest::Approx(av_ra.result()).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(80.0).epsilon(0.05));
//        }
//
//    }
//   TEST_CASE("check_solar_radiation_slope_45s"){
//        parameter p;
//        response r;
//        p.albedo = 0.2;
//        p.turbidity = 1.0;
//        calculator<parameter,response> rad(p);
//        calendar utc_cal;
//        double lat = 44.0;
//        utctime t;
//        // checking for horizontal surface Eugene, OR, p.64, fig.1d
//        // 24h  average radiation
//        double slope = 45*shyft::core::radiation::pi/180; // 45 S
//        double proj = sin(slope);
//        double aspect = 0.0*shyft::core::radiation::pi/180;// facing south
//        arma::vec surface_normal({proj*cos(aspect),proj*sin(aspect),cos(slope)});
//        utctime ta;
//        trapezoidal_average av_rahor;
//        trapezoidal_average av_ra;
//        trapezoidal_average av_rso;
//        trapezoidal_average av_rs;
//        std::uniform_real_distribution<double> ur(100.0, 390.0);
//        std::default_random_engine gen;
//        std::cout << "========= Slope 45S =======" << std::endl;
//        SUBCASE("June") {
//            std::cout << "========= June ========" << std::endl;
//            ta = utc_cal.time(1970, 06, 21, 00, 00, 0, 0);
//            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
//            av_ra.initialize(rad.ra_radiation(), 0.0);
//            av_rso.initialize(r.psw_radiation, 0.0);
//            av_rs.initialize(r.tsw_radiation, 0.0);
//            for (int h = 1; h < 24; ++h) {
//                t = utc_cal.time(1970, 06, 21, h, 00, 0, 0); // June
//                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//                av_rahor.add(rad.ra_radiation_hor(), h);
//                av_ra.add(rad.ra_radiation(), h);
//                av_rso.add(r.psw_radiation, h);
//                av_rs.add(r.tsw_radiation, h);
//            }
//
//            std::cout << "rahor: " << av_rahor.result() << std::endl;
//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rso: " << av_rso.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
//            std::cout << "sun_rise: " << rad.sun_rise() << std::endl;
//            std::cout << "sun_set: " << rad.sun_set() << std::endl;
//                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(390.0).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(310.0).epsilon(0.05));
//        }
//        SUBCASE("January") {
//            std::cout << "========= January ========" << std::endl;
//            ta = utc_cal.time(1970, 01, 1, 00, 00, 0, 0);
//            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
//            av_ra.initialize(rad.ra_radiation(), 0.0);
//            av_rso.initialize(r.psw_radiation, 0.0);
//            av_rs.initialize(r.tsw_radiation, 0.0);
//            for (int h = 1; h < 24; ++h) {
//                t = utc_cal.time(1970, 01, 1, h, 00, 0, 0); // June
//                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//                av_rahor.add(rad.ra_radiation_hor(), h);
//                av_ra.add(rad.ra_radiation(), h);
//                av_rso.add(r.psw_radiation, h);
//                av_rs.add(r.tsw_radiation, h);
//            }
//            std::cout << "rahor: " << av_rahor.result() << std::endl;
//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rso: " << av_rso.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
//            std::cout << "sun_rise: " << rad.sun_rise() << std::endl;
//            std::cout << "sun_set: " << rad.sun_set() << std::endl;
//                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(390.0).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(200.0).epsilon(0.05));
//        }
//        SUBCASE("December") {
//            std::cout << "========= December ========" << std::endl;
//            ta = utc_cal.time(1970, 12, 12, 00, 00, 0, 0);
//            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
//            av_ra.initialize(rad.ra_radiation(), 0.0);
//            av_rso.initialize(r.psw_radiation, 0.0);
//            av_rs.initialize(r.tsw_radiation, 0.0);
//            for (int h = 1; h < 24; ++h) {
//                t = utc_cal.time(1970, 12, 12, h, 00, 0, 0); // June
//                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//                av_rahor.add(rad.ra_radiation_hor(), h);
//                av_ra.add(rad.ra_radiation(), h);
//                av_rso.add(r.psw_radiation, h);
//                av_rs.add(r.tsw_radiation, h);
//            }
//            std::cout << "rahor: " << av_rahor.result() << std::endl;
//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rso: " << av_rso.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
//                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(390.0).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(200.0).epsilon(0.05));
//        }
//
//    }
//    TEST_CASE("check_solar_radiation_slope_90S"){
//        parameter p;
//        response r;
//        p.albedo = 0.05;
//        p.turbidity = 1.0;
//        calculator<parameter,response> rad(p);
//        calendar utc_cal;
//        double lat = 44.0;
//        utctime t;
//        // checking for horizontal surface Eugene, OR, p.64, fig.1e
//        double slope = 90*shyft::core::radiation::pi/180; // 90 S
//        double proj = sin(slope);
//        double aspect = 0.0;// facing south
//        arma::vec surface_normal({proj*cos(aspect),proj*sin(aspect),cos(slope)});
//
//        utctime ta;
//        trapezoidal_average av_rahor;
//        trapezoidal_average av_ra;
//        trapezoidal_average av_rso;
//        trapezoidal_average av_rs;
//        std::uniform_real_distribution<double> ur(0.0, 250.0);
//        std::default_random_engine gen;
//        std::cout << "========= Slope 90S =======" << std::endl;
//        SUBCASE("June") {
//            std::cout << "========= June ========" << std::endl;
//            ta = utc_cal.time(1970, 06, 21, 00, 00, 0, 0);
//            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
//            av_ra.initialize(rad.ra_radiation(), 0.0);
//            av_rso.initialize(r.psw_radiation, 0.0);
//            av_rs.initialize(r.tsw_radiation, 0.0);
//            for (int h = 1; h < 24; ++h) {
//                t = utc_cal.time(1970, 06, 21, h, 00, 0, 0); // June
//                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//                av_rahor.add(rad.ra_radiation_hor(), h);
//                av_ra.add(rad.ra_radiation(), h);
//                av_rso.add(r.psw_radiation, h);
//                av_rs.add(r.tsw_radiation, h);
//            }
//            std::cout << "rahor: " << av_rahor.result() << std::endl;
//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rso: " << av_rso.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
//                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(110.0).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(85.0).epsilon(0.05));
//        }
//        SUBCASE("January") {
//            std::cout << "========= January ========" << std::endl;
//            ta = utc_cal.time(1970, 01, 1, 00, 00, 0, 0);
//            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
//            av_ra.initialize(rad.ra_radiation(), 0.0);
//            av_rso.initialize(r.psw_radiation, 0.0);
//            av_rs.initialize(r.tsw_radiation, 0.0);
//            for (int h = 1; h < 24; ++h) {
//                t = utc_cal.time(1970, 01, 1, h, 00, 0, 0); // June
//                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//                av_rahor.add(rad.ra_radiation_hor(), h);
//                av_ra.add(rad.ra_radiation(), h);
//                av_rso.add(r.psw_radiation, h);
//                av_rs.add(r.tsw_radiation, h);
//            }
//            std::cout << "rahor: " << av_rahor.result() << std::endl;
//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rso: " << av_rso.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
//                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(410.0).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(200.0).epsilon(0.05));
//        }
//        SUBCASE("December") {
//            std::cout << "========= December ========" << std::endl;
//            ta = utc_cal.time(1970, 12, 12, 00, 00, 0, 0);
//            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
//            av_ra.initialize(rad.ra_radiation(), 0.0);
//            av_rso.initialize(r.psw_radiation, 0.0);
//            av_rs.initialize(r.tsw_radiation, 0.0);
//            for (int h = 1; h < 24; ++h) {
//                t = utc_cal.time(1970, 12, 12, h, 00, 0, 0); // June
//                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//                av_rahor.add(rad.ra_radiation_hor(), h);
//                av_ra.add(rad.ra_radiation(), h);
//                av_rso.add(r.psw_radiation, h);
//                av_rs.add(r.tsw_radiation, h);
//            }
//            std::cout << "rahor: " << av_rahor.result() << std::endl;
//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rso: " << av_rso.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
//                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(410.0).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(200.0).epsilon(0.05));
//        }
//
//    }
//    TEST_CASE("check_solar_radiation_slope_90N"){
//        parameter p;
//        response r;
//        p.albedo = 0.05;
//        p.turbidity = 1.0;
//        calculator<parameter,response> rad(p);
//        calendar utc_cal;
//        double lat = 44.0;
//        utctime t;
//        // checking for horizontal surface Eugene, OR, p.64, fig.1f
//        double slope = 90*shyft::core::radiation::pi/180; // 90 N
//        double proj = sin(slope);
//        double aspect = 180*shyft::core::radiation::pi/180;// facing north
//        arma::vec surface_normal({proj*cos(aspect),proj*sin(aspect),cos(slope)});
//        double ra_sum = 0.0;
//        double rso_sum = 0.0;
//        double rahor_sum = 0.0;
//        utctime ta;
//        trapezoidal_average av_rahor;
//        trapezoidal_average av_ra;
//        trapezoidal_average av_rso;
//        trapezoidal_average av_rs;
//        std::uniform_real_distribution<double> ur(10.0, 100.0);
//        std::default_random_engine gen;
//        std::cout << "========= Slope 90N =======" << std::endl;
//        SUBCASE("June") {
//            std::cout << "========= June ========" << std::endl;
//            ta = utc_cal.time(1970, 06, 21, 00, 00, 0, 0);
//            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
//            av_ra.initialize(rad.ra_radiation(), 0.0);
//            av_rso.initialize(r.psw_radiation, 0.0);
//            av_rs.initialize(r.tsw_radiation, 0.0);
//            for (int h = 1; h < 24; ++h) {
//                t = utc_cal.time(1970, 06, 21, h, 00, 0, 0); // June
//                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//                av_rahor.add(rad.ra_radiation_hor(), h);
//                av_ra.add(rad.ra_radiation(), h);
//                av_rso.add(r.psw_radiation, h);
//                av_rs.add(r.tsw_radiation, h);
//            }
//            std::cout << "rahor: " << av_rahor.result() << std::endl;
//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rso: " << av_rso.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
//                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(105.0).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(50.0).epsilon(0.05));
//        }
//        SUBCASE("January") {
//            std::cout << "========= January ========" << std::endl;
//            ta = utc_cal.time(1970, 01, 1, 00, 00, 0, 0);
//            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
//            av_ra.initialize(rad.ra_radiation(), 0.0);
//            av_rso.initialize(r.psw_radiation, 0.0);
//            av_rs.initialize(r.tsw_radiation, 0.0);
//            for (int h = 1; h < 24; ++h) {
//                t = utc_cal.time(1970, 01, 1, h, 00, 0, 0); // June
//                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//                av_rahor.add(rad.ra_radiation_hor(), h);
//                av_ra.add(rad.ra_radiation(), h);
//                av_rso.add(r.psw_radiation, h);
//                av_rs.add(r.tsw_radiation, h);
//            }
//            std::cout << "rahor: " << av_rahor.result() << std::endl;
//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rso: " << av_rso.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
//                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(0.0).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(9.0).epsilon(0.05));
//        }
//        SUBCASE("December") {
//            std::cout << "========= December ========" << std::endl;
//            ta = utc_cal.time(1970, 12, 12, 00, 00, 0, 0);
//            rad.tsw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
//            av_ra.initialize(rad.ra_radiation(), 0.0);
//            av_rso.initialize(r.psw_radiation, 0.0);
//            av_rs.initialize(r.tsw_radiation, 0.0);
//            for (int h = 1; h < 24; ++h) {
//                t = utc_cal.time(1970, 12, 12, h, 00, 0, 0); // June
//                rad.tsw_radiation(r, lat, t, surface_normal, 20.0, 50.0, 150.0, ur(gen));
//                av_rahor.add(rad.ra_radiation_hor(), h);
//                av_ra.add(rad.ra_radiation(), h);
//                av_rso.add(r.psw_radiation, h);
//                av_rs.add(r.tsw_radiation, h);
//            }
//            std::cout << "rahor: " << av_rahor.result() << std::endl;
//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rso: " << av_rso.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
//                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(0.0).epsilon(0.05));
//                    FAST_CHECK_EQ(av_rso.result(), doctest::Approx(9.0).epsilon(0.05));
//        }
//
//    }
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
