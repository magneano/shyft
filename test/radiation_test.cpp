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
const double pi = boost::math::constants::pi<double>();

struct radiation_calculator {
    calendar utc;
    double doy;
    double hour;
    double s = 0.0; // horizontal slope
    double delta = 0.1; // declination of the earth, positive for northern hemisphere summer
    double gamma = 0.0; // surface aspect angle, 0 -- due south, -pi/2 -- due east, pi/2 -- due west, +-pi -- due north
    double phi = 0.0 ;// latitude
    double ra = 0; // extraterrestrial solar radiation [W/m2]
    double rso = 0; // clear sky solar radiation [W/m2]
    
    /** \tparam V is typically arma::vec */
    template <class V>
    double real_radiation(double latitude, utctime t, const V & surface_normal) {
        doy=utc.day_of_year(t);
        hour = utc.calendar_units(t).hour;
        return hour;
    }
    const double gsc = 1367; // W/m2 -- solar constant
    /** \brief computes instantaneous extraterrestrial solar radiation*/
    template <class V>
    double ra_radiation(double latitude, utctime t, const V & surface_normal){
        phi = latitude;
        doy = utc.day_of_year(t);
        delta = 23.45*pi/180*sin(2*pi*(284+doy)/36.25); // declination angle based on eq.(3.1) https://www.itacanet.org/the-sun-as-a-source-of-energy/part-3-calculating-solar-angles/
        hour = utc.calendar_units(t).hour;
        compute_costt(hour);
        ra = gsc*cos_tt*(1+0.0033*cos(doy*2*pi/365)); // eq.(1)
        return ra;
    }
    /*double rso_cs_radiation(){
        return ra*tau_swo;
    }*/
private:
    double a, b, c, cos_tt;
    /** \brief computes necessary geometric parameters
     * \param omega -- hour angle
     * */
    void compute_costt(double omega){
        a = sin(delta)*cos(phi)*sin(s)*cos(gamma) - sin(delta)*sin(phi)*cos(s); // eq.11a
        b = cos(delta)*cos(phi)*cos(s) + cos(delta)*sin(phi)*sin(s)*cos(gamma);//eq 11b
        c = cos(delta)*sin(s)*sin(gamma);
        cos_tt = -a + b*cos(omega) + c*sin(omega); // eq.14
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

namespace shyft::test {
    using shyft::core::geo_cell_data;
    struct cell {
        geo_cell_data geo;
    };
}

TEST_SUITE("radiation") {
    using shyft::core::radiation_calculator;
    using shyft::core::surface_normal;
    using shyft::core::calendar;
    using shyft::test::cell;
    using std::vector;
    
    TEST_CASE("compute_real_radiation") {
        radiation_calculator r;
        FAST_CHECK_EQ(r.real_radiation(1.0,calendar::HOUR*25,1),doctest::Approx(1));
    }
    
    TEST_CASE("surface_normal_from_cells") {
        vector<cell> cells;
        auto r= surface_normal(cells);
    }
}
