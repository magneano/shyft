#include "test_pch.h"
#include "core/utctime_utilities.h"
#include "core/geo_cell_data.h"
#include <armadillo>
#include <vector>

namespace shyft::core {
using std::vector;


struct radiation_calculator {
    calendar utc;
    
    /** \tparam V is typically arma::vec */
    template <class V>
    double real_radiation(double latitude, utctime t, const V & surface_normal) {
        auto doy=utc.day_of_year(t);
        auto hour = utc.calendar_units(t).hour;
        return hour;
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
