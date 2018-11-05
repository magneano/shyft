#include "test_pch.h"
#include "core/utctime_utilities.h"
#include "core/geo_cell_data.h"
#include <armadillo>
#include <vector>
#include <chrono>
#include <boost/math/constants/constants.hpp>
#include "core/radiation_dingman.h"

namespace shyft::test {
    using shyft::core::geo_cell_data;
    struct cell {
        geo_cell_data geo;
    };
}

TEST_SUITE("radiation_dingman") {
    using shyft::core::radiation_dingman::parameter;
    using shyft::core::radiation_dingman::calculator;
    using shyft::core::radiation_dingman::surface_normal;
    using shyft::core::calendar;
    using shyft::test::cell;
    using std::vector;
    using shyft::core::utctime;
    // test basics: creation, etc

    TEST_CASE("check_radiation"){
        /** \brief Check Earth Hour Angle
         * based on ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, no EOT correction provided*/
        parameter p;
        calculator ra(p);
        calendar utc_cal;
        double lat = 45.0;
        double rh = 40;
        double temp = 273.16;
        double slope = 40;
        double aspect = 0.0;
        utctime t1, t2, t3;
        //arma::vec surface_normal({0.0,0.0,1.0});
        t1 = utc_cal.time(1970, 03, 31, 10, 30, 0, 0);
        ra.compute_ra_radiation(lat,rh,temp,t1);
                FAST_CHECK_EQ(ra.rcs_radiation_, doctest::Approx(0.0));// clear sky radiation
                FAST_CHECK_EQ(ra.slope_factor_, doctest::Approx(1.0));// slope factor
        ra.compute_ra_radiation(lat,rh,temp,t1,slope,aspect );
                FAST_CHECK_EQ(ra.rcs_radiation_, doctest::Approx(0.0));// clear sky radiation
                FAST_CHECK_EQ(ra.slope_factor_, doctest::Approx(1.0));// slope factor

    }

}
