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
    using shyft::core::radiation_dingman::response;
    using shyft::core::radiation_dingman::calculator;
    using shyft::core::radiation_dingman::surface_normal;
    using shyft::core::calendar;
    using shyft::test::cell;
    using std::vector;
    using shyft::core::utctime;
    // test basics: creation, etc

    TEST_CASE("points_slope_factor"){
        /** \brief Check Earth Hour Angle
         * based on ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, no EOT correction provided*/
        parameter p;
        response r;
        calculator ra(p,r);
        calendar utc_cal;
        double lat = 45.0;
        double rh = 40;
        double temp = 273.16;
        double slope = 40;
        double aspect = 0.0;
        utctime t1, t2, t3;
        //arma::vec surface_normal({0.0,0.0,1.0});
        t1 = utc_cal.time(1970, 03, 31, 10, 30, 0, 0);
        ra.compute_ra_radiation(r,lat,rh,temp,t1);
                FAST_CHECK_EQ(r.slope_factor, doctest::Approx(1.0));// slope factor
        ra.compute_ra_radiation(r, lat,rh,temp,t1,slope,aspect );
                FAST_CHECK_EQ(r.slope_factor, doctest::Approx(1.35));// slope factor

    }
    TEST_CASE("regression"){
        parameter p;
        response r;
        calculator ra(p,r);
        calendar utc_cal;
        double lat = 45.0;
        double rh = 40;
        double temp = 273.16;
        double slope = 0.0;
        double aspect = 0.0;
        utctime t1, t2, t3;
        //arma::vec surface_normal({0.0,0.0,1.0});
        t1 = utc_cal.time(1970, 03, 31, 10, 30, 0, 0);
        double r_prev = 0.0;
        double r_cur = 0.0;
        for(double lat=60.0; lat > 20.0; lat -= 1.0) {
            ra.compute_ra_radiation(r,lat,rh,temp,t1);
            r_cur = r.db_radiation;
            std::cout<<"Kdb: " << r_cur<<std::endl;
            TS_ASSERT(r_prev < r_cur);//the radiation should increase closer to equator
            r_prev = r_cur;
        }
    }

}
