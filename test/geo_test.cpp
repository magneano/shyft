#include "test_pch.h"
#include <vector>
#include <algorithm>
#include "core/geo_cell_data.h"

using shyft::core::geo_cell_data;
using shyft::core::land_type_fractions;

TEST_SUITE("geo_cell_data") {
TEST_CASE("land_type_fractions") {
    land_type_fractions a;
    FAST_CHECK_EQ(a.glacier(),doctest::Approx(0.0));
    FAST_CHECK_EQ(a.lake(),doctest::Approx(0.0));
    FAST_CHECK_EQ(a.reservoir(),doctest::Approx(0.0));
    FAST_CHECK_EQ(a.forest(),doctest::Approx(0.0));
    FAST_CHECK_EQ(a.unspecified(),doctest::Approx(1.0));
    FAST_CHECK_EQ(a.snow_storage(),doctest::Approx(1.0));
    land_type_fractions b;
    FAST_CHECK_EQ(a,b);
    land_type_fractions c(0.1,0.11,0.123,0.1234,0.5436);
    FAST_CHECK_EQ(c.glacier(),doctest::Approx(0.1));
    FAST_CHECK_EQ(c.lake(),doctest::Approx(0.11));
    FAST_CHECK_EQ(c.reservoir(),doctest::Approx(0.123));
    FAST_CHECK_EQ(c.forest(),doctest::Approx(0.1234));
    FAST_CHECK_EQ(c.unspecified(),doctest::Approx(0.5436));
    FAST_CHECK_EQ(c.snow_storage(),doctest::Approx(1-(0.11+0.123)));
    CHECK(a!=c);
    
    
}
}
