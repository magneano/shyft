#include "test_pch.h"

#include <cmath>
#include <vector>

#include "core/utctime_utilities.h"
#include "core/time_axis.h"
#include "core/time_series_dd.h"

namespace shyft::time_series::dd {
extern vector<double> bucket_fix(vector<double> const& v, size_t i0, size_t n);
}

TEST_SUITE("time_series") {

    using shyft::core::no_utctime;
	using shyft::core::utctime;
    using std::numeric_limits;
	using std::chrono::seconds;
    //const double eps = numeric_limits<double>::epsilon();
    using shyft::time_series::dd::apoint_ts;
    using shyft::time_series::dd::bit_decoder;
    using shyft::time_series::dd::bitmask;
	using shyft::time_series::dd::bucket_ts;
	using shyft::time_series::dd::bucket_parameter;
	using shyft::core::calendar;
    using shyft::time_axis::generic_dt;
    using shyft::time_series::ts_point_fx;
    using std::vector;
    using std::make_shared;
    using std::isfinite;
    vector<double> data{ 0, 0, 0, 0, 0, 0,
                         4, 4, 4, 4, 4, 4,
                         1, 1, 1, 1, 1, 1,
                         9, 9, 9, 9, 9, 9,
                         9, 9, 9, 9, 9, 9,
                         -6,-6,-6,-6,-6,-6,
                         6, 6, 6, 6, 6, 6,
                         12,12,12,12,12,12,
                         12,12,12,12,12,12,
                         12,12,12,12,12,12,
                         12,12,12,12,12,12,
                         12,12,12,12,12,12};

    TEST_CASE("bucket_ts_normal") {
        generic_dt ta(seconds(0), calendar::HOUR, data.size());
        apoint_ts ts_raw(ta, data, ts_point_fx::POINT_AVERAGE_VALUE);
        bucket_ts bts(ts_raw, bucket_parameter{ seconds(0),-100.0 });
        auto vv = bts.values();
        vector<double> expected{0, 0, 0, 0, 0, 0,
                                3, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                6, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                2, 0, 0, 0, 0, 0,
                                1, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0 };
        //std::cout << "\nexpected:\n";
        //for (const auto v : expected) std::cout << v << " ";
        //std::cout << "\ngot:\n";
        //for (const auto v : vv) std::cout << v << " ";
        //std::cout << "\nvalue:\n";
        //for (int i = 0; i < vv.size(); ++i) std::cout << i << " " << bts.value(i) << "\n";
        for (int i = 0; i < vv.size(); ++i) {
            FAST_CHECK_EQ(vv[i], doctest::Approx(expected[i]));
            FAST_CHECK_EQ(bts.value(i),doctest::Approx(expected[i]));
            FAST_CHECK_EQ(bts.value_at(ta.time(i)),doctest::Approx(expected[i]));
        }
    }
	TEST_CASE("bucket_ts_offset") {
            using shyft::nan;
			generic_dt ta(seconds(0), calendar::HOUR, data.size());
			apoint_ts ts_raw(ta, data, ts_point_fx::POINT_AVERAGE_VALUE);
			bucket_ts bts(ts_raw,bucket_parameter{seconds(4*3600),-1000.0});
			vector<double> ee{nan,nan,nan,nan,0, 0,
                                3, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                6, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                2, 0, 0, 0, 0, 0,
                                1, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
            };
			auto vv = bts.values();
            //TODO: consider if we should clip front/tail time-axis
            //std::cout<<"\nexpected:\n";
            //for (const auto v : ee) std::cout << v << " ";
            //std::cout<<"\ngot:\n";
			//for (const auto v : vv) std::cout << v << " ";
            
			for (int i = 4; i < vv.size(); ++i) {
				FAST_CHECK_EQ(vv[i], doctest::Approx(ee[i]));
                FAST_CHECK_EQ(bts.value(i),doctest::Approx(ee[i]));
                FAST_CHECK_EQ(bts.value_at(ta.time(i)),doctest::Approx(ee[i]));
            }
   }
}
