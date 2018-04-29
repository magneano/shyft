#include "test_pch.h"

#include <cmath>
#include <vector>

#include "core/utctime_utilities.h"
#include "core/time_axis.h"
#include "core/time_series_dd.h"



TEST_SUITE("time_series") {

    using shyft::core::no_utctime;
    using std::numeric_limits;
    const double eps = numeric_limits<double>::epsilon();
    using shyft::time_series::dd::apoint_ts;
    using shyft::time_series::dd::bit_decoder;
    using shyft::time_series::dd::bitmask;
    using shyft::time_axis::generic_dt;
    using shyft::time_series::ts_point_fx;
    using std::vector;
    using std::make_shared;
    using std::isfinite;
    
       TEST_CASE("ts_decode") {
        

        //
        // 0. just to ensure that the bitmask function do the right thing:
        //
        uint64_t expected[]={
            0x00000000,0x00000001,0x00000003,0x00000007,0x0000000f,0x0000001f,0x0000003f,0x0000007f,0x000000ff,0x000001ff,0x000003ff,0x000007ff,
            0x00000fff,0x00001fff,0x00003fff,0x00007fff,0x0000ffff,0x0001ffff,0x0003ffff,0x0007ffff,0x000fffff,0x001fffff,0x003fffff,0x007fffff,
            0x00ffffff,0x01ffffff,0x03ffffff,0x07ffffff,0x0fffffff,0x1fffffff,0x3fffffff,0x7fffffff,0xffffffff,0x1ffffffff,

        };
        for(unsigned n_bits=0;n_bits<33;++n_bits) {
            FAST_CHECK_EQ(bitmask<uint64_t>(n_bits),expected[n_bits]);
        }
        //
        // 1. checkout fast the 10 lowest bits, 
        //
        bit_decoder b(0,1);
        for(uint64_t i=0;i< (1UL<<10);++i) {
            double x=i;
            FAST_CHECK_EQ(b.decode(x), doctest::Approx( (i&1?1.0:0.0)));
        }
        //
        // 2. some really large numbers
        //
        double x0=(uint64_t(1)<<51);
        double x1=x0+1;
        FAST_CHECK_EQ(b.decode(x0), doctest::Approx(0.0));
        FAST_CHECK_EQ(b.decode(x1), doctest::Approx(1.0));
        //
        // 3. test it can extract 0..n-bits from some place within the integer
        //
        for(unsigned int start_bit=0;start_bit<10;++start_bit){
            for(unsigned int n_bits=1;n_bits<10;++n_bits) {
                uint64_t max_val= (1UL<<n_bits)-1;
                uint64_t one_val=1;
                uint64_t min_val=0;
                bit_decoder bd{start_bit,n_bits};
                FAST_CHECK_EQ(bd.decode( double(max_val<<start_bit)),doctest::Approx(double(max_val)));
                FAST_CHECK_EQ(bd.decode( double(min_val<<start_bit)),doctest::Approx(double(min_val)));
                FAST_CHECK_EQ(bd.decode( double(one_val<<start_bit)),doctest::Approx(double(one_val)));
            }
        }
        //
        // 4. test it can extract isolated bit-groups from same integer
        //
        
        bit_decoder b_0_1(0,1);
        bit_decoder b_1_3(1,3);
        bit_decoder b_9_1(9,1);
        double x_1_2_1= (1<<9) + (2<<1) + (1);
        double x_0_7_0= (0<<9) + (7<<1) + (0);
        double x_1_5_1= (1<<9) + (5<<1) + (1);
        FAST_CHECK_EQ(b_0_1.decode(x_1_2_1), doctest::Approx(1));
        FAST_CHECK_EQ(b_1_3.decode(x_1_2_1), doctest::Approx(2));
        FAST_CHECK_EQ(b_9_1.decode(x_1_2_1), doctest::Approx(1));

        FAST_CHECK_EQ(b_0_1.decode(x_0_7_0), doctest::Approx(0));
        FAST_CHECK_EQ(b_1_3.decode(x_0_7_0), doctest::Approx(7));
        FAST_CHECK_EQ(b_9_1.decode(x_0_7_0), doctest::Approx(0));

        FAST_CHECK_EQ(b_0_1.decode(x_1_5_1), doctest::Approx(1));
        FAST_CHECK_EQ(b_1_3.decode(x_1_5_1), doctest::Approx(5));
        FAST_CHECK_EQ(b_9_1.decode(x_1_5_1), doctest::Approx(1));
        
        //
        // 5. test it gives nan for conditions not well defined
        //
        FAST_CHECK_UNARY_FALSE(isfinite(b_0_1.decode(shyft::nan))); // nan -> nan
        FAST_CHECK_UNARY_FALSE(isfinite(b_0_1.decode(-2.3))); // neg. value -> nan
        FAST_CHECK_UNARY_FALSE(isfinite(b_0_1.decode(double(uint64_t(1)<<52) +1.0))); // not precise representative value -> nan
    }

}
