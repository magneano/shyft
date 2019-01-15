#include "test_pch.h"

#include <cmath>
#include <vector>

#include "core/utctime_utilities.h"
#include "core/time_axis.h"
#include "core/time_series_dd.h"

static inline shyft::core::utctime _t(int64_t t1970s) {return shyft::core::utctime{shyft::core::seconds(t1970s)};}

TEST_SUITE("time_series") {

    using shyft::core::no_utctime;
    using shyft::core::seconds;
    using std::numeric_limits;
    const double eps = numeric_limits<double>::epsilon();
    using shyft::time_series::dd::apoint_ts;
    using shyft::time_axis::generic_dt;
    using shyft::time_series::ts_point_fx;
    using std::vector;
    using std::make_shared;
    using std::isfinite;
    using shyft::time_series::dd::qac_parameter;
    using shyft::time_series::dd::qac_ts;
    using namespace shyft::time_series::dd;
    TEST_CASE("qac_is_ok") {
        FAST_CHECK_UNARY(qac::is_ok_quality_lower_bound(0.0,1.0));
        FAST_CHECK_UNARY(!qac::is_ok_quality_lower_bound(1.0,0.0));
        FAST_CHECK_UNARY(qac::is_ok_quality_upper_bound(0.0,1.0));
        FAST_CHECK_UNARY(!qac::is_ok_quality_upper_bound(1.0,0.0));
        FAST_CHECK_UNARY(qac::is_ok_quality_in_range(0.0,0.5,1.0));
        FAST_CHECK_UNARY(!qac::is_ok_quality_in_range(0.0,1.5,1.0));
    }
    TEST_CASE("qac_is_repeated_once") {
        qac_parameter p;
        auto a=make_shared<apoint_ts>(generic_dt{seconds(0),seconds(10),5},vector<double>{1.0,1.0,3.0,4.0,5.0},ts_point_fx::POINT_AVERAGE_VALUE);
        FAST_CHECK_UNARY(qac::is_repeated_once(a,1,a->value(1),p));
        FAST_CHECK_UNARY(!qac::is_repeated_once(a,1,shyft::nan,p));// inf/nan should ->false
        
        FAST_CHECK_UNARY(!qac::is_repeated_once(a,0,a->value(0),p));
        FAST_CHECK_UNARY(!qac::is_repeated_once(a,2,a->value(2),p));
        p.repeat_allowed.push_back(1.0001);// slighly off, within repeat_tolerance
        FAST_CHECK_UNARY(!qac::is_repeated_once(a,1,a->value(1),p));
        p.repeat_tolerance=1e-9;// tighten repeat_tolerance
        FAST_CHECK_UNARY(qac::is_repeated_once(a,1,a->value(1),p));
    }
    TEST_CASE("qac_find_last_valid_repeat") {
        qac_parameter p;//                           last valid repeat is here                 v 
        auto a=make_shared<apoint_ts>(generic_dt{seconds(0),seconds(10),5},vector<double>{1.01,1.02,1.01,1.01,5.0},ts_point_fx::POINT_AVERAGE_VALUE);
        p.repeat_tolerance=0.1;
        p.repeat_timespan=seconds(9);FAST_CHECK_EQ(0,qac::find_last_valid_repeat(a,3,1.01,p));
        p.repeat_timespan=seconds(12);FAST_CHECK_EQ(1,qac::find_last_valid_repeat(a,3,1.01,p));
        p.repeat_timespan=seconds(20);FAST_CHECK_EQ(2,qac::find_last_valid_repeat(a,3,1.01,p));
        p.repeat_tolerance=0.0015;
        p.repeat_timespan=seconds(9);FAST_CHECK_EQ(2,qac::find_last_valid_repeat(a,3,1.01,p));
        p.repeat_timespan=seconds(12);FAST_CHECK_EQ(3,qac::find_last_valid_repeat(a,3,1.01,p));
        p.repeat_timespan=seconds(20);FAST_CHECK_EQ(3,qac::find_last_valid_repeat(a,3,1.01,p));
    }
    TEST_CASE("qac_find_left_ok_value") {
        qac_parameter p;
        auto a=make_shared<apoint_ts>(generic_dt{seconds(0),seconds(10),6},vector<double>{1.01,1.01,1.01,-999.0,shyft::nan,5.0},ts_point_fx::POINT_AVERAGE_VALUE);
        p.min_x=0.0;
        p.repeat_timespan=seconds(9);FAST_CHECK_EQ(0,qac::find_left_ok_value(a,4,p));
        p.repeat_timespan=seconds(10);FAST_CHECK_EQ(1,qac::find_left_ok_value(a,4,p));
        p.repeat_timespan=seconds(20);FAST_CHECK_EQ(2,qac::find_left_ok_value(a,4,p));
        p.repeat_timespan=seconds(30);FAST_CHECK_EQ(2,qac::find_left_ok_value(a,4,p));
        p.repeat_allowed.push_back(1.01);
        p.repeat_timespan=seconds(10);FAST_CHECK_EQ(2,qac::find_left_ok_value(a,4,p));
        p.min_x=-1000.0;
        p.repeat_timespan=seconds(30);FAST_CHECK_EQ(3,qac::find_left_ok_value(a,4,p));        
    }
    TEST_CASE("qac_find_first_ok_value_right") {
        qac_parameter p;
        auto a=make_shared<apoint_ts>(generic_dt{seconds(0),seconds(10),6},vector<double>{1.01,1.01,1.01,1.01,shyft::nan,5.0},ts_point_fx::POINT_AVERAGE_VALUE);
        auto b=make_shared<apoint_ts>(generic_dt{seconds(0),seconds(10),6},vector<double>{1.01,1.01,1.01,1.01,1.01,1.01},ts_point_fx::POINT_AVERAGE_VALUE);
        p.min_x=0.0;
        p.max_timespan=seconds(10*10);// allow to find any value
        p.repeat_timespan=seconds(10);FAST_CHECK_EQ(5,qac::find_first_ok_value_right(a,2,a->value(2),a->time(2),p));
        p.max_timespan=seconds(10*2);FAST_CHECK_EQ(2,qac::find_first_ok_value_right(a,2,a->value(2),a->time(2),p));
        p.max_timespan=seconds(10*2);FAST_CHECK_EQ(2,qac::find_first_ok_value_right(b,2,b->value(2),b->time(2),p));
    }
    TEST_CASE("qac_find_first_ok_value_right_no_repeat") {
        qac_parameter p;
        auto a=make_shared<apoint_ts>(generic_dt{seconds(0),seconds(10),6},vector<double>{1.01,1.01,1.01,1.01,shyft::nan,5.0},ts_point_fx::POINT_AVERAGE_VALUE);
        p.repeat_timespan=seconds(10);
        p.max_timespan=seconds(10*10);FAST_CHECK_EQ(3,qac::find_first_ok_value_right_no_repeat(a,2,a->time(2),p));
        p.max_timespan=seconds(10*2);FAST_CHECK_EQ(5,qac::find_first_ok_value_right_no_repeat(a,3,a->time(3),p));
        p.max_timespan=seconds(10);FAST_CHECK_EQ(3,qac::find_first_ok_value_right_no_repeat(a,3,a->time(3),p));
    }
    TEST_CASE("qac_legal_repeat") {
        qac_parameter p;
        FAST_CHECK_EQ(false,qac::legal_repeat(p,3.5));
        p.repeat_allowed.push_back(3.5);
        FAST_CHECK_EQ(true,qac::legal_repeat(p,3.5));
        p.repeat_allowed.push_back(5.0);
        FAST_CHECK_EQ(true,qac::legal_repeat(p,5.0));
        p.repeat_tolerance=0.11;
        FAST_CHECK_EQ(true,qac::legal_repeat(p,3.6));
        FAST_CHECK_EQ(true,qac::legal_repeat(p,3.4));
    }
    TEST_CASE("qac_repeated_value") {
        qac_parameter p;
        p.repeat_tolerance=0.1;
        FAST_CHECK_EQ(true,qac::repeated_value(p,1.1,1.1));
        FAST_CHECK_EQ(false,qac::repeated_value(p,1.0,1.100001));
        FAST_CHECK_EQ(false,qac::repeated_value(p,1.0,0.899999));
        FAST_CHECK_EQ(true,qac::repeated_value(p,1.0,0.91));
    }
    TEST_CASE("qac_fill_fx_for_bad_values") {
        qac_parameter p;        
        generic_dt ta{seconds(0),seconds(10),7};
        vector<double> v{1.01,1.01,1.01,9.0,shyft::nan,2.0,2.0};
        SUBCASE("fill_min_max_only") {
            p.min_x=0.0;
            p.max_x=3.0;
            size_t n_fills=0;
            qac::fill_fx_for_bad_values(ta,v,p,[&n_fills](size_t,size_t i)->double {++n_fills;return double(i);});
            FAST_CHECK_EQ(n_fills,2);
            FAST_CHECK_EQ(v[3],doctest::Approx(double(3)));
            FAST_CHECK_EQ(v[4],doctest::Approx(double(4)));
        }
        SUBCASE("fill_repeated_and_min_max") {
            p.min_x=0.0;
            p.max_x=3.0;
            size_t n_fills=0;
            p.repeat_timespan=seconds(10);
            qac::fill_fx_for_bad_values(ta,v,p,[&n_fills](size_t l,size_t i)->double {++n_fills;return double(i)+(l==std::string::npos?0.0:double(l)/10.0);});
            FAST_CHECK_EQ(n_fills,3);
            FAST_CHECK_EQ(v[2],doctest::Approx(double(2.0)));
            FAST_CHECK_EQ(v[3],doctest::Approx(double(3.0)));// because repeated sequence,not ok to extend
            FAST_CHECK_EQ(v[4],doctest::Approx(double(4.1)));
        }
    }
    TEST_CASE("qac_fill_constant_for_missing_values") {
        vector<double> v{1.0,shyft::nan,3.0};
        qac::fill_constant_for_missing_values(v,2.0);
        for(size_t i=0;i<v.size();++i)
            FAST_CHECK_EQ(v[i],doctest::Approx(double(i+1.0)));
    }
    TEST_CASE("qac_fill_ts_for_missing_values") {
        vector<double> v{1.0,shyft::nan,3.0};
        auto cts=make_shared<apoint_ts>(generic_dt{seconds(0),seconds(10),3},vector<double>{1.0,2.0,3.0},ts_point_fx::POINT_AVERAGE_VALUE);
        qac::fill_ts_for_missing_values(cts->time_axis(),v,cts->sts());
        for(size_t i=0;i<v.size();++i)
            FAST_CHECK_EQ(v[i],doctest::Approx(double(i+1.0)));
    }
    TEST_CASE("qac_fill_linear_for_missing_values") {
        vector<double> v{1.0,shyft::nan,3.0,shyft::nan,shyft::nan,6.0};
        generic_dt ta{seconds(0),seconds(10),6};
        SUBCASE("can_fill_all") {
            qac_parameter p;p.max_timespan=seconds(30);
            qac::fill_linear_for_missing_values(ta,v,p);
            for(size_t i=0;i<v.size();++i)
                FAST_CHECK_EQ(v[i],doctest::Approx(double(i+1.0)));
        }
        SUBCASE("can_fill_some") {
            qac_parameter p;p.max_timespan=seconds(20);
            qac::fill_linear_for_missing_values(ta,v,p);
            for(size_t i=0;i<3;++i)
                FAST_CHECK_EQ(v[i],doctest::Approx(double(i+1.0)));
            FAST_CHECK_EQ(v[5],doctest::Approx(6.0));
            FAST_CHECK_UNARY(!isfinite(v[4]));
        }
        SUBCASE("can_fill_none") {
            qac_parameter p;p.max_timespan=seconds(19);
            qac::fill_linear_for_missing_values(ta,v,p);
            for(auto i:vector<size_t>{0,2,5})
                FAST_CHECK_EQ(v[i],doctest::Approx(double(i+1.0)));
            for(auto i:vector<size_t>{1,3,4})
                FAST_CHECK_UNARY(!isfinite(v[i]));
        }
        SUBCASE("bad_start") {
            v[0]=shyft::nan;
            qac_parameter p;p.max_timespan=seconds(30);
            qac::fill_linear_for_missing_values(ta,v,p);
            for(auto i:vector<size_t>{2,3,4,5})
                FAST_CHECK_EQ(v[i],doctest::Approx(double(i+1.0)));
            for(auto i:vector<size_t>{0,1})
                FAST_CHECK_UNARY(!isfinite(v[i]));
            
        }
    }
    TEST_CASE("qac_parameter") {
        using shyft::time_series::dd::qac::is_ok_quality;
        using namespace shyft::time_series::dd;
        qac_parameter q;

        SUBCASE("no limits set, allow all values, except nan") {
            FAST_CHECK_EQ(is_ok_quality(q,shyft::nan),false);
            FAST_CHECK_EQ(is_ok_quality(q,1.0),true);
        }

        SUBCASE("min/max abs limits") {
            q.max_x=1.0;
            FAST_CHECK_EQ(is_ok_quality(q,1.0),true);
            FAST_CHECK_EQ(is_ok_quality(q,1.0+eps),false);
            q.min_x=-1.0;
            FAST_CHECK_EQ(is_ok_quality(q,-1.0),true);
            FAST_CHECK_EQ(is_ok_quality(q,-1.0-eps),false);
            
        }
    }

    TEST_CASE("qac_ts") {
        generic_dt ta{_t(0),seconds(10),5};
        //                 0    1       2     3    4
        vector<double> v {0.0,1.0,shyft::nan,3.0,-20.1};
        vector<double>ev {0.0,1.0,       2.0,3.0,-20.1};
        vector<double>cv {1.0,0.0,      -1.0,3.0,-20.1};
        apoint_ts src(ta,v,ts_point_fx::POINT_AVERAGE_VALUE);
        apoint_ts cts{ta,cv,ts_point_fx::POINT_AVERAGE_VALUE};
        qac_parameter qp;
        auto ts = make_shared<qac_ts>(src,qp);


        // verify simple min-max limit cases
        FAST_CHECK_UNARY(ts.get()!=nullptr);
        FAST_CHECK_EQ(ts->value(2),doctest::Approx(1.0));
        FAST_CHECK_EQ(ts->value_at(ts->time(2)),doctest::Approx(1.0));
        FAST_CHECK_EQ(ts->value_at(ts->time(2)+seconds(1)),doctest::Approx(1.0));
        src.set_point_interpretation(ts_point_fx::POINT_INSTANT_VALUE);
        FAST_CHECK_EQ(ts->value_at(ts->time(2)+ seconds(1) ),doctest::Approx(2.1));
        FAST_CHECK_EQ(ts->value_at(ts->time(2)- seconds(1) ),doctest::Approx(1.9));
        ts->p.min_x = 0.0;
        FAST_CHECK_UNARY(!isfinite(ts->value_at(ts->time(3)+ seconds(1) )));
        ts->p.min_x = -40.0;
        FAST_CHECK_UNARY(isfinite(ts->value_at(ts->time(3)+ seconds(1) )));
        ts->p.max_x = 2.9; // clip 3.0 out of range
        FAST_CHECK_EQ(ts->value(2),
                      doctest::Approx(
                         v[1] + 10*(v[4] - v[1])/(40 - 10)
                      )
        );

        FAST_CHECK_EQ(ts->value(3),doctest::Approx(
                         v[1] + 20*(v[4] - v[1])/(40 - 10)
                      )
        );

        ts->p.max_x = shyft::nan;
        auto qv=ts->values();
        for(size_t i=0;i<ts->size();++i)
            FAST_CHECK_EQ(qv[i], doctest::Approx(ev[i]));
        src.set(0,shyft::nan);
        src.set(1,shyft::nan);
        FAST_CHECK_UNARY(!isfinite(ts->value(0)));
        FAST_CHECK_UNARY(!isfinite(ts->value(1)));
        FAST_CHECK_UNARY(!isfinite(ts->value(2)));
        FAST_CHECK_UNARY(isfinite(ts->value(3)));
        FAST_CHECK_UNARY(isfinite(ts->value(4)));
        auto tsv=ts->values();
        FAST_CHECK_UNARY(!isfinite(tsv[0]));
        FAST_CHECK_UNARY(!isfinite(tsv[1]));
        FAST_CHECK_UNARY(!isfinite(tsv[2]));
        FAST_CHECK_UNARY(isfinite(tsv[3]));
        FAST_CHECK_UNARY(isfinite(tsv[4]));
        
        ts->cts=cts.ts;// now set in replacement values as time-series
        FAST_CHECK_EQ(ts->value(0),doctest::Approx(cts.value(0)));
        FAST_CHECK_EQ(ts->value(1),doctest::Approx(cts.value(1)));
        FAST_CHECK_EQ(ts->value(2),doctest::Approx(cts.value(2)));
        FAST_CHECK_EQ(ts->value(3),doctest::Approx(v[3]));// own value
        FAST_CHECK_EQ(ts->value(4),doctest::Approx(v[4]));// own value
        tsv=ts->values();
        FAST_CHECK_EQ(tsv[0],doctest::Approx(cts.value(0)));
        FAST_CHECK_EQ(tsv[1],doctest::Approx(cts.value(1)));
        FAST_CHECK_EQ(tsv[2],doctest::Approx(cts.value(2)));
        FAST_CHECK_EQ(tsv[3],doctest::Approx(v[3]));// own value
        FAST_CHECK_EQ(tsv[4],doctest::Approx(v[4]));// own value
        
        ts->p.min_x = 0.0;// clip out neg values
        FAST_CHECK_EQ(ts->value(4),doctest::Approx(cts.value(4)));// own value replaces -20.1
        tsv=ts->values();
        FAST_CHECK_EQ(tsv[4],doctest::Approx(cts.value(4)));// own value replaces -20.1
        ts->p.max_timespan=seconds(7);// the time-step is 10 seconds, so effectively, we have only valid values at specific time-points
        // NOW checking what happens when we examine the time-series (currently linear interpretation)
        // asking for the i'th value just returns that value at the beginning of interval (regardless, so it's valid).
        FAST_CHECK_EQ(ts->value(0),doctest::Approx(cts.value(0)));
        FAST_CHECK_EQ(ts->value(1),doctest::Approx(cts.value(1)));
        FAST_CHECK_EQ(ts->value(2),doctest::Approx(cts.value(2)));
        FAST_CHECK_EQ(ts->value(3),doctest::Approx(v[3]));// own value
        FAST_CHECK_EQ(ts->value(4),doctest::Approx(cts.value(4)));// own value replaces -20.1
        tsv=ts->values();
        FAST_CHECK_EQ(tsv[0],doctest::Approx(cts.value(0)));
        FAST_CHECK_EQ(tsv[1],doctest::Approx(cts.value(1)));
        FAST_CHECK_EQ(tsv[2],doctest::Approx(cts.value(2)));
        FAST_CHECK_EQ(tsv[3],doctest::Approx(v[3]));// own value
        FAST_CHECK_EQ(tsv[4],doctest::Approx(cts.value(4)));// own value replaces -20.1
        
        for(size_t i=0;i<ts->time_axis().size();++i) {
            FAST_CHECK_UNARY(isfinite(ts->value_at(ts->time_axis().time(i))));
            FAST_CHECK_UNARY(!isfinite(ts->value_at(ts->time_axis().time(i)+seconds(1))));
        }
        // switch to stair-case, and re-test:
        ts->set_point_interpretation(ts_point_fx::POINT_AVERAGE_VALUE);
        FAST_CHECK_EQ(ts->value(0),doctest::Approx(cts.value(0)));
        FAST_CHECK_EQ(ts->value(1),doctest::Approx(cts.value(1)));
        FAST_CHECK_EQ(ts->value(2),doctest::Approx(cts.value(2)));
        FAST_CHECK_EQ(ts->value(3),doctest::Approx(v[3]));// own value
        FAST_CHECK_EQ(ts->value(4),doctest::Approx(cts.value(4)));// own value replaces -20.1
        tsv=ts->values();
        FAST_CHECK_EQ(tsv[0],doctest::Approx(cts.value(0)));
        FAST_CHECK_EQ(tsv[1],doctest::Approx(cts.value(1)));
        FAST_CHECK_EQ(tsv[2],doctest::Approx(cts.value(2)));
        FAST_CHECK_EQ(tsv[3],doctest::Approx(v[3]));// own value
        FAST_CHECK_EQ(tsv[4],doctest::Approx(cts.value(4)));// own value replaces -20.1
        for(size_t i=0;i<ts->time_axis().size();++i) {
            FAST_CHECK_UNARY(isfinite(ts->value_at(ts->time_axis().time(i))));
            FAST_CHECK_UNARY(isfinite(ts->value_at(ts->time_axis().time(i)+seconds(7))));// YES! its valid for 7 seconds(inclusive)
            FAST_CHECK_UNARY(!isfinite(ts->value_at(ts->time_axis().time(i)+seconds(8))));// And, at the 8't second, it is invalid
        }
    }
    TEST_CASE("qac_ts_repeat_linear") {
        using shyft::nan;
        //                      *
        //                 0    1       2     3    4
        vector<double> vv {1.0,1.0,1.0,nan, 2.0};
        vector<double> cc {1.0,2.0,3.0,4.0, 5.0};
        vector<double> eav{1.0,1.0,1+1/3.0,1+2/3.0, 2.0};
        vector<double> ebv{1.0,1.0,3.0,4.0, 2.0};
        vector<double> ecv{1.0,1.0,-1.0,-1.0, 2.0};
        generic_dt ta{_t(0),seconds(10),vv.size()};

        apoint_ts src{ta,vv,ts_point_fx::POINT_INSTANT_VALUE};
        apoint_ts cts{ta,cc,ts_point_fx::POINT_INSTANT_VALUE};
        qac_parameter p1{seconds(30),0.0,9.0,seconds(10),0.1,-1.0};
        qac_parameter p2{seconds(0),0.0,9.0,seconds(10),0.1,-1.0};
        auto a = make_shared<qac_ts>(src,p1);//linear interpolation
        auto b = make_shared<qac_ts>(src,p1,cts);// replace from ts
        auto c = make_shared<qac_ts>(src,p2);// replace by fill = -1.0
        
        auto av=a->values();
        for(size_t i=0;i<av.size();++i) {
            FAST_CHECK_EQ(av[i],doctest::Approx(a->value(i)));
            FAST_CHECK_EQ(av[i],doctest::Approx(eav[i]));
        }
        auto bv=b->values();
        for(size_t i=0;i<bv.size();++i) {
            FAST_CHECK_EQ(bv[i],doctest::Approx(b->value(i)));
            FAST_CHECK_EQ(bv[i],doctest::Approx(ebv[i]));
        }
        auto cv=c->values();
        for(size_t i=0;i<cv.size();++i) {
            FAST_CHECK_EQ(cv[i],doctest::Approx(c->value(i)));
            FAST_CHECK_EQ(cv[i],doctest::Approx(ecv[i]));
        }

    }
    TEST_CASE("qac_ts_repeat_stair_case") {
        using shyft::nan;
        //                      *
        //                 0    1       2     3    4
        vector<double> vv {1.0,1.0,1.0,nan, 2.0};
        vector<double> cc {1.0,2.0,3.0,4.0, 5.0};
        vector<double> eav{1.0,1.0,nan,nan, 2.0};
        vector<double> ebv{1.0,1.0,3.0,4.0, 2.0};
        vector<double> ecv{1.0,1.0,-1.0,-1.0, 2.0};
        generic_dt ta{_t(0),seconds(10),vv.size()};

        apoint_ts src{ta,vv,ts_point_fx::POINT_AVERAGE_VALUE};
        apoint_ts cts{ta,cc,ts_point_fx::POINT_AVERAGE_VALUE};
        qac_parameter p1{seconds(30),0.0,9.0,seconds(10),0.1,-1.0};
        qac_parameter p2{seconds(0),0.0,9.0,seconds(10),0.1,-1.0};
        auto a = make_shared<qac_ts>(src,p1);//stair-case extrapolate interpolation
        auto b = make_shared<qac_ts>(src,p1,cts);// replace from ts
        auto c = make_shared<qac_ts>(src,p2);// replace by fill = -1.0
        
        auto av=a->values();
        for(size_t i=0;i<av.size();++i) {
            if(isfinite(eav[i])) {
                FAST_CHECK_EQ(av[i],doctest::Approx(a->value(i)));
                FAST_CHECK_EQ(av[i],doctest::Approx(eav[i]));
            } else {
                FAST_CHECK_UNARY(!isfinite(a->value(i)));
                FAST_CHECK_UNARY(!isfinite(av[i]));
            }
        }
        auto bv=b->values();
        for(size_t i=0;i<bv.size();++i) {
            FAST_CHECK_EQ(bv[i],doctest::Approx(b->value(i)));
            FAST_CHECK_EQ(bv[i],doctest::Approx(ebv[i]));
        }
        auto cv=c->values();
        for(size_t i=0;i<cv.size();++i) {
            FAST_CHECK_EQ(cv[i],doctest::Approx(c->value(i)));
            FAST_CHECK_EQ(cv[i],doctest::Approx(ecv[i]));
        }

    }
    TEST_CASE("qac_linear_fx") {
        using shyft::time_series::dd::qac::linear_fx;
        using shyft::core::utctime;
        using shyft::core::seconds;
        using shyft::core::to_seconds;
        double a=1.2;
        double b=10.0;
        utctime t0 = seconds(0);
        utctime t1 = seconds(1);
        double x0 = a*to_seconds(t0)+b;
        double x1 = a*to_seconds(t1)+b;
        
        for(utctime t=seconds(-10);t<seconds(10);t+=seconds(1)) {
            FAST_CHECK_EQ(linear_fx(t,t0,x0,t1,x1),doctest::Approx(a*to_seconds(t)+b) );
        }
    }
    TEST_CASE("qac_ts_repeating") {
        
        SUBCASE("No repetitions") {
            /* Test that in a series without any repetitions no values are replaced.
             */

            std::vector<double> values{0., 1., 3., 2., 4., 2., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
            generic_dt ta{ _t(0), seconds(10), values.size() };

            std::vector<double> expected = values;
        
            apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

            qac_parameter params = qac_parameter(seconds(0),shyft::nan,shyft::nan,seconds(60), 1e-2,100.0);
            auto ts = qac_ts(data_ts, params);
            auto result = ts.values();

            for ( std::size_t i = 0; i < values.size(); ++i ) {
                FAST_CHECK_EQ(result[i], expected[i]);
            }
        }
        SUBCASE("One repetition") {
            SUBCASE("Exact detection period span") {
                /* Test that a single repetition sequence is detected when the detection period is exactly as long as the repetition.
                 */

                //                                         *     *     *
                std::vector<double> values  {0., 1., 3.,   2.,   2.,   2., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
                std::vector<double> expected{0., 1., 3.,   2.,   2., 100., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params{ //= qac_parameter{//::create_repeating_constant_fill_parameters(
                    seconds(0),shyft::nan,shyft::nan, seconds(20)-utctime{1}, 1e-2, 100.
                };
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }

            SUBCASE("Shorter detection period span") {
                /* Test that a single repetition sequence is detected when the detection period is shorter than the repetition sequence.
                 */

                //                                         *     *     *
                std::vector<double> values  {0., 1., 3.,   2.,   2.,   2., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
                std::vector<double> expected{0., 1., 3.,   2., 100., 100., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params {//create_repeating_constant_fill_parameters(
                    seconds(0),shyft::nan,shyft::nan, seconds(9), 1e-2, 100.
                };
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }


            SUBCASE("Longer detection period span") {
                /* Test that a single repetition sequence is not detected when the detection period is longer than the repetition sequence.
                 */

                //                                       *   *   *
                std::vector<double> values  {0., 1., 3., 2., 2., 2., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
                std::vector<double> expected{0., 1., 3., 2., 2., 2., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params{// = qac_parameter::create_repeating_constant_fill_parameters(
                    seconds(0),shyft::nan,shyft::nan, seconds(40), 1e-2, 100.
                };
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }
        
            SUBCASE("All repetition") {
                /* Test that in a series with only repetition the entire series is flagged.
                 */

                std::vector<double> values  {  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.};
                std::vector<double> expected{ 0., 100., 100., 100., 100., 100., 100., 100., 100., 100., 100., 100.};
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params = qac_parameter{//create_repeating_constant_fill_parameters(
                    seconds(0),shyft::nan,shyft::nan, seconds(9), 1e-2, 100.
                };
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }
        }

        SUBCASE("Adjacent repetitions") {
            /* Test that in a series with adjacent repeting seqences values are replaced as expected.
             */

            //                                         *     *     *     *     *     *     *
            std::vector<double> values  {0., 1., 3.,   2.,   2.,   2.,   3.,   3.,   3.,   3., 4., 2., 5., 2., 1. };
            std::vector<double> expected{0., 1., 3.,   2.,   2., 100.,   3.,   3., 100., 100., 4., 2., 5., 2., 1. };
            generic_dt ta{ _t(0), seconds(10), values.size() };

            apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

            qac_parameter params = qac_parameter{//::create_repeating_constant_fill_parameters(
                seconds(0),shyft::nan, shyft::nan, seconds(20)-utctime{1}, 1e-2, 100.
            };
            auto ts = qac_ts(data_ts, params);
            auto result = ts.values();

            for ( std::size_t i = 0; i < values.size(); ++i ) {
                FAST_CHECK_EQ(result[i], expected[i]);
            }
        }
        SUBCASE("Leading repetition") {
            /* Test that a leading repetition sequence is detected.
             */

            //                             *     *     *     *     *     *
            std::vector<double> values  {  1.,   1.,   1.,   1.,   1.,   1., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
            std::vector<double> expected{  1.,   1.,   1., 100., 100., 100., 3., 4., 5., 7., 4., 2., 5., 2., 1. };
            generic_dt ta{ _t(0), seconds(10), values.size() };

            apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

            qac_parameter params = qac_parameter{//::create_repeating_constant_fill_parameters(
                seconds(0),shyft::nan,shyft::nan, seconds(30)-utctime{1}, 1e-2, 100.
            };
            auto ts = qac_ts(data_ts, params);
            auto result = ts.values();

            for ( std::size_t i = 0; i < values.size(); ++i ) {
                FAST_CHECK_EQ(result[i], expected[i]);
            }
        }

        SUBCASE("End repetition") {
            /* Test that a repetition sequence at the end of the series is detected.
             */

            //                                                                 *     *     *     *     *     *
            std::vector<double> values  {3., 4., 5., 7., 4., 2., 5., 2., 3.,   1.,   1.,   1.,   1.,   1.,   1. };
            std::vector<double> expected{3., 4., 5., 7., 4., 2., 5., 2., 3.,   1.,   1.,   1., 100., 100., 100. };
            generic_dt ta{ _t(0), seconds(10), values.size() };

            apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

            qac_parameter params = qac_parameter{//::create_repeating_constant_fill_parameters(
                seconds(0),shyft::nan,shyft::nan, seconds(30)-utctime{1}, 1e-2, 100.
            };
            auto ts = qac_ts(data_ts, params);
            auto result = ts.values();

            for ( std::size_t i = 0; i < values.size(); ++i ) {
                FAST_CHECK_EQ(result[i], expected[i]);
            }
        }

        SUBCASE("Repetition tolerance") {
            SUBCASE("Repetition around zero") {
                /* Test that repetitions around zero with a detection tolerance is replaced correctly.
                 */

                //                                            *     *     *     *     *
                std::vector<double> values  {1.1, 1.4, 2.5,  0.1, -0.2, -0.1,  0.0, -0.1, 2.5, 1.9, 2.1, 2.3, 3.1};
                std::vector<double> expected{1.1, 1.4, 2.5,  0.1, -0.2, -0.1, 100., 100., 2.5, 1.9, 2.1, 2.3, 3.1};
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params{// = qac_parameter::create_repeating_constant_fill_parameters(
                    seconds(0),shyft::nan,shyft::nan, seconds(30)-utctime{1}, 0.5, 100.
                };
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }

            SUBCASE("Repetition around value") {
                /* Test that repetitions around a value with a detection tolerance is replaced correctly.
                 */

                //                                            *     *     *     *     *
                std::vector<double> values  {1.1, 1.4, 2.5,  5.1,  4.9,  4.8,  5.0,  4.9, 2.5, 1.9, 2.1, 2.3, 3.1};
                std::vector<double> expected{1.1, 1.4, 2.5,  5.1,  4.9,  4.8, 100., 100., 2.5, 1.9, 2.1, 2.3, 3.1};
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params{// = qac_parameter::create_repeating_constant_fill_parameters(
                    seconds(0),shyft::nan,shyft::nan, seconds(30)-utctime{1}, 0.5, 100.
                };
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }

            SUBCASE("Repetition just inside tolerance") {
                /* Test that repetitions exactly within the tolerance threshold is replaced.
                 */

                const double eps = std::numeric_limits<double>::epsilon();

                std::vector<double> values  {1.1, 1.4, 2.5, 4.0 + 10*eps, 4.0 - 10*eps,  4.0 + 10*eps,  4.0 - 10*eps,  4.0 + 10*eps, 2.5, 1.9, 2.1, 2.3, 3.1};
                std::vector<double> expected{1.1, 1.4, 2.5, 4.0 + 10*eps,         100.,          100.,          100.,          100., 2.5, 1.9, 2.1, 2.3, 3.1};
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params {// qac_parameter::create_repeating_constant_fill_parameters(
                    seconds(0),shyft::nan,shyft::nan, seconds(10)-utctime{1}, 20*eps, 100.
                };
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }
            SUBCASE("Repetition just outside tolerance") {
                /* Test that repetitions exactly outside the tolerance threshold is not replaced.
                */

                const double eps = std::numeric_limits<double>::epsilon();

                std::vector<double> values  {1.1, 1.4, 2.5, 4.0 + 10*eps, 4.0 - 10*eps,  4.0 + 10*eps,  4.0 - 10*eps,  4.0 + 10*eps, 2.5, 1.9, 2.1, 2.3, 3.1};
                std::vector<double> expected{1.1, 1.4, 2.5, 4.0 + 10*eps, 4.0 - 10*eps,  4.0 + 10*eps,  4.0 - 10*eps,  4.0 + 10*eps, 2.5, 1.9, 2.1, 2.3, 3.1};
                generic_dt ta{ _t(0), seconds(10), values.size() };

                apoint_ts data_ts{ ta, values, ts_point_fx::POINT_AVERAGE_VALUE };

                qac_parameter params{// = qac_parameter::create_repeating_constant_fill_parameters(
                    seconds(0),shyft::nan,shyft::nan, seconds(10), 9*eps, 100.
                };
                auto ts = qac_ts(data_ts, params);
                auto result = ts.values();

                for ( std::size_t i = 0; i < values.size(); ++i ) {
                    FAST_CHECK_EQ(result[i], expected[i]);
                }
            }
        }
    }

}
