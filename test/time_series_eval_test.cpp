#include "test_pch.h"
#include "core/time_series_dd.h"

using namespace shyft::core;
using namespace shyft::time_series::dd;
using std::vector;
using std::make_shared;
using shyft::time_series::ts_point_fx;

/** helper function for testing */
static bool eval_ts_same_as_expr(apoint_ts expr) {
    using std::cout;
    eval_ctx c;
    c.evaluate(expr.ts);c.ref_count_mode=false;
    apoint_ts eval_expr(c.evaluate(expr.ts));
    if (eval_expr==expr)
        return true;
    if(expr.size() != eval_expr.size() ) {
        cout<<"differs:e.size()="<<expr.size()<<"eval.size()="<<eval_expr.size()<<"\n";
        return false;
    }
    if(expr.point_interpretation() != eval_expr.point_interpretation() ) {
        cout<<"differs in point interpretation\n";
        return false;
    }
    if(expr.time_axis() != eval_expr.time_axis()) {
        cout<<".time_axis differs\n";
        return false;
        
    }
    if(expr.values()!=eval_expr.values()) {
        cout<<".values differs\n";
        return false;
    }
    cout<<"something else differs\n";
    for(size_t i=0;i<expr.size();++i) {
        auto a=expr.value(i);
        auto b=eval_expr.value(i);
        if(fabs(a-b)>1e-8)
            cout<<"the "<<i<<"th value differs: "<<a<<"!=" <<b<<"\n";
    }
    return false;
    
    
}

TEST_SUITE("time_series") {
    
TEST_CASE("eval_ctx") {
    eval_ctx c;
    FAST_CHECK_EQ(c.ref_count_mode,true);
    FAST_CHECK_EQ(c.evaluated.size(),0u);
    FAST_CHECK_EQ(c.ref_count.size(),0u);
    calendar utc;
    auto t0=utc.time(2018,1,1);
    auto dt=seconds(3600);
    size_t n=5;
    apoint_ts a(gta_t(t0,dt,n),vector<double>(n,1.0),ts_point_fx::POINT_AVERAGE_VALUE);
    apoint_ts b(gta_t(t0,dt,n),vector<double>(n,2.0),ts_point_fx::POINT_AVERAGE_VALUE);
    auto aa= a/b;
    auto bb = b/a;
    auto e= a+b+aa+bb;
    auto f= 2*e+2*bb;
    auto d= 2*e + f;
    c.evaluate(d.ts);
    FAST_CHECK_EQ(c.ref_count.size(),10);// a/b ,b/a, a+b, a+b+aa,a+b+aa+b, 2*e, 2*b, 2*e+2*bb, 2*e, 2*e+f
    //FAST_CHECK_EQ(c.ref_count[a.ts.get()],2);
    //FAST_CHECK_EQ(c.ref_count[b.ts.get()],2);
    FAST_CHECK_EQ(c.ref_count[e.ts.get()],2);
    FAST_CHECK_EQ(c.ref_count[f.ts.get()],1);
    FAST_CHECK_EQ(c.ref_count[d.ts.get()],1);
    FAST_CHECK_EQ(c.ref_count[aa.ts.get()],1);// ensure that aa is only ref'd once,since it's a sibling of e only
    FAST_CHECK_EQ(c.ref_count[bb.ts.get()],2);// ensure that bb is 2, since it's ref'd from e, and also through f (bypassing e)
    c.ref_count_mode=false;
    auto d_v=c.evaluate(d.ts);
    //FAST_CHECK_EQ(c.is_evaluated(a.ts.get()),true);
    //FAST_CHECK_EQ(c.is_evaluated(b.ts.get()),true);
    FAST_CHECK_EQ(c.is_evaluated(e.ts.get()),true);// because it's ref'd twice,we remember it
    FAST_CHECK_EQ(c.is_evaluated(bb.ts.get()),true);// because it's ref'd twice,we remember it
    FAST_CHECK_EQ(c.is_evaluated(f.ts.get()),false);
    FAST_CHECK_EQ(c.is_evaluated(d.ts.get()),false);
    FAST_CHECK_EQ(c.is_evaluated(aa.ts.get()),false);
    
}

TEST_CASE("ts_evaluate") {
    calendar utc;
    auto t0=utc.time(2018,1,1);
    auto dt=seconds(3600);
    size_t n=5*24;
    gta_t tah(t0,dt,n);
    gta_t tad(t0,dt*24,n/24);
    apoint_ts a(tah,vector<double>(n,1.0),ts_point_fx::POINT_AVERAGE_VALUE);
    apoint_ts b(tah,vector<double>(n,2.0),ts_point_fx::POINT_AVERAGE_VALUE);
    double c=3.0;
    auto e = a+4*b;
                        
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a+b));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a-b));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a*b));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a/b));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a.pow(b)));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a.max(b)));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a.min(b)));
    
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(c+b));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(c-b));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(c*b));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(c/b));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a.pow(c)));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a.max(c)));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a.min(c)));

    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a+c));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a-c));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a*c));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a/c));
    
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.average(tad)));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.accumulate(tad)));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.integral(tad)));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.derivative()));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.time_shift(dt)));

    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.log()));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.abs()));
    
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.convolve_w(vector<double>{0.1,0.2,0.4,0.2,0.1},
                                                         shyft::time_series::convolve_policy::USE_NEAREST |
                                                         shyft::time_series::convolve_policy::BACKWARD)));
    {
        using shyft::time_series::rating_curve_function;
        using shyft::time_series::rating_curve_segment;
        using shyft::time_series::rating_curve_parameters;
        std::array<std::pair<utctime, rating_curve_function>, 3> curve_data{
            std::make_pair(t0+deltahours(1u),  std::vector<rating_curve_segment>{
                rating_curve_segment{  0., 1., 0., 1. },
                rating_curve_segment{  4., 2., 0., 1. }
            }),

            std::make_pair(t0+deltahours(9u),  std::vector<rating_curve_segment>{
                rating_curve_segment{  0., 3., 0., 1. },
                rating_curve_segment{  5., 4., 0., 1. },
                rating_curve_segment{ 10., 5., 0., 1. }
            }),

            std::make_pair(t0+deltahours(20u), std::vector<rating_curve_segment>{
                rating_curve_segment{  0., 6., 0., 1. },
                rating_curve_segment{  6., 7., 0., 1. },
                rating_curve_segment{  9., 8., 0., 1. }
            })
        };

        rating_curve_parameters rcp{ curve_data.cbegin(), curve_data.cend() };
        FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.rating_curve(rcp)));

    }
    {
        shyft::time_series::ice_packing_parameters p;
        auto ice_ts=e.ice_packing(p,shyft::time_series::ice_packing_temperature_policy::ALLOW_ANY_MISSING);
        FAST_CHECK_EQ(true,eval_ts_same_as_expr(ice_ts));
        ice_packing_recession_parameters rp;        
        FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.ice_packing_recession(ice_ts,rp)));
    }
    qac_parameter qc1;
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.quality_and_self_correction(qc1)));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.quality_and_ts_correction(qc1,a)));
 
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.inside(0.5,2.8,shyft::nan,1.0,0.0)));

    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.decode(0,1)));
    
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(b.merge_points(a)));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.use_time_axis_from(a.average(tad))));
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(a.slice(2,10)));
    
    FAST_CHECK_EQ(true,eval_ts_same_as_expr(e.krls_interpolation(seconds(3600),1e-3,0.1,100)));

    {
        vector<double> bucket_data{ 0, 0, 0, 0, 0, 0,
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
        gta_t tab(seconds(0), calendar::HOUR, bucket_data.size());
        apoint_ts bts_raw(tab, bucket_data, ts_point_fx::POINT_AVERAGE_VALUE);
        FAST_CHECK_EQ(true,eval_ts_same_as_expr(bts_raw.bucket_to_hourly(0,-100.0)));

    }
    
}

}
