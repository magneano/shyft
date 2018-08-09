#include "test_pch.h"

#include "core/dtss_db.h"

#include <boost/filesystem.hpp>
#include <cstdio>
namespace  fs=boost::filesystem;

using namespace shyft::dtss;
using namespace shyft::core;
using namespace shyft;
using shyft::time_series::POINT_AVERAGE_VALUE;
using std::exception;
using std::fopen;
using std::fclose;

#define TEST_SECTION(x) 
//{std::cout<<"section:"<<x<<"\n";}
static inline shyft::core::utctime _t(int64_t t1970s) {return shyft::core::utctime{shyft::core::seconds(t1970s)};}

template <class T>
void test_time_io() {
        auto fpath = fs::unique_path();
		string cp = fpath.string();
        auto *fp=fopen(cp.c_str(),"wb");
        vector<utctime> tv; tv.push_back(seconds(1));tv.push_back(seconds(2));
        utctime t{seconds(4)};
        T::write(fp,tv);
        T::write(fp,t);
        fclose(fp);
        fp=fopen(cp.c_str(),"rb");
        vector<utctime> tvr;tvr.resize(tv.size());
        utctime tr;
        T::read(fp,tvr);
        T::read(fp,tr);
        fclose(fp);
        fs::remove(fpath);
        FAST_CHECK_EQ(tv,tvr);
        FAST_CHECK_EQ(t,tr);
}

void test_dtss_db_basics(bool time_in_micro_seconds) {
        using namespace shyft::dtss;
        using namespace shyft::time_series::dd;
        using time_series::point_ts;

        std::shared_ptr<core::calendar> utc = std::make_shared<core::calendar>();
        std::shared_ptr<core::calendar> osl = std::make_shared<core::calendar>("Europe/Oslo");

        core::utctime t = utc->time(2016, 1, 1);
        core::utctimespan dt = core::deltahours(1);
        core::utctimespan dt_half = core::deltaminutes(30);
        std::size_t n = 24 * 365 * 2;//24*365*5;

        // construct time-axis that we want to test.
        time_axis::fixed_dt fta(t, dt, n);
        time_axis::calendar_dt cta1(utc,t,dt,n);
        time_axis::calendar_dt cta2(osl,t,dt,n);

        vector<utctime> tp;for(std::size_t i=0;i<fta.size();++i)tp.push_back(fta.time(i));
        time_axis::point_dt pta(tp,fta.total_period().end);
        auto dirname = "ts.db.test." +(time_in_micro_seconds?string(".s."):string(".us."))+ std::to_string(to_seconds64(core::utctime_now()) );
        auto tmpdir = (fs::temp_directory_path()/dirname);
        ts_db db(tmpdir.string());
        db.time_format_micro_seconds=time_in_micro_seconds;
        TEST_SECTION("store_fixed_dt") {
            gts_t o(gta_t(fta),10.0,time_series::ts_point_fx::POINT_AVERAGE_VALUE);
            o.set(0, 100.); o.set(o.size()-1, 100.);
            std::string fn("measurements/tssf.db");  // verify we can have path-parts
            db.save(fn,o);

            // read all
            auto r1 = db.read(fn, utcperiod{});
            FAST_CHECK_EQ(o.point_interpretation(), r1.point_interpretation());
            FAST_CHECK_EQ(o.time_axis(), r1.time_axis());

            // read inner slice
            core::utctime tb = t + 3u*dt/2u;
            core::utctime te = t + (2u*n - 3u)*dt/2u;
            auto r2 = db.read(fn, utcperiod{ tb, te });
            FAST_CHECK_EQ(o.point_interpretation(), r2.point_interpretation());
            FAST_CHECK_EQ(r2.time_axis(), time_axis::generic_dt(t + dt, dt, n - 2u));
            FAST_CHECK_EQ(r2.value(0), o.value(1));  // dropped first value of o
            FAST_CHECK_EQ(r2.value(r2.size() - 1), o.value(o.size() - 2));  // dropped last value of o

            auto fr = db.find(string("measurements/.*\\.db")); // should match our ts.
            FAST_CHECK_EQ(fr.size(), 1 );

            db.remove(fn);
            fr = db.find(string("measurements/.*\\.db")); // should match our ts.
            FAST_CHECK_EQ(fr.size(),0);
        }

        TEST_SECTION("store_calendar_utc_dt") {
            gts_t o(gta_t(cta1), 10.0, time_series::ts_point_fx::POINT_AVERAGE_VALUE);
            o.set(0, 100.); o.set(o.size()-1, 100.);
            std::string fn("tssf1.db");
            db.save(fn, o);

            // read all
            auto r=db.read(fn, utcperiod{ });
            FAST_CHECK_EQ(o.point_interpretation(), r.point_interpretation());
            FAST_CHECK_EQ(o.time_axis(), r.time_axis());

            // read inner slice
            core::utctime tb = cta1.cal->add(t, dt_half, 3);
            core::utctime te = cta1.cal->add(t, dt_half, 2*n - 3);
            auto r2 = db.read(fn, utcperiod{ tb, te });

            FAST_CHECK_EQ(o.point_interpretation(), r2.point_interpretation());
            FAST_CHECK_EQ(r2.time_axis(), time_axis::generic_dt{ cta1.cal, cta1.cal->add(t, dt, 1), dt, n - 2u });
            FAST_CHECK_EQ(r2.value(0), o.value(1));  // dropped first value of o
            FAST_CHECK_EQ(r2.value(r2.size() - 1), o.value(o.size() - 2));  // dropped last value of o

            auto i = db.get_ts_info(fn);
            FAST_CHECK_EQ(i.name,fn);
            FAST_CHECK_EQ(i.data_period,o.total_period());
            FAST_CHECK_EQ(i.point_fx, o.point_interpretation());
            FAST_CHECK_LE( i.modified, utctime_now());

            auto fr = db.find(string(".ss.1\\.db")); // should match our ts.
            FAST_CHECK_EQ(fr.size(), 1 );

            db.remove(fn);
            try {
                auto rx=db.read(fn+".not.there", utcperiod{});
                FAST_CHECK_UNARY(rx.size()==3);
            } catch(const exception&) {
                FAST_CHECK_UNARY(true);
            }

        }

        TEST_SECTION("store_calendar_osl_dt") {
            gts_t o(gta_t(cta2),10.0,time_series::ts_point_fx::POINT_AVERAGE_VALUE);
            string fn("tssf2.db");
            db.save(fn,o);
            auto r=db.read(fn,utcperiod{});
            FAST_CHECK_EQ(o.point_interpretation(),r.point_interpretation());
            FAST_CHECK_EQ(o.time_axis(),r.time_axis());
            db.remove(fn);
        }
        TEST_SECTION("store_point_dt") {
            gts_t o(gta_t(pta), 10.0, time_series::ts_point_fx::POINT_INSTANT_VALUE);
            o.set(0, 100.); o.set(o.size()-1, 100.);
            string fn("tssf3.db");
            db.save(fn, o);
            // read all
            auto r = db.read(fn, utcperiod{});
            FAST_CHECK_EQ(o.point_interpretation(), r.point_interpretation());
            FAST_CHECK_EQ(o.time_axis(), r.time_axis());

            // read slice
            auto tb = _t( (to_seconds64(pta.time(2)) + to_seconds64(pta.time(1)) ) / 2 );
            auto te = _t( (to_seconds64(pta.time(pta.size() - 1)) + to_seconds64(pta.time(pta.size() - 2)) ) / 2);
            auto r2 = db.read(fn, utcperiod{ tb, te });
            FAST_CHECK_EQ(o.point_interpretation(), r2.point_interpretation());
            time_axis::point_dt exp{ std::vector<core::utctime>( &o.ta.p.t[1], &o.ta.p.t[o.ta.p.t.size()-1] ), o.ta.p.t[o.ta.p.t.size()-1] };
            FAST_CHECK_EQ(r2.time_axis(), gta_t(exp));
            FAST_CHECK_EQ(r2.value(0), o.value(1));  // dropped first value of o
            FAST_CHECK_EQ(r2.value(r2.size() - 1), o.value(o.size() - 2));  // dropped last value of o

            auto i = db.get_ts_info(fn);
            FAST_CHECK_EQ(i.name,fn);
            FAST_CHECK_EQ(i.data_period,o.total_period());
            FAST_CHECK_EQ(i.point_fx, o.point_interpretation());
            FAST_CHECK_LE( i.modified, utctime_now());
            db.remove(fn);
        }

        TEST_SECTION("dtss_db_speed") {
            int n_ts = 120;
            vector<gts_t> tsv; tsv.reserve(n_ts);
            double fv = 1.0;
            for (int i = 0; i < n_ts; ++i)
                tsv.emplace_back(gta_t(fta), fv += 1.0,shyft::time_series::ts_point_fx::POINT_AVERAGE_VALUE);
            FAST_CHECK_EQ(n_ts,tsv.size());

            auto t0 = timing::now();
            for(std::size_t i=0;i<tsv.size();++i) {
                std::string fn("ts."+std::to_string(i)+".db");
                db.save(fn,tsv[i]);
            }
            auto t1 = timing::now();
            vector<gts_t> rv;
            for(std::size_t i=0;i<tsv.size();++i) {
                std::string fn("ts."+std::to_string(i)+".db");
                rv.push_back(db.read(fn,utcperiod{}));
            }
            auto t2= timing::now();
            auto w_mb_s= n_ts*n/double(elapsed_ms(t0,t1))/1000.0;
            auto r_mb_s= n_ts*n/double(elapsed_ms(t1,t2))/1000.0;
            // on windows(before workaround): ~ 6 mpts/sec write, 162 mpts/sec read (slow close->workaround with thread?)
            // on linux: ~ 120 mpts/sec write, 180 mpts/sec read
            std::cout<<"write Mpts/s = "<<w_mb_s<<", read Mpts/s = "<<r_mb_s<<" pts = "<<n_ts*n<<", roundtrip ms="<< double(elapsed_ms(t0,t2)) <<"\n";
            //std::cout << "open_ms:" << db.t_open << ", write_ms:" << db.t_write << ", t_close_ms:" << db.t_close << std::endl;
            FAST_CHECK_EQ(rv.size(),tsv.size());
            //fs::remove_all("*.db");
        }
#ifdef _WIN32
        for (int i = 0; i<10; ++i) {
            std::this_thread::sleep_for(std::chrono::duration<int, std::milli>(1000));
            try {
                fs::remove_all(tmpdir);
                break;
            }
            catch (...) {
                std::cout <<"Try #"<<i+1<< ":Failed to remove " << tmpdir << "\n";
            }
        }
#else
        fs::remove_all(tmpdir);
#endif
   
}

void test_dtss_db_store_merge_write(bool time_in_micro_seconds ) {
    namespace core = shyft::core;
    namespace dtss = shyft::dtss;
    namespace ta = shyft::time_axis;
    namespace ts = shyft::time_series;
    using shyft::time_series::dd::gta_t;
    // setup db
    auto dirname = "ts.db.test." + std::to_string(to_seconds64(core::utctime_now()));
    auto tmpdir = (fs::temp_directory_path()/dirname);
    dtss::ts_db db(tmpdir.string());
    db.time_format_micro_seconds=time_in_micro_seconds;
    SUBCASE("error_handling") {
        SUBCASE("extending with different ta") {
            // data
            std::shared_ptr<core::calendar> utc_ptr = std::make_shared<core::calendar>();
            const core::utctime t0 = core::utctime_now();
            const core::utctimespan dt_h = core::calendar::HOUR;
            const core::utctimespan dt_d = core::calendar::DAY;
            const std::size_t n = 1000;
            // -----
            ta::fixed_dt f_ta_h{ t0, dt_h, n };
            ta::calendar_dt c_ta_d{ utc_ptr, t0, dt_d, n };
            ts::point_ts<ta::generic_dt> pts_h{ gta_t(f_ta_h), 0. };
            ts::point_ts<ta::generic_dt> pts_d{ gta_t(c_ta_d), 0. };
            // -----
            std::string fn("dtss_save_merge/ext_diff_ta.db");

            // save initital data
            db.save(fn, pts_d, false);
            auto find_res = db.find(string("dtss_save_merge/ext_diff_ta\\.db"));
            FAST_CHECK_EQ(find_res.size(), 1);

            // add data to the same path
            CHECK_THROWS_AS_MESSAGE(db.save(fn, pts_h, false),
                std::runtime_error, "dtss_store: cannot merge with different ta type");

            // cleanup
            db.remove(fn);
            find_res = db.find(string("dtss_save_merge/ext_diff_ta\\.db")); // should match our ts.
            FAST_CHECK_EQ(find_res.size(), 0);
        }
        SUBCASE("extending with different point interpretation") {
            // data
            std::shared_ptr<core::calendar> utc_ptr = std::make_shared<core::calendar>();
            const core::utctime t0 = core::utctime_now();
            const core::utctimespan dt_h = core::calendar::HOUR;
            const core::utctimespan dt_d = core::calendar::DAY;
            const std::size_t n = 1000;
            // -----
            ta::fixed_dt f_ta_h{ t0, dt_h, n };
            ta::fixed_dt f_ta_d{ t0, dt_d, n };
            ts::point_ts<ta::generic_dt> pts_h{ gta_t(f_ta_h), 0., ts::POINT_INSTANT_VALUE };
            ts::point_ts<ta::generic_dt> pts_d{ gta_t(f_ta_d), 0., ts::POINT_AVERAGE_VALUE };
            // -----
            std::string fn("dtss_save_merge/ext_diff_fx.db");

            // save initital data
            db.save(fn, pts_d, false);
            auto find_res = db.find(string("dtss_save_merge/ext_diff_fx\\.db"));
            FAST_REQUIRE_EQ(find_res.size(), 1);

            // add data to the same path
            CHECK_THROWS_AS_MESSAGE(db.save(fn, pts_h, false),
                std::runtime_error, "dtss_store: cannot merge with different point interpretation");

            // cleanup
            db.remove(fn);
            find_res = db.find(string("dtss_save_merge/ext_diff_fx\\.db")); // should match our ts.
            FAST_CHECK_EQ(find_res.size(), 0);
        }
        SUBCASE("fixed_dt old_dt != new_dt") {
            // data
            const core::utctime t0 = core::utctime_now();
            const core::utctimespan dt_h = core::calendar::HOUR;
            const core::utctimespan dt_d = core::calendar::DAY;
            const std::size_t n = 1000;
            // -----
            ta::fixed_dt f_ta_h{ t0, dt_h, n };
            ta::fixed_dt f_ta_d{ t0, dt_d, n };
            ts::point_ts<ta::generic_dt> pts_h{ gta_t(f_ta_h), 0. };
            ts::point_ts<ta::generic_dt> pts_d{ gta_t(f_ta_d), 0. };
            // -----
            std::string fn("dtss_save_merge/ext_fixed_diff_dt.db");  // verify we can have path-parts

            // save initital data
            db.save(fn, pts_d, false);
            auto find_res = db.find(string("dtss_save_merge/ext_fixed_diff_dt\\.db"));
            FAST_CHECK_EQ(find_res.size(), 1);

            // add data to the same path
            CHECK_THROWS_AS_MESSAGE(db.save(fn, pts_h, false),
                std::runtime_error, "dtss_store: cannot merge unaligned fixed_dt");

            // cleanup
            db.remove(fn);
            find_res = db.find(string("dtss_save_merge/ext_fixed_diff_dt\\.db")); // should match our ts.
            FAST_CHECK_EQ(find_res.size(), 0);
        }
        SUBCASE("fixed_dt unaligned axes") {
            // data
            core::calendar utc{};
            const core::utctime t0_1 = utc.time(2000, 1, 1, 0, 0);
            const core::utctime t0_2 = utc.time(2000, 1, 1, 0, 13);  // 13 minutes shifted
            const core::utctimespan dt_h = core::calendar::HOUR;
            const std::size_t n = 1000;
            // -----
            ta::fixed_dt f_ta_h{ t0_1, dt_h, n };
            ta::fixed_dt f_ta_d{ t0_2, dt_h, n };
            ts::point_ts<ta::generic_dt> pts_h{ gta_t(f_ta_h), 0. };
            ts::point_ts<ta::generic_dt> pts_d{ gta_t(f_ta_d), 0. };
            // -----
            std::string fn("dtss_save_merge/ext_fixed_unaligned.db");  // verify we can have path-parts

                                                     // save initital data
            db.save(fn, pts_d, false);
            auto find_res = db.find(string("dtss_save_merge/ext_fixed_unaligned\\.db"));
            FAST_CHECK_EQ(find_res.size(), 1);

            // add data to the same path
            CHECK_THROWS_AS_MESSAGE(db.save(fn, pts_h, false),
                std::runtime_error, "dtss_store: cannot merge unaligned fixed_dt");

            // cleanup
            db.remove(fn);
            find_res = db.find(string("dtss_save_merge/ext_fixed_unaligned\\.db")); // should match our ts.
            FAST_CHECK_EQ(find_res.size(), 0);
        }
        SUBCASE("calendar_dt old_dt != new_dt") {
            // data
            std::shared_ptr<core::calendar> utc_ptr = std::make_shared<core::calendar>();
            const core::utctime t0 = core::utctime_now();
            const core::utctimespan dt_h = core::calendar::HOUR;
            const core::utctimespan dt_d = core::calendar::DAY;
            const std::size_t n = 1000;
            // -----
            ta::calendar_dt f_ta_h{ utc_ptr, t0, dt_h, n };
            ta::calendar_dt f_ta_d{ utc_ptr, t0, dt_d, n };
            ts::point_ts<ta::generic_dt> pts_h{ gta_t(f_ta_h), 0. };
            ts::point_ts<ta::generic_dt> pts_d{ gta_t(f_ta_d), 0. };
            // -----
            std::string fn("dtss_save_merge/ext_cal_diff_dt.db");  // verify we can have path-parts

           // save initital data
            db.save(fn, pts_d, false);
            auto find_res = db.find(string("dtss_save_merge/ext_cal_diff_dt\\.db"));
            FAST_CHECK_EQ(find_res.size(), 1);

            // add data to the same path
            CHECK_THROWS_AS_MESSAGE(db.save(fn, pts_h, false),
                std::runtime_error, "dtss_store: cannot merge unaligned calendar_dt");

            // cleanup
            db.remove(fn);
            find_res = db.find(string("dtss_save_merge/ext_cal_diff_dt\\.db")); // should match our ts.
            FAST_CHECK_EQ(find_res.size(), 0);
        }
        SUBCASE("calendar_dt unaligned axes") {
            // data
            std::shared_ptr<core::calendar> utc_ptr = std::make_shared<core::calendar>();
            const core::utctime t0_1 = utc_ptr->time(2000, 1, 1, 0, 0);
            const core::utctime t0_2 = utc_ptr->time(2000, 1, 1, 0, 13);  // 13 minutes shifted
            const core::utctimespan dt_h = core::calendar::HOUR;
            const std::size_t n = 1000;
            // -----
            ta::calendar_dt f_ta_h{ utc_ptr, t0_1, dt_h, n };
            ta::calendar_dt f_ta_d{ utc_ptr, t0_2, dt_h, n };
            ts::point_ts<ta::generic_dt> pts_h{ gta_t(f_ta_h), 0. };
            ts::point_ts<ta::generic_dt> pts_d{ gta_t(f_ta_d), 0. };
            // -----
            std::string fn("dtss_save_merge/ext_cal_unaligned.db");  // verify we can have path-parts

            // save initital data
            db.save(fn, pts_d, false);
            auto find_res = db.find(string("dtss_save_merge/ext_cal_unaligned\\.db"));
            FAST_CHECK_EQ(find_res.size(), 1);

            // add data to the same path
            CHECK_THROWS_AS_MESSAGE(db.save(fn, pts_h, false),
                std::runtime_error, "dtss_store: cannot merge unaligned calendar_dt");

            // cleanup
            db.remove(fn);
            find_res = db.find(string("dtss_save_merge/ext_cal_unaligned\\.db")); // should match our ts.
            FAST_CHECK_EQ(find_res.size(), 0);
        }
        SUBCASE("calendar_dt different calendars") {
            // data
            std::shared_ptr<core::calendar> utc_ptr = std::make_shared<core::calendar>();
            std::shared_ptr<core::calendar> osl_ptr = std::make_shared<core::calendar>("Europe/Oslo");
            const core::utctime t0_1 = utc_ptr->time(2000, 1, 1, 0, 0);
            const core::utctime t0_2 = osl_ptr->time(2000, 1, 1, 0, 13);  // 13 minutes shifted
            const core::utctimespan dt_h = core::calendar::HOUR;
            const std::size_t n = 1000;
            // -----
            ta::calendar_dt f_ta_h{ utc_ptr, t0_1, dt_h, n };
            ta::calendar_dt f_ta_d{ osl_ptr, t0_2, dt_h, n };
            ts::point_ts<ta::generic_dt> pts_h{ gta_t(f_ta_h), 0. };
            ts::point_ts<ta::generic_dt> pts_d{ gta_t(f_ta_d), 0. };
            // -----
            std::string fn("dtss_save_merge/ext_cal_diff_cal.db");  // verify we can have path-parts

            // save initital data
            db.save(fn, pts_d, false);
            auto find_res = db.find(string("dtss_save_merge/ext_cal_diff_cal\\.db"));
            FAST_CHECK_EQ(find_res.size(), 1);

            // add data to the same path
            CHECK_THROWS_AS_MESSAGE(db.save(fn, pts_h, false),
                std::runtime_error, "dtss_store: cannot merge calendar_dt with different calendars");

            // cleanup
            db.remove(fn);
            find_res = db.find(string("dtss_save_merge/ext_cal_diff_cal\\.db")); // should match our ts.
            FAST_CHECK_EQ(find_res.size(), 0);
        }
    }
    SUBCASE("merging time-series") {
        std::shared_ptr<core::calendar> utc_ptr = std::make_shared<core::calendar>();
        std::shared_ptr<core::calendar> osl_ptr = std::make_shared<core::calendar>("Europe/Oslo");

        SUBCASE("fixed_dt") {
            SUBCASE("exact") {
                // data
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0 = utc_ptr->time(2016, 1, 1);
                // -----
                ta::fixed_dt f_ta{ t0, dt, n };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(f_ta), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(f_ta), 10. };
                // -----
                std::string fn("dtss_save_merge/fixed_exact.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/fixed_exact\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/fixed_exact.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::FIXED);
                FAST_CHECK_EQ(res.total_period().start, f_ta.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, f_ta.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 10.);
                FAST_CHECK_EQ(res.v.at(n - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/fixed_exact\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new contained in old") {
                // data
                const std::size_t drop = 10u;  // points to drop from start/end of old
                // -----
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old + dt*drop;
                // -----
                ta::fixed_dt f_ta_old{ t0_old, dt, n };
                ta::fixed_dt f_ta_new{ t0_new, dt, n - 2*drop };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(f_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(f_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/fixed_new_in_old.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/fixed_new_in_old\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/fixed_new_in_old.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::FIXED);
                FAST_CHECK_EQ(res.total_period().start, f_ta_old.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, f_ta_old.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 1.);
                FAST_CHECK_EQ(res.v.at(drop - 1), 1.);
                FAST_CHECK_EQ(res.v.at(drop), 10.);
                FAST_CHECK_EQ(res.v.at(n - drop - 1), 10.);
                FAST_CHECK_EQ(res.v.at(n - drop), 1.);
                FAST_CHECK_EQ(res.v.at(n - 1), 1.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/fixed_new_in_old\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("old contained in new") {
                // data
                const std::size_t extra = 10u;  // points to drop from start/end of old
                // -----
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old - dt*extra;
                // -----
                ta::fixed_dt f_ta_old{ t0_old, dt, n };
                ta::fixed_dt f_ta_new{ t0_new, dt, n + extra };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(f_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(f_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/fixed_old_in_new.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/fixed_old_in_new\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/fixed_old_in_new.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::FIXED);
                FAST_CHECK_EQ(res.total_period().start, f_ta_new.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, f_ta_new.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 10.);
                FAST_CHECK_EQ(res.v.at(n - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/fixed_old_in_new\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new overlap start of old") {
                // data
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old - dt*n/2;
                // -----
                ta::fixed_dt f_ta_old{ t0_old, dt, n };
                ta::fixed_dt f_ta_new{ t0_new, dt, n };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(f_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(f_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/fixed_new_over_start.db");

                FAST_REQUIRE_LT(t0_new, t0_old);

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/fixed_new_over_start\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/fixed_new_over_start.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::FIXED);
                FAST_CHECK_EQ(res.total_period().start, f_ta_new.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, f_ta_old.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 10.);
                FAST_CHECK_EQ(res.v.at(n - 1), 10.);
                FAST_CHECK_EQ(res.v.at(n), 1.);
                FAST_CHECK_EQ(res.v.at(n + n/2 - 1), 1.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/fixed_new_over_start\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new overlap end of old") {
                // data
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old + dt*n/2;
                // -----
                ta::fixed_dt f_ta_old{ t0_old, dt, n };
                ta::fixed_dt f_ta_new{ t0_new, dt, n };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(f_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(f_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/fixed_new_over_end.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/fixed_new_over_end\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/fixed_new_over_end.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::FIXED);
                FAST_CHECK_EQ(res.total_period().start, f_ta_old.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, f_ta_new.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 1.);
                FAST_CHECK_EQ(res.v.at(n/2 - 1), 1.);
                FAST_CHECK_EQ(res.v.at(n/2), 10.);
                FAST_CHECK_EQ(res.v.at(n + n/2 - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/fixed_new_over_end\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("consecutive without gap") {
                // data
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old + dt*n;
                // -----
                ta::fixed_dt f_ta_old{ t0_old, dt, n };
                ta::fixed_dt f_ta_new{ t0_new, dt, n };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(f_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(f_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/fixed_consec.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/fixed_consec\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/fixed_consec.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::FIXED);
                FAST_CHECK_EQ(res.total_period().start, f_ta_old.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, f_ta_new.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 1.);
                FAST_CHECK_EQ(res.v.at(n - 1), 1.);
                FAST_CHECK_EQ(res.v.at(n), 10.);
                FAST_CHECK_EQ(res.v.at(2 * n - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/fixed_consec\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new after old with gap") {
                // data
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old + 2*n*dt;
                // -----
                ta::fixed_dt f_ta_old{ t0_old, dt, n };
                ta::fixed_dt f_ta_new{ t0_new, dt, n };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(f_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(f_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/fixed_gap_after.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/fixed_gap_after\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/fixed_gap_after.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::FIXED);
                FAST_CHECK_EQ(res.ta.size(), 3*n);
                FAST_CHECK_EQ(res.total_period().start, f_ta_old.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, f_ta_new.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 1.);
                FAST_CHECK_EQ(res.v.at(n - 1), 1.);
                FAST_CHECK_UNARY(std::isnan(res.v.at(n)));
                FAST_CHECK_UNARY(std::isnan(res.v.at(2*n - 1)));
                FAST_CHECK_EQ(res.v.at(2 * n), 10.);
                FAST_CHECK_EQ(res.v.at(3 * n - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/fixed_gap_after\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new before old with gap") {
                // data
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old - 2*n*dt;
                // -----
                ta::fixed_dt f_ta_old{ t0_old, dt, n };
                ta::fixed_dt f_ta_new{ t0_new, dt, n };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(f_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(f_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/fixed_gap_before.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/fixed_gap_before\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/fixed_gap_before.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::FIXED);
                FAST_CHECK_EQ(res.ta.size(), 3*n);
                FAST_CHECK_EQ(res.total_period().start, f_ta_new.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, f_ta_old.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 10.);
                FAST_CHECK_EQ(res.v.at(n - 1), 10.);
                FAST_CHECK_UNARY(std::isnan(res.v.at(n)));
                FAST_CHECK_UNARY(std::isnan(res.v.at(2*n - 1)));
                FAST_CHECK_EQ(res.v.at(2 * n), 1.);
                FAST_CHECK_EQ(res.v.at(3 * n - 1), 1.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/fixed_gap_before\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
        }
        SUBCASE("calendar_dt") {
            SUBCASE("exact") {
                // data
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0 = utc_ptr->time(2016, 1, 1);
                // -----
                ta::calendar_dt c_ta{ utc_ptr, t0, dt, n };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(c_ta), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(c_ta), 10. };
                // -----
                std::string fn("dtss_save_merge/calendar_exact.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/calendar_exact\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/calendar_exact.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::CALENDAR);
                FAST_CHECK_EQ(res.total_period().start, c_ta.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, c_ta.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 10.);
                FAST_CHECK_EQ(res.v.at(n - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/calendar_exact\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new contained in old") {
                // data
                const std::size_t drop = 10u;  // points to drop from start/end of old
                // -----
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old + dt*drop;
                // -----
                ta::calendar_dt c_ta_old{ utc_ptr, t0_old, dt, n };
                ta::calendar_dt c_ta_new{ utc_ptr, t0_new, dt, n - 2*drop };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(c_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(c_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/calendar_new_in_old.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/calendar_new_in_old\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/calendar_new_in_old.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::CALENDAR);
                FAST_CHECK_EQ(res.total_period().start, c_ta_old.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, c_ta_old.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 1.);
                FAST_CHECK_EQ(res.v.at(drop - 1), 1.);
                FAST_CHECK_EQ(res.v.at(drop), 10.);
                FAST_CHECK_EQ(res.v.at(n - drop - 1), 10.);
                FAST_CHECK_EQ(res.v.at(n - drop), 1.);
                FAST_CHECK_EQ(res.v.at(n - 1), 1.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/calendar_new_in_old\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("old contained in new") {
                // data
                const std::size_t extra = 10u;  // points to drop from start/end of old
                // -----
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old - dt*extra;
                // -----
                ta::calendar_dt c_ta_old{ utc_ptr, t0_old, dt, n };
                ta::calendar_dt c_ta_new{ utc_ptr, t0_new, dt, n + extra };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(c_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(c_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/calendar_old_in_new.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/calendar_old_in_new\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/calendar_old_in_new.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::CALENDAR);
                FAST_CHECK_EQ(res.total_period().start, c_ta_new.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, c_ta_new.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 10.);
                FAST_CHECK_EQ(res.v.at(n - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/calendar_old_in_new\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new overlap start of old") {
                // data
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old - dt*n/2;
                // -----
                ta::calendar_dt c_ta_old{ utc_ptr, t0_old, dt, n };
                ta::calendar_dt c_ta_new{ utc_ptr, t0_new, dt, n };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(c_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(c_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/calendar_new_over_start.db");

                FAST_REQUIRE_LT(t0_new, t0_old);

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/calendar_new_over_start\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/calendar_new_over_start.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::CALENDAR);
                FAST_CHECK_EQ(res.total_period().start, c_ta_new.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, c_ta_old.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 10.);
                FAST_CHECK_EQ(res.v.at(n - 1), 10.);
                FAST_CHECK_EQ(res.v.at(n), 1.);
                FAST_CHECK_EQ(res.v.at(n + n/2 - 1), 1.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/calendar_new_over_start\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new overlap end of old") {
                // data
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old + dt*n/2;
                // -----
                ta::calendar_dt c_ta_old{ utc_ptr, t0_old, dt, n };
                ta::calendar_dt c_ta_new{ utc_ptr, t0_new, dt, n };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(c_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(c_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/calendar_new_over_end.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/calendar_new_over_end\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/calendar_new_over_end.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::CALENDAR);
                FAST_CHECK_EQ(res.total_period().start, c_ta_old.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, c_ta_new.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 1.);
                FAST_CHECK_EQ(res.v.at(n/2 - 1), 1.);
                FAST_CHECK_EQ(res.v.at(n/2), 10.);
                FAST_CHECK_EQ(res.v.at(n + n/2 - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/calendar_new_over_end\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("consecutive without gap") {
                // data
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old + dt*n;
                // -----
                ta::calendar_dt c_ta_old{ utc_ptr, t0_old, dt, n };
                ta::calendar_dt c_ta_new{ utc_ptr, t0_new, dt, n };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(c_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(c_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/calendar_consec.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/calendar_consec\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/calendar_consec.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::CALENDAR);
                FAST_CHECK_EQ(res.total_period().start, c_ta_old.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, c_ta_new.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 1.);
                FAST_CHECK_EQ(res.v.at(n - 1), 1.);
                FAST_CHECK_EQ(res.v.at(n), 10.);
                FAST_CHECK_EQ(res.v.at(2 * n - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/calendar_consec\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new after old with gap") {
                // data
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old + 2*n*dt;
                // -----
                ta::calendar_dt c_ta_old{ utc_ptr, t0_old, dt, n };
                ta::calendar_dt c_ta_new{ utc_ptr, t0_new, dt, n };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(c_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(c_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/calendar_gap_after.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/calendar_gap_after\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/calendar_gap_after.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::CALENDAR);
                FAST_CHECK_EQ(res.ta.size(), 3*n);
                FAST_CHECK_EQ(res.total_period().start, c_ta_old.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, c_ta_new.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 1.);
                FAST_CHECK_EQ(res.v.at(n - 1), 1.);
                FAST_CHECK_UNARY(std::isnan(res.v.at(n)));
                FAST_CHECK_UNARY(std::isnan(res.v.at(2*n - 1)));
                FAST_CHECK_EQ(res.v.at(2 * n), 10.);
                FAST_CHECK_EQ(res.v.at(3 * n - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/calendar_gap_after\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new before old with gap") {
                // data
                const core::utctimespan dt = core::calendar::DAY;
                const std::size_t n = 100;
                const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
                const core::utctime t0_new = t0_old - 2*n*dt;
                // -----
                ta::calendar_dt c_ta_old{ utc_ptr, t0_old, dt, n };
                ta::calendar_dt c_ta_new{ utc_ptr, t0_new, dt, n };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(c_ta_old), 1. };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(c_ta_new), 10. };
                // -----
                std::string fn("dtss_save_merge/calendar_gap_before.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/calendar_gap_before\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/calendar_gap_before.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::CALENDAR);
                FAST_CHECK_EQ(res.ta.size(), 3*n);
                FAST_CHECK_EQ(res.total_period().start, c_ta_new.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, c_ta_old.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 10.);
                FAST_CHECK_EQ(res.v.at(n - 1), 10.);
                FAST_CHECK_UNARY(std::isnan(res.v.at(n)));
                FAST_CHECK_UNARY(std::isnan(res.v.at(2*n - 1)));
                FAST_CHECK_EQ(res.v.at(2 * n), 1.);
                FAST_CHECK_EQ(res.v.at(3 * n - 1), 1.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/calendar_gap_before\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
        }
        SUBCASE("point_dt") {
            // data
            const core::utctimespan dt = core::calendar::DAY;
            const std::size_t n = 100;
            const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
            const core::utctimespan offset = core::deltaminutes(25);  // to misalign old and new timepoints
            // -----
            std::vector<core::utctime> old_timepoints{};
            const std::vector<double> old_values(n, 1.);
            // -----
            std::vector<core::utctime> new_timepoints{};
            std::vector<double> new_values(n, 10.);
            // -----
            old_timepoints.reserve(n + 1);
            for ( std::size_t i = 0; i <= n; ++i ) {
                old_timepoints.emplace_back(t0_old + dt*i);
            }

            SUBCASE("exact") {
                // data
                ta::point_dt p_ta{ old_timepoints };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(p_ta), old_values };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(p_ta), new_values };
                // -----
                std::string fn("dtss_save_merge/point_exact.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/point_exact\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/point_exact.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::POINT);
                FAST_CHECK_EQ(res.total_period().start, p_ta.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, p_ta.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 10.);
                FAST_CHECK_EQ(res.v.at(n - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/point_exact\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new contained in old") {
                // data
                const std::size_t drop = 10u;  // points to drop from start/end of old
                const core::utctime t0_new = t0_old + offset + dt*drop;
                // -----
                new_timepoints.reserve(n - 2 * drop);
                for ( std::size_t i = 0; i <= n - 2 * drop; ++i ) {
                    new_timepoints.emplace_back(t0_new + i * dt);
                }
                new_values.resize(new_timepoints.size() - 1);
                // -----
                ta::point_dt p_ta_old{ old_timepoints };
                ta::point_dt p_ta_new{ new_timepoints };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(p_ta_old), old_values };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(p_ta_new), new_values };
                // -----
                std::string fn("dtss_save_merge/point_new_in_old.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/point_new_in_old\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/point_new_in_old.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::POINT);
                FAST_CHECK_EQ(res.total_period().start, p_ta_old.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, p_ta_old.total_period().end);
                // values
                //  - NB! (n+1) because the point_dt's are unaligned, so a new point is introduced
                FAST_CHECK_EQ(res.v.at(0), 1.);
                FAST_CHECK_EQ(res.v.at(drop), 1.);  // this is the extra new point
                FAST_CHECK_EQ(res.v.at(drop + 1), 10.);
                FAST_CHECK_EQ(res.v.at((n+1) - drop - 1), 10.);
                FAST_CHECK_EQ(res.v.at((n+1) - drop), 1.);
                FAST_CHECK_EQ(res.v.at((n+1) - 1), 1.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/point_new_in_old\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("old contained in new") {
                // data
                const std::size_t extra = 10u;  // points to drop from start/end of old
                const core::utctime t0_new = t0_old + offset - dt*extra;
                // -----
                new_timepoints.reserve(n + 2 * extra);
                for ( std::size_t i = 0; i <= n + 2 * extra; ++i ) {
                    new_timepoints.emplace_back(t0_new + i * dt);
                }
                new_values.resize(new_timepoints.size() - 1, 10.);
                // -----
                ta::point_dt p_ta_old{ old_timepoints };
                ta::point_dt p_ta_new{ new_timepoints };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(p_ta_old), old_values };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(p_ta_new), new_values };
                // -----
                std::string fn("dtss_save_merge/point_old_in_new.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/point_old_in_new\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/point_old_in_new.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::POINT);
                FAST_CHECK_EQ(res.total_period().start, p_ta_new.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, p_ta_new.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 10.);
                FAST_CHECK_EQ(res.v.at(n - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/point_old_in_new\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new overlap start of old") {
                // data
                const core::utctime t0_new = t0_old + offset - dt*n/2;
                for ( std::size_t i = 0; i <= n; ++i ) {
                    new_timepoints.emplace_back(t0_new + i * dt);
                }
                // -----
                ta::point_dt p_ta_old{ old_timepoints };
                ta::point_dt p_ta_new{ new_timepoints };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(p_ta_old), old_values };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(p_ta_new), new_values };
                // -----
                std::string fn("dtss_save_merge/point_new_over_start.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/point_new_over_start\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/point_new_over_start.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::POINT);
                FAST_CHECK_EQ(res.total_period().start, p_ta_new.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, p_ta_old.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 10.);
                FAST_CHECK_EQ(res.v.at(n - 1), 10.);
                FAST_CHECK_EQ(res.v.at(n), 1.);
                FAST_CHECK_EQ(res.v.at(n + n/2 - 1), 1.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/point_new_over_start\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new overlap end of old") {
                // data
                const core::utctime t0_new = t0_old + offset + dt*n/2;
                for ( std::size_t i = 0; i <= n; ++i ) {
                    new_timepoints.emplace_back(t0_new + i * dt);
                }
                // -----
                ta::point_dt p_ta_old{ old_timepoints };
                ta::point_dt p_ta_new{ new_timepoints };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(p_ta_old), old_values };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(p_ta_new), new_values };
                // -----
                std::string fn("dtss_save_merge/point_new_over_end.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/point_new_over_end\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/point_new_over_end.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::POINT);
                FAST_CHECK_EQ(res.total_period().start, p_ta_old.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, p_ta_new.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 1.);
                FAST_CHECK_EQ(res.v.at(n/2), 1.);
                FAST_CHECK_EQ(res.v.at(n/2 + 1), 10.);
                FAST_CHECK_EQ(res.v.at(n + n/2), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/point_new_over_end\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("consecutive without gap") {
                // data
                const core::utctime t0_new = t0_old + n * dt;
                for ( std::size_t i = 0; i <= n; ++i ) {
                    new_timepoints.emplace_back(t0_new + i * dt);
                }
                // -----
                ta::point_dt p_ta_old{ old_timepoints };
                ta::point_dt p_ta_new{ new_timepoints };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(p_ta_old), old_values };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(p_ta_new), new_values };
                // -----
                std::string fn("dtss_save_merge/point_consec.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/point_consec\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/point_consec.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::POINT);
                FAST_CHECK_EQ(res.total_period().start, p_ta_old.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, p_ta_new.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 1.);
                FAST_CHECK_EQ(res.v.at(n - 1), 1.);
                FAST_CHECK_EQ(res.v.at(n), 10.);
                FAST_CHECK_EQ(res.v.at(2 * n - 1), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/point_consec\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new after old with gap") {
                // data
                const core::utctime t0_new = t0_old + n * dt + offset;
                for ( std::size_t i = 0; i <= n; ++i ) {
                    new_timepoints.emplace_back(t0_new + i * dt);
                }
                // -----
                ta::point_dt p_ta_old{ old_timepoints };
                ta::point_dt p_ta_new{ new_timepoints };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(p_ta_old), old_values };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(p_ta_new), new_values };
                // -----
                std::string fn("dtss_save_merge/point_gap_after.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/point_gap_after\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/point_gap_after.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::POINT);
                FAST_CHECK_EQ(res.ta.size(), 2 * n + 1);
                FAST_CHECK_EQ(res.total_period().start, p_ta_old.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, p_ta_new.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 1.);
                FAST_CHECK_EQ(res.v.at(n - 1), 1.);
                FAST_CHECK_UNARY(std::isnan(res.v.at(n)));
                FAST_CHECK_EQ(res.v.at(n + 1), 10.);
                FAST_CHECK_EQ(res.v.at(2 * n), 10.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/point_gap_after\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
            SUBCASE("new before old with gap") {
                // data
                const core::utctime t0_new = t0_old - n * dt - offset;
                for ( std::size_t i = 0; i <= n; ++i ) {
                    new_timepoints.emplace_back(t0_new + i * dt);
                }
                // -----
                ta::point_dt p_ta_old{ old_timepoints };
                ta::point_dt p_ta_new{ new_timepoints };
                ts::point_ts<ta::generic_dt> pts_old{ gta_t(p_ta_old), old_values };
                ts::point_ts<ta::generic_dt> pts_new{ gta_t(p_ta_new), new_values };
                // -----
                std::string fn("dtss_save_merge/point_gap_before.db");

                // save initital data
                db.save(fn, pts_old, false);
                auto find_res = db.find("dtss_save_merge/point_gap_before\\.db");
                FAST_CHECK_EQ(find_res.size(), 1);

                // add data to the same path
                db.save(fn, pts_new, false);

                // check merged data
                ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/point_gap_before.db", core::utcperiod{ });
                // time-axis
                FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::POINT);
                FAST_CHECK_EQ(res.ta.size(), 2 * n + 1);
                FAST_CHECK_EQ(res.total_period().start, p_ta_new.total_period().start);
                FAST_CHECK_EQ(res.total_period().end, p_ta_old.total_period().end);
                // values
                FAST_CHECK_EQ(res.v.at(0), 10.);
                FAST_CHECK_EQ(res.v.at(n - 1), 10.);
                FAST_CHECK_UNARY(std::isnan(res.v.at(n)));
                FAST_CHECK_EQ(res.v.at(n + 1), 1.);
                FAST_CHECK_EQ(res.v.at(2 * n), 1.);

                // cleanup
                db.remove(fn);
                find_res = db.find(string("dtss_save_merge/point_gap_before\\.db"));
                FAST_CHECK_EQ(find_res.size(), 0);
            }
        }
        SUBCASE("force overwrite") {
            // data
            const std::size_t drop = 10u;  // points to drop from start/end of old
            // -----
            const core::utctimespan dt = core::calendar::DAY;
            const std::size_t n = 100;
            const core::utctime t0_old = utc_ptr->time(2016, 1, 1);
            const core::utctime t0_new = t0_old + dt*drop;
            // -----
            ta::fixed_dt f_ta_old{ t0_old, dt, n };
            ta::fixed_dt f_ta_new{ t0_new, dt, n - 2*drop };
            ts::point_ts<ta::generic_dt> pts_old{ gta_t(f_ta_old), 1. };
            ts::point_ts<ta::generic_dt> pts_new{ gta_t(f_ta_new), 10. };
            // -----
            std::string fn("dtss_save_merge/force_overwrite.db");

            // save initital data
            db.save(fn, pts_old, false);
            auto find_res = db.find("dtss_save_merge/force_overwrite\\.db");
            FAST_CHECK_EQ(find_res.size(), 1);

            // add data to the same path with overwrite
            db.save(fn, pts_new, true);

            // check merged data
            ts::point_ts<ta::generic_dt> res = db.read("dtss_save_merge/force_overwrite.db", core::utcperiod{ });
            // time-axis
            FAST_CHECK_EQ(res.ta.gt, time_axis::generic_dt::FIXED);
            FAST_CHECK_EQ(res.total_period().start, f_ta_new.total_period().start);
            FAST_CHECK_EQ(res.total_period().end, f_ta_new.total_period().end);
            // values
            FAST_CHECK_EQ(res.v.at(0), 10.);
            FAST_CHECK_EQ(res.v.at(n - 2 * drop - 1), 10.);

            // cleanup
            db.remove(fn);
            find_res = db.find(string("dtss_save_merge/force_overwrite\\.db"));
            FAST_CHECK_EQ(find_res.size(), 0);
        }
    }
}

TEST_SUITE("dtss_db") {
    TEST_CASE("dtss_db_seconds_time_io") {
        test_time_io<seconds_time_io>();
    }
    TEST_CASE("dtss_db_native_time_io") {
        test_time_io<native_time_io>();
    }
    TEST_CASE("dtss_db_basics_micro_seconds") {
        test_dtss_db_basics(true);
    }
    TEST_CASE("dtss_db_basics_seconds") {
        test_dtss_db_basics(false);
    }
   TEST_CASE("dtss_db_merge_write_micro_seconds") {
        test_dtss_db_store_merge_write(true);
    }
    TEST_CASE("dtss_db_merge_write_seconds") {
        test_dtss_db_store_merge_write(false);
    }
    TEST_CASE("dtss_db_merge_us_to_seconds_throws") {
        using namespace shyft::dtss;
        using namespace shyft::time_series::dd;
        using time_series::point_ts;

        std::shared_ptr<core::calendar> utc = std::make_shared<core::calendar>();
        std::shared_ptr<core::calendar> osl = std::make_shared<core::calendar>("Europe/Oslo");

        core::utctime t = utc->time(2016, 1, 1);
        core::utctimespan dt = core::deltahours(1);
        std::size_t n = 24 * 365 * 2;//24*365*5;

        // construct time-axis that we want to test.
        time_axis::fixed_dt fta(t, dt, n);
        constexpr bool time_in_micro_seconds=false;// use seconds when writing startts.
        vector<utctime> tp;for(std::size_t i=0;i<fta.size();++i)tp.push_back(fta.time(i));
        time_axis::point_dt pta(tp,fta.total_period().end);
        auto dirname = "ts.db.test." +(time_in_micro_seconds?string(".s."):string(".us."))+ std::to_string(to_seconds64(core::utctime_now()) );
        auto tmpdir = (fs::temp_directory_path()/dirname);
        ts_db db(tmpdir.string());
        db.time_format_micro_seconds=time_in_micro_seconds;
        // 1. store a break-point ts, with seconds resolution
        gts_t o(gta_t(pta),10.0,time_series::ts_point_fx::POINT_AVERAGE_VALUE);
        o.set(0, 100.); o.set(o.size()-1, 100.);
        std::string fn("measurements/tssf.db");
        db.save(fn,o);

        // 2. try to merge another ts, with seconds resolution (should work)
        vector<utctime> tp2;tp2.push_back(tp.back());
        time_axis::point_dt pta2(tp2,tp2.back()+dt);
        gts_t x(gta_t(pta2),9.0,time_series::ts_point_fx::POINT_AVERAGE_VALUE);
        db.save(fn,x,false);// merge store
        // 3. try to merge another ts, with us resolution (should throw).
        tp2[0] += utctime{5000};// add some microseconds
        gts_t us(gta_t(tp2,tp2.back()+dt),7.0,time_series::ts_point_fx::POINT_AVERAGE_VALUE);
        try {
            db.save(fn,us,false);
            FAST_CHECK_UNARY(false);
        } catch( const runtime_error&) {
            FAST_CHECK_UNARY(true);
        }
        db.remove(fn);
    }
}
