#include "test_pch.h"

#include "core/dtss.h"
#include "core/dtss_cache.h"
#include "core/dtss_client.h"

#include "core/utctime_utilities.h"
#include "core/time_axis.h"
#include "core/time_series_merge.h"
#include "core/time_series_dd.h"


#include <future>
#include <mutex>
#include <regex>
#include <boost/filesystem.hpp>
#include <cstdint>
#ifdef _WIN32
#include <io.h>
#else
#include <sys/io.h>
#define O_BINARY 0
#define O_SEQUENTIAL 0
#include <sys/stat.h>
#endif
#include <fcntl.h>
#include <string_view>

namespace  fs=boost::filesystem;
#include <armadillo>

using namespace std;
using namespace shyft;
using namespace shyft::core;
using shyft::time_series::dd::apoint_ts;
using shyft::time_series::dd::gta_t;
static inline shyft::core::utctime _t(int64_t t1970s) {return shyft::core::utctime{shyft::core::seconds(t1970s)};}

apoint_ts mk_expression(utctime t, utctimespan dt, int n) {

    std::vector<double> x; x.reserve(n);
    for (int i = 0; i < n; ++i)
        x.push_back(-double(n) / 2.0 + i);
    apoint_ts aa(gta_t(t, dt, n), x);
    auto a = aa*3.0 + aa;
    return a;
}

dlib::logger dlog("dlib.log");

#define TEST_SECTION(x)

TEST_SUITE("dtss") {

TEST_CASE("dtss_lru_cache") {
    using shyft::dtss::lru_cache;
    using std::map;
    using std::list;
    using std::vector;
    using std::string;
    using std::back_inserter;
    using shyft::time_series::dd::apoint_ts;
    using shyft::time_series::dd::gta_t;
    const auto stair_case=shyft::time_series::POINT_AVERAGE_VALUE;
    lru_cache<string, apoint_ts, map > c(2);

    apoint_ts r;
    vector<string> mru;
    gta_t ta(_t(0), seconds(1), 10);

    TEST_SECTION("empty_cache") {
        FAST_CHECK_UNARY_FALSE(c.try_get_item("a", r));
    }
    TEST_SECTION("add_one_item") {
        c.add_item("a", apoint_ts(ta, 1.0, stair_case));
        FAST_CHECK_UNARY(c.try_get_item("a", r));
        FAST_CHECK_EQ(ta.size(), r.time_axis().size());
        FAST_CHECK_UNARY_FALSE(c.try_get_item("b", r));
    }
    TEST_SECTION("add_second_item") {
        c.add_item("b", apoint_ts(ta, 2.0, stair_case));
        FAST_CHECK_UNARY(c.try_get_item("a", r));
        FAST_CHECK_UNARY(c.try_get_item("b", r));
        c.get_mru_keys(back_inserter(mru));
        FAST_CHECK_EQ(string("b"), mru[0]);
        FAST_CHECK_EQ(string("a"), mru[1]);
    }
    TEST_SECTION("mru_item_in_front") {
        c.try_get_item("a", r);
        mru.clear(); c.get_mru_keys(back_inserter(mru));
        FAST_CHECK_EQ(string("a"), mru[0]);
        FAST_CHECK_EQ(string("b"), mru[1]);
    }
    TEST_SECTION("excessive_lru_item_evicted_when_adding") {
        c.add_item("c", apoint_ts(ta, 3.0, stair_case));
        FAST_CHECK_UNARY_FALSE(c.try_get_item("b", r));
        FAST_CHECK_UNARY(c.try_get_item("c", r));
        FAST_CHECK_UNARY(c.try_get_item("a", r));
    }
    TEST_SECTION("remove_item") {
        c.remove_item("a");
        FAST_CHECK_UNARY_FALSE(c.try_get_item("a", r));
    }
    TEST_SECTION("ensure_items_added_are_first") {
        c.add_item("d", apoint_ts(ta, 4.0, stair_case));
        mru.clear(); c.get_mru_keys(back_inserter(mru));
        FAST_CHECK_EQ(string("d"), mru[0]);
        FAST_CHECK_EQ(string("c"), mru[1]);
    }
    TEST_SECTION("update_existing") {
        c.try_get_item("c", r);//just to ensure "c" is in first position
        c.add_item("d", apoint_ts(ta, 4.2, stair_case)); //update "d"
        c.try_get_item("d", r);
        FAST_CHECK_GT(r.value(0), 4.1);
        mru.clear(); c.get_mru_keys(back_inserter(mru));
        FAST_CHECK_EQ(string("d"), mru[0]);
        FAST_CHECK_EQ(string("c"), mru[1]);
    }
}
TEST_CASE("dtss_ts_cache") {
    using std::vector;
    using std::string;
    using shyft::core::utctime;
    using shyft::core::deltahours;
    using shyft::dtss::cache_stats;
    using shyft::time_series::dd::apoint_ts;
    using shyft::time_series::dd::gta_t;
    using shyft::dtss::apoint_ts_frag;
    using dtss_cache=shyft::dtss::cache<apoint_ts_frag,apoint_ts>;
    const auto stair_case=shyft::time_series::POINT_AVERAGE_VALUE;
    size_t max_ids=10;

    dtss_cache c(max_ids);

    apoint_ts x;
    utcperiod p{_t(0),_t(10)};
    FAST_CHECK_EQ(false,c.try_get("a",p,x));

    utctime t0{_t(5)};
    utctime t1{_t(10)};

    utctimespan dt{seconds(1)};
    size_t n{3};
    apoint_ts ts_a{gta_t{t0,dt,n},1.0,stair_case};
    apoint_ts ts_a2{gta_t{t1,dt,n},1.0,stair_case};

    c.add("a",ts_a);
    FAST_REQUIRE_EQ(true,c.try_get("a",utcperiod{t0,t0+dt},x));
    FAST_CHECK_EQ(1.0,x.value(0));

    c.add("a",ts_a2);// test we can add twice, a fragment (replace)

    FAST_REQUIRE_EQ(false,c.try_get("a",utcperiod{t0,t0+ 10*dt},x));

    auto s=c.get_cache_stats();
    FAST_CHECK_EQ(s.hits,2);
    FAST_CHECK_EQ(s.misses,1);
    FAST_CHECK_EQ(s.coverage_misses,1);
    FAST_CHECK_EQ(s.point_count,6);
    FAST_CHECK_EQ(s.fragment_count,2);
    FAST_CHECK_EQ(s.id_count,1);

    c.clear_cache_stats();
    c.flush();
    s = c.get_cache_stats();
    FAST_CHECK_EQ(s.hits, 0);
    FAST_CHECK_EQ(s.misses, 0);
    FAST_CHECK_EQ(s.coverage_misses, 0);
    FAST_CHECK_EQ(s.point_count, 0);
    FAST_CHECK_EQ(s.fragment_count, 0);
    FAST_CHECK_EQ(s.id_count, 0);

    c.add("a", ts_a2);

    c.remove("a");
    FAST_REQUIRE_EQ(false,c.try_get("a",utcperiod{t0,t0+dt},x));

    //-- test vector operations
    // arrange n_ts
    vector<string> ids;
    vector<apoint_ts> tss;
    size_t n_ts = 3;
    gta_t mta{ t0,dt,n };
    for (size_t i = 0; i<n_ts; ++i) {
        ids.push_back(to_string(i));
        tss.emplace_back(mta, double(i), stair_case);
    }

    c.add(ids, tss); // add a vector of ids|tss

    auto mts = c.get(ids, mta.total_period());// get a vector of  ids back as map[id]=ts
    FAST_REQUIRE_EQ(n_ts, mts.size());
    for (size_t i = 0; i<n_ts; ++i) {
        FAST_REQUIRE_UNARY(mts.find(ids[i])!=mts.end());
        FAST_CHECK_EQ(mts[ids[i]].value(0), double(i)); // just check one value unique for ts.
    }

    auto ids2 = ids; ids2.push_back("not there");// ask for something that's not there
    auto mts2 = c.get(ids2, mta.total_period());
    FAST_REQUIRE_EQ(n_ts, mts.size());
    for (size_t i = 0; i<n_ts; ++i) {
        FAST_REQUIRE_UNARY(mts2.find(ids[i]) != mts2.end());
        FAST_CHECK_EQ(mts2[ids[i]].value(0), double(i)); // just check one value unique for ts.
    }

    c.remove(ids2); // remove by vector (even with elem not there)
    s = c.get_cache_stats();
    FAST_CHECK_EQ(s.point_count, 0);
    FAST_CHECK_EQ(s.fragment_count, 0);
    FAST_CHECK_EQ(s.id_count, 0);

}
TEST_CASE("dtss_mini_frag") {
    using std::vector;
    using std::string;
    using std::min;
    using std::max;
    using shyft::core::utcperiod;
    using shyft::core::utctime;
    using shyft::dtss::mini_frag;
    using shyft::time_axis::continuous_merge;

    struct tst_frag {
        utcperiod p;
        int id{0};
        tst_frag(utctime f,utctime u):p(f,u){}
        tst_frag(utctime f,utctime u,int id):p(f,u),id(id){}
        utcperiod total_period() const {return p;}
        size_t size() const {return size_t(to_seconds64(p.timespan()) );}
        tst_frag merge(const tst_frag&o) const {
            if(!continuous_merge(p,o.p))
                throw runtime_error("Wrong merge op attempted");
            return tst_frag{min(o.p.start,p.start),max(o.p.end,p.end)};
        }
    };

    mini_frag<tst_frag> m;

    FAST_CHECK_EQ(m.count_fragments(),0);
    FAST_CHECK_EQ(m.estimate_size(),0);
    FAST_CHECK_EQ(m.get_ix(utcperiod(_t(0),_t(10))),string::npos);

    m.add(tst_frag(_t(5),_t(10)) );// add first
    FAST_CHECK_EQ(m.count_fragments(),1);
    FAST_CHECK_EQ(m.estimate_size(),5);
    FAST_CHECK_EQ(m.get_ix(utcperiod{_t(5),_t(10)}),0);

    m.add(tst_frag(_t(4),_t(5) ));// add merge element
    FAST_CHECK_EQ(m.count_fragments(),1);
    FAST_CHECK_EQ(m.estimate_size(),6);
    FAST_CHECK_EQ(m.get_ix(utcperiod{_t(4),_t(6)}),0);

    m.add(tst_frag(_t(1),_t(3) ));// add before first
    FAST_CHECK_EQ(m.count_fragments(),2);
    FAST_CHECK_EQ(m.estimate_size(),8);
    FAST_CHECK_EQ(m.get_ix(utcperiod{_t(4),_t(6)}),1);
    FAST_CHECK_EQ(m.get_ix(utcperiod{_t(1),_t(2)}),0);

    m.add(tst_frag(_t(3),_t(4) ));// add a piece that merge [0] [1]
    FAST_CHECK_EQ(m.count_fragments(),1);
    FAST_CHECK_EQ(m.estimate_size(),9);
    FAST_CHECK_EQ(m.get_ix(utcperiod{_t(4),_t(6)}),0);

    m.add(tst_frag{_t(11),_t(12) }); // append a frag at the end
    FAST_CHECK_EQ(m.count_fragments(),2);
    FAST_CHECK_EQ(m.estimate_size(),10);
    FAST_CHECK_EQ(m.get_ix(utcperiod{_t(4),_t(6)}),0);
    FAST_CHECK_EQ(m.get_ix(utcperiod{_t(11),_t(12)}),1);
    FAST_CHECK_EQ(m.get_ix(utcperiod{_t(11),_t(13)}),string::npos);
    FAST_CHECK_EQ(m.get_ix(utcperiod{_t(1),_t(12)}),string::npos);

    m.add(tst_frag{_t(13),_t(15) }); // append a frag at the end
    FAST_CHECK_EQ(m.count_fragments(),3);
    FAST_CHECK_EQ(m.estimate_size(),12);
    FAST_CHECK_EQ(m.get_ix(utcperiod{_t(13),_t(14)}),2);
    FAST_CHECK_EQ(m.get_ix(utcperiod{_t(11),_t(12)}),1);
    FAST_CHECK_EQ(m.get_ix(utcperiod{_t(5),_t(6)}),0);

    m.add(tst_frag{_t(11), _t(14)});// append a frag that melts [1]..[2] into one
    FAST_CHECK_EQ(m.count_fragments(),2);
    FAST_CHECK_EQ(m.get_ix(utcperiod{11,15}),1);

    m.add(tst_frag{_t(0), _t(20)}); // append a frag that melts all into one
    FAST_CHECK_EQ(m.count_fragments(),1);
    FAST_CHECK_EQ(m.estimate_size(),20);

    m.add(tst_frag{_t(0),_t(20),1}); // append a marked frag that is exactly equal
    FAST_CHECK_EQ(m.count_fragments(),1);
    FAST_CHECK_EQ(m.estimate_size(),20);
    FAST_CHECK_EQ(m.get_ix(utcperiod{0,10}),0);
    FAST_CHECK_EQ(m.get_by_ix(0).id,1);// verify we got the marked in place

    m.add(tst_frag{_t(21),_t(23)});// two frags
    m.add(tst_frag{_t(25),_t(27)});// three
    m.add(tst_frag{_t(0),_t(23)});// exactly cover second frag
    FAST_CHECK_EQ(m.count_fragments(),2);

    m.add(tst_frag{_t(1),_t(24)});
    FAST_CHECK_EQ(m.count_fragments(),2);

    m.add(tst_frag{_t(2),_t(27)});// parts p1, p2.end
    FAST_CHECK_EQ(m.count_fragments(),1);

    m.add(tst_frag{_t(-1),_t(27)});
    FAST_CHECK_EQ(m.count_fragments(),1);

}

TEST_CASE("dlib_server_basics") {
    dlog.set_level(dlib::LALL);
    dlib::set_all_logging_output_streams(std::cout);
    dlib::set_all_logging_levels(dlib::LALL);
    dlog << dlib::LINFO << "Starting dtss test";
    using namespace shyft::dtss;
    try {
        // Tell the server to begin accepting connections.
        calendar utc;
        auto t = utc.time(2016, 1, 1);
        auto dt = deltahours(1);
        auto dt24 = deltahours(24);
        int n = 240;
        int n24 = n / 24;
        time_axis::fixed_dt ta(t, dt, n);
        gta_t ta24(t, dt24, n24);
        bool throw_exception = false;
        read_call_back_t cb = [ta, &throw_exception](id_vector_t ts_ids, core::utcperiod p)
            ->ts_vector_t {
            ts_vector_t r; r.reserve(ts_ids.size());
            double fv = 1.0;
            for (size_t i = 0; i < ts_ids.size(); ++i)
                r.emplace_back(ta, fv += 1.0);
            if (throw_exception) {
                dlog << dlib::LINFO << "Throw from inside dtss executes!";
                throw std::runtime_error("test exception");
            }
            return r;
        };
        std::vector<std::string> ts_names = {
            string("a.prod.mw"),
            string("b.prod.mw")
        };
        find_call_back_t fcb = [&ts_names](std::string search_expression)
            ->ts_info_vector_t {
            ts_info_vector_t r;
            dlog << dlib::LINFO << "find-callback with search-string:" << search_expression;
            std::regex re(search_expression);
            auto match_end = std::sregex_iterator();
            for (auto const&tsn : ts_names) {
                if (std::sregex_iterator(tsn.begin(), tsn.end(), re)!= match_end) {
                    ts_info tsi; tsi.name = tsn;
                    r.push_back(tsi);
                }
            }
            return r;
        };

        server<standard_dtss_dispatcher> our_server(cb,fcb);

        // set up the server object we have made
        our_server.set_listening_ip("127.0.0.1");
        //int port_no = 20000;
        //our_server.set_listening_port(port_no);
        int port_no = our_server.start_server();
        dlog<<dlib::LINFO<<"Serving connection @"<<port_no;
        {
            string host_port = string("localhost:") + to_string(port_no);
            dlog << dlib::LINFO << "sending an expression ts to " << host_port ;
            std::vector<apoint_ts> tsl;
            for (size_t kb = 4;kb < 16;kb += 2)
                tsl.push_back(mk_expression(t, dt, kb * 1000)*apoint_ts(string("netcdf://group/path/ts") + std::to_string(kb)));
            client dtss(host_port);
            auto ts_b = dtss.evaluate(tsl, ta.total_period(),false,false);
            dlog << dlib::LINFO << "Got vector back, size= " << ts_b.size();
            for (const auto& ts : ts_b)
                dlog << dlib::LINFO << "ts.size()" << ts.size();
            dlog << dlib::LINFO << "testing 2 time:";
            FAST_REQUIRE_UNARY(our_server.is_running());
            dtss.evaluate(tsl, ta.period(0),false,false);
            dlog << dlib::LINFO << "done second test";
            // test search functions
            dlog << dlib::LINFO << "test .find function";
            auto found_ts = dtss.find(string("a.*"));
            FAST_REQUIRE_EQ(found_ts.size(), 1);
            FAST_CHECK_EQ(found_ts[0].name, ts_names[0]);
            dlog << dlib::LINFO << "test .find function done";

            throw_exception = true;// verify server-side exception gets back here.
            TS_ASSERT_THROWS_ANYTHING(dtss.evaluate(tsl, ta.period(0),false,false));
            dlog<<dlib::LINFO << "exceptions done,testing ordinary evaluate after exception";
            throw_exception = false;// verify server-side exception gets back here.
            dtss.evaluate(tsl, ta.period(0),false,false); // verify no exception here, should still work ok
            FAST_REQUIRE_UNARY(our_server.is_running());
            dlog << dlib::LINFO << "ok, -now testing percentiles";
            std::vector<int64_t> percentile_spec{ 0,25,50,-1,75,100 };
            auto percentiles = dtss.percentiles(tsl, ta.total_period(), ta24, percentile_spec,false,false);
            FAST_CHECK_EQ(percentiles.size(), percentile_spec.size());
            FAST_CHECK_EQ(percentiles[0].size(), ta24.size());
            dlog << dlib::LINFO << "done with percentiles, stopping localhost server";
            dtss.close();
            our_server.clear();
            dlog << dlib::LINFO << "done";
        }
    } catch (exception& e) {
        cout << e.what() << endl;
    }
    dlog << dlib::LINFO << "done";
}

TEST_CASE("dlib_server_cache_update") {
    //
    // This test is the the very special case that custom-made reader-callbacks returns back
    // time-series with different resolution etc. if the time-series is not found (instead of fail/exception)
    // This was failing in those cases where the next sequent cached write operation ended up in
    // a cache merge, with incompatible time-series.
    // This should however work, given that the user passes the overwrite flag, meaning replace
    // the whole time-series. The cache should then be replaced for those items.
    //
    // this test was constructed so that the initial case reported was failing, 
    // then fix added, and now working ok.
    //
    using namespace shyft::dtss;
    try {
        // Tell the server to begin accepting connections.
        calendar utc;
        auto t = utc.time(2016, 1, 1);
        auto dt = deltahours(1);
        auto dt24 = deltahours(24);
        int n = 240;
        int n24 = n / 24;
        time_axis::fixed_dt ta(t, dt, n);
        gta_t ta24(t, dt24, n24);
        bool return_nan_24 = true;
        map<string,ts_vector_t::value_type> memory_store;

        read_call_back_t cb = [ta,ta24, &return_nan_24](id_vector_t ts_ids, core::utcperiod p)
            ->ts_vector_t {
            ts_vector_t r; r.reserve(ts_ids.size());
            double fv = 1.0;
            auto dt=return_nan_24?deltahours(24):deltahours(1);
            gta_t rta{p.start,dt,static_cast<size_t>(1u+p.timespan()/dt)};
            for (size_t i = 0; i < ts_ids.size(); ++i){
                if(return_nan_24) r.emplace_back(rta,shyft::nan);
                else r.emplace_back(rta, fv += 1.0);
            }
            return r;
        };
        find_call_back_t fcb = [&memory_store](std::string search_expression)
            ->ts_info_vector_t {
            ts_info_vector_t r;
            std::regex re(search_expression);
            auto match_end = std::sregex_iterator();
            for (auto const&kv : memory_store) {
                string tsn{kv.first};
                if (std::sregex_iterator(tsn.begin(), tsn.end(), re)!= match_end) {
                    ts_info tsi; tsi.name = tsn;
                    r.push_back(tsi);
                }
            }
            return r;
        };
        store_call_back_t scb=[&memory_store](const ts_vector_t& tsv)->void {
            for(auto const& ts:tsv) 
                memory_store[ts.id()]=ts;
        };
        
        server<standard_dtss_dispatcher> our_server(cb,fcb,scb);
        our_server.set_listening_ip("127.0.0.1");
        our_server.set_auto_cache(true);// special for this case
        int port_no = our_server.start_server();
        {
            string host_port = string("localhost:") + to_string(port_no);
            std::vector<apoint_ts> tsl;
            std::vector<apoint_ts> tsv_store;
            for (size_t kb = 4;kb < 16;kb += 2) {
                auto ts_url=string("netcdf://group/path/ts") + std::to_string(kb);
                tsl.push_back(mk_expression(t, dt, kb * 1000)*apoint_ts(ts_url));
                tsv_store.push_back(apoint_ts(ts_url,apoint_ts(ta,double(kb),shyft::time_series::POINT_AVERAGE_VALUE)));
            }
            client dtss(host_port);
            auto p_r=ta.period(0);
            p_r.start -= deltahours(24);// subtract 1 day so that we read slighly before, next write will cause merge(not replace)
            auto ts_b = dtss.evaluate(tsl, p_r,true,false);
            dtss.store_ts(tsv_store, true, true); // store these time-series, cache on store
            FAST_CHECK_EQ(memory_store.size(),tsv_store.size());
            return_nan_24=false;// return other values, not nans
            dtss.evaluate(tsl, ta.period(0),true,false);
            FAST_REQUIRE_UNARY(our_server.is_running());
            dtss.close();
            our_server.clear();
        }
    } catch (exception& e) {
        cout << e.what() << endl;
        CHECK_MESSAGE(false,"The server throw exception, test failed");
    }
}


TEST_CASE("dlib_multi_server_basics") {
    dlog.set_level(dlib::LALL);
    dlib::set_all_logging_output_streams(std::cout);
    dlib::set_all_logging_levels(dlib::LALL);
    dlog << dlib::LINFO << "Starting dtss multi-test";
    using namespace shyft::dtss;
    try {
        calendar utc;
        auto t = utc.time(2016, 1, 1);
        auto dt = deltahours(1);
        auto dt24 = deltahours(24);
        int n = 240;
        int n24 = n / 24;
        time_axis::fixed_dt ta(t, dt, n);
        gta_t ta24(t, dt24, n24);
        read_call_back_t rcb = [ta](id_vector_t ts_ids, core::utcperiod p)
            ->ts_vector_t {
            ts_vector_t r; r.reserve(ts_ids.size());
            double fv = 1.0;
            for (size_t i = 0; i < ts_ids.size(); ++i)
                r.emplace_back(ta, fv += 1.0);
            return r;
        };

        size_t n_servers=2;
        vector<unique_ptr<server<standard_dtss_dispatcher>>> servers;
        vector<string> host_ports;
        for(size_t i=0;i<n_servers;++i) {
            auto srv = make_unique<server<standard_dtss_dispatcher>>(rcb);
            srv->set_listening_ip("127.0.0.1");
            //srv->set_listening_port(base_port +i);
            int port_no=srv->start_server();
            servers.emplace_back(move(srv));
            host_ports.push_back(string("localhost:") + to_string(port_no));
        }
        {
            dlog << dlib::LINFO << "sending permutation of an expression sizes/vectors ts to "<<n_servers<< " hosts" ;
            client c(host_ports,false,1000);
            for(size_t n_ts=1;n_ts<10;n_ts+=1) {
                std::vector<apoint_ts> tsl;
                for (size_t kb = 1;kb <= n_ts;kb += 1) {
                    tsl.push_back(mk_expression(t, dt, (kb+1) * 10)*apoint_ts(string("netcdf://group/path/ts") + std::to_string(kb)));
                }
                //dlog << dlib::LINFO<< "Try with size " << tsl.size();
                auto ts_b = c.evaluate(tsl, ta.total_period(),false,false);
                //dlog << dlib::LINFO << "Got vector back, size= " << ts_b.size();
                //for (const auto& ts : ts_b)
                //    dlog << dlib::LINFO << "ts.size()" << ts.size();
                //dlog << dlib::LINFO << "-now testing percentiles";
                std::vector<int64_t> percentile_spec{ 0,25,50,-1,75,100 };
                auto percentiles = c.percentiles(tsl, ta.total_period(), ta24, percentile_spec,false,false);
                //dlog<<dlib::LINFO<<"done calc, verify";
                FAST_CHECK_EQ(percentiles.size(), percentile_spec.size());
                FAST_CHECK_EQ(percentiles[0].size(), ta24.size());
            }
            dlog << dlib::LINFO << "done with multi-server requests, stopping localhost server";
            c.close();
            for(size_t i =0;i<n_servers;++i) {
                dlog<<dlib::LINFO<<"Terminating server "<<i;
                servers[i]->clear();
            }
            dlog << dlib::LINFO << "done";
        }
    } catch (exception& e) {
        cout << e.what() << endl;
    }
    dlog << dlib::LINFO << "done";
}

TEST_CASE("dlib_server_performance") {
    dlog.set_level(dlib::LALL);
    dlib::set_all_logging_output_streams(std::cout);
    dlib::set_all_logging_levels(dlib::LALL);
    dlog << dlib::LINFO << "Starting dtss server performance test";
    using namespace shyft::dtss;
    try {
        // Tell the server to begin accepting connections.
        calendar utc;
        auto t = utc.time(2016, 1, 1);
        auto dt = deltahours(1);
        auto dt24 = deltahours(24);
        int n = 24 * 365 * 5;// 5years of hourly data
        int n24 = n / 24;
        int n_ts = 10;//83;
        gta_t ta(t, dt, n);
        gta_t ta24(t, dt24, n24);
        bool throw_exception = false;
        ts_vector_t from_disk; from_disk.reserve(n_ts);
        double fv = 1.0;
        for (int i = 0; i < n_ts; ++i)
            from_disk.emplace_back(ta, fv += 1.0,shyft::time_series::ts_point_fx::POINT_AVERAGE_VALUE);

        read_call_back_t cb = [&from_disk, &throw_exception](id_vector_t ts_ids, core::utcperiod p)
            ->ts_vector_t {
            if (throw_exception) {
                dlog << dlib::LINFO << "Throw from inside dtss executes!";
                throw std::runtime_error("test exception");
            }
            return from_disk;
        };
        server<standard_dtss_dispatcher> our_server(cb);

        // set up the server object we have made
        our_server.set_listening_ip("127.0.0.1");
        int port_no = our_server.start_server();
        size_t n_threads = 1;

        vector<future<void>> clients;
        for (size_t i = 0;i < n_threads;++i) {
            clients.emplace_back(
                    async(launch::async, [port_no,ta,ta24,i,n_ts]()         /** thread this */ {
                    string host_port = string("localhost:") + to_string(port_no);
                    dlog << dlib::LINFO << "sending an expression ts to " << host_port;
                    std::vector<apoint_ts> tsl;
                    for (int x = 1; x <= n_ts; ++x) {// just make a  very thin request, that get loads of data back
#if 0
                        auto ts_expr = apoint_ts(string("netcdf://group/path/ts_") + std::to_string(x));
                        tsl.push_back(ts_expr);
#else
                        auto ts_expr = 10.0 + 3.0*apoint_ts(string("netcdf://group/path/ts_") + std::to_string(x));
                        if (x > 1) {
                            ts_expr = ts_expr - 3.0*apoint_ts(string("netcdf://group/path/ts_") + std::to_string(x - 1));
                        }
                        tsl.push_back(ts_expr.average(ta));
#endif
                    }

                    client dtss(host_port);
                    auto t0 = timing::now();
                    size_t eval_count = 0;
                    int test_duration_ms = 5000;
                    int kilo_points= tsl.size()*ta.size()/1000;
                    while (elapsed_ms(t0, timing::now()) < test_duration_ms) {
                        // burn cpu server side, save time on serialization
                        std::vector<int64_t> percentile_spec{ -1 };
                        auto percentiles = dtss.percentiles(tsl, ta.total_period(), ta24, percentile_spec,false,false);
                        //slower due to serialization:
                        //auto ts_b = dtss.evaluate(tsl, ta.total_period());
                        ++eval_count;
                    }
                    auto total_ms = double(elapsed_ms(t0, timing::now()));
                    dlog << dlib::LINFO << "Done testing " << i << ": #= " << eval_count
                    << " e/s = " << 1000 * double(eval_count) / total_ms << " [e/s]\n"
                    << " netw throughput = "<< kilo_points*eval_count*8/(total_ms)<< " Mb/sec";
                    dtss.close();
                }

                )
            );
        }
        dlog << dlib::LINFO << "waiting for test to complete";

        for (auto &f : clients) f.get();
        our_server.clear();
        dlog << dlib::LINFO << "done";
    } catch (exception& e) {
        cout << e.what() << endl;
        dlog<<dlib::LERROR<<"exception:"<<e.what();
    }
    dlog << dlib::LINFO << "done";
}

TEST_CASE("dtss_store") { /*
    This test simply create and host a dtss on port 20000,
    then uses shyft:// prefix to test
    all internal operations that involve mapping to the
    shyft ts-db-store.
    */
    using namespace shyft::dtss;
    using namespace shyft::time_series::dd;
    using time_series::point_ts;
    using time_series::ts_point_fx;

    auto utc=make_shared<calendar>();
    auto t = utc->time(2016, 1, 1);
    auto dt = deltahours(1);
    int n = 24 * 365 * 2;//24*365*5;

    // make dtss server
    auto tmpdir = fs::temp_directory_path()/("shyft.c.test"+std::to_string(std::hash<std::thread::id>()(std::this_thread::get_id())));
    server<standard_dtss_dispatcher> our_server{};
    string tc{"tc"};
    our_server.add_container(tc,tmpdir.string());
    our_server.set_listening_ip("127.0.0.1");
    int port_no=our_server.start_server();
    string host_port = string("localhost:") + to_string(port_no);
    // make corresponding client that we will use for the test.
    client dtss(host_port);
    SUBCASE("save_find_read") {
        size_t n_ts=10;
        time_axis::fixed_dt fta(t, dt, n);
        time_axis::generic_dt gta{t,dt*24,size_t(n/24)};
        const auto stair_case=ts_point_fx::POINT_AVERAGE_VALUE;
        ts_vector_t tsv;
        vector<point_ts<time_axis::fixed_dt>> ftsv;

        for(size_t i=0;i<n_ts;++i) {
            tsv.emplace_back(
                    shyft_url(tc,to_string(i)),
                    apoint_ts{fta,i*10.0,stair_case}
            );
            ftsv.emplace_back(fta,i*10.0,stair_case);
        }
        auto f0 = dtss.find(shyft_url(tc,".*"));
        FAST_CHECK_EQ(f0.size(),0);// expect zero to start with
        auto t0=timing::now();
        dtss.store_ts(tsv, true, false);
        auto t1=timing::now();
        auto f1 = dtss.find(shyft_url(tc,".*"));
        FAST_CHECK_EQ(f1.size(),tsv.size());
        ts_vector_t ev;
        for(size_t i=0;i<tsv.size();++i)
            ev.push_back(
                         3.0*apoint_ts(shyft_url(tc,to_string(i)))
                           //+ apoint_ts(shyft_url(tc,to_string(i>0?i-1:i)))
                         );
        // activate auto-cache, to prepare for next
        our_server.set_auto_cache(true);
        vector<int> pc{10,50,90};
        auto t2 = timing::now();
        //auto er= dtss.percentiles(ev,fta.total_period(),gta,pc);//uncached read
        auto er= dtss.evaluate(ev,fta.total_period(),true,true);//uncached read
        auto t3 = timing::now();
        //auto ec = dtss.percentiles(ev, fta.total_period(),gta,pc);
        auto ec = dtss.evaluate(ev, fta.total_period(),true,true);
        auto t4 = timing::now();// cached read.
        //-- establish benchmark
        vector<vector<double>> bmr;bmr.reserve(n_ts);
        for(const auto &ts:ftsv) {
            auto calc= 3.0* ts;
            vector<double> r;r.reserve(calc.size());
            for(size_t i=0;i<calc.size();++i)
                r.emplace_back(calc.value(i));
            bmr.emplace_back(move(r));
        }
        auto t5 =timing::now();
        FAST_CHECK_EQ(bmr.size(),n_ts);
//        FAST_CHECK_EQ(er.size(),pc.size());
//        FAST_CHECK_EQ(ec.size(), pc.size());
        FAST_CHECK_EQ(er.size(),ev.size());
        FAST_CHECK_EQ(ec.size(), ev.size());
        std::cout<<"store mpts/s "<<double(n_ts*n)/(double(elapsed_us(t0,t1))/1000000.0)/1e6<<"\n";
        std::cout<<"evalr mpts/s "<<double(n_ts*n)/(double(elapsed_us(t2,t3))/1000000.0)/1e6<<"\n";
        std::cout<<"evalc mpts/s "<<double(n_ts*n)/(double(elapsed_us(t3,t4))/1000000.0)/1e6<<"\n";
        std::cout<<"bench mpts/s "<<double(n_ts*n)/(double(elapsed_us(t4,t5))/1000000.0)/1e6<<"\t time :"<<double(elapsed_ms(t4,t5))<<"\n";
        auto cs = our_server.get_cache_stats();
        std::cout<<"cache stats(hits,misses,cover_misses,id_count,frag_count,point_count):\n "<<cs.hits<<","<<cs.misses<<","<<cs.coverage_misses<<","<<cs.id_count<<","<<cs.fragment_count<<","<<cs.point_count<<")\n";
    }

    our_server.clear();
#ifdef _WIN32
    for (int i = 0; i<10; ++i) {
        this_thread::sleep_for(chrono::duration<int, std::milli>(1000));
        try {
            fs::remove_all(tmpdir);
            break;
        }
        catch (...) {
            std::cout<<"Try #"<<i+1 << ": Failed to remove " << tmpdir << "\n";
        }
    }
#else
    fs::remove_all(tmpdir);
#endif

}


TEST_CASE("dtss_baseline") {
    using namespace shyft::dtss;
    using namespace shyft::time_series::dd;
    using time_series::point_ts;
    using time_series::ts_point_fx;
    using std::cout;
    auto utc=make_shared<calendar>();
    auto t = utc->time(2016, 1, 1);
    auto dt = deltahours(1);
    const int n = 24 * 365 * 5/3;//24*365*5;

    vector<point_ts<time_axis::fixed_dt>> ftsv;
    const size_t n_ts=10*83;
    arma::mat a_mat(n,n_ts);
    time_axis::fixed_dt fta(t, dt, n);
    //time_axis::generic_dt gta{t,dt*24,size_t(n/24)};
    const auto stair_case=ts_point_fx::POINT_AVERAGE_VALUE;
    ts_vector_t tsv;
    for(size_t i=0;i<n_ts;++i) {
        tsv.emplace_back(to_string(i),apoint_ts(fta,i*10.0,stair_case));
        ftsv.emplace_back(fta,i*10.0,stair_case);
        for(size_t t =0;t<n;++t)
            a_mat(t,i) = i*10.0;
    }
    tsv = 3.0*tsv;


    //-- establish benchmark core-ts
    auto t0 = timing::now();

    vector<vector<double>> bmr;bmr.reserve(n_ts);
    for(const auto &ts:ftsv) {
        auto calc= 3.0* ts;
        vector<double> r;r.reserve(calc.size());
        for(size_t i=0;i<calc.size();++i)
            r.emplace_back(calc.value(i));
        bmr.emplace_back(move(r));
    }
    auto t1 =timing::now();

    //-- establish benchmark armadillo
    vector<vector<double>> amr;amr.reserve(n_ts);
    auto a_res= (a_mat*3.0).eval();
    for(size_t i=0;i<n_ts;++i) {
        amr.emplace_back(arma::conv_to<vector<double>>::from(a_res.col(i)) );
    }
    auto t2 = timing::now();

    //-- establish timing for apoint_ts eval.
    vector<vector<double>> xmr;xmr.reserve(n_ts);
    //auto xtsv = deflate_ts_vector<point_ts<time_axis::generic_dt>>(tsv);
    for(const auto &ts:tsv) {
        xmr.emplace_back(move(ts.values()));
    }
    auto t3 = timing::now();

    FAST_CHECK_EQ(bmr.size(),n_ts);
    FAST_CHECK_EQ(amr.size(),n_ts);
    FAST_CHECK_EQ(xmr.size(),n_ts);

    cout<<"core-ts base-line n_ts= "<<n_ts<<", n="<<n<<", time="<<double(elapsed_us(t0,t1))/1000.0<<"ms ->"
    << double(n*n_ts)/(elapsed_us(t0,t1)/1e6)/1e6<<" mops/s \n";

    cout<<"api -ts base-line n_ts= "<<n_ts<<", n="<<n<<", time="<<double(elapsed_us(t2,t3))/1000.0<<"ms ->"
    << double(n*n_ts)/(elapsed_us(t2,t3)/1e6)/1e6<<" mops/s \n";

    cout<<"armavec base-line n_ts= "<<n_ts<<", n="<<n<<", time="<<double(elapsed_us(t1,t2))/1000.0<<"ms ->"
    << double(n*n_ts)/(elapsed_us(t1,t2)/1e6)/1e6<<" mops/s \n";

}

TEST_CASE("dtss_ltm") {
    // this is basically just for performance study of
    // for api type of ts-expressions,
    using namespace shyft::dtss;
    using namespace shyft::time_series::dd;
    using shyft::time_series::point_ts;
    using time_series::ts_point_fx;
    using std::cout;
    auto utc=make_shared<calendar>();
    auto t = utc->time(2016, 1, 1);
    auto dt = deltahours(1);
    const int n = 24 * 365 * 5/3;//24*365*5;

    const size_t n_scn=83;
    const size_t n_obj =2;
    const size_t n_ts=n_obj*2*n_scn;

    vector<point_ts<time_axis::fixed_dt>> ftsv;
    arma::mat a_mat(n,n_ts);
    time_axis::fixed_dt fta(t, dt, n);
    time_axis::generic_dt gta{t,dt*24,size_t(n/24)};
    const auto stair_case=ts_point_fx::POINT_AVERAGE_VALUE;
    map<string,apoint_ts> rtsv;
    ts_vector_t stsv;
    for(size_t i=0;i<n_ts;++i) {
        rtsv[to_string(i)] = apoint_ts(fta,i*10.0,stair_case);
        stsv.emplace_back(to_string(i));
        ftsv.emplace_back(fta,i*10.0,stair_case);
        for(size_t t =0;t<n;++t)
            a_mat(t,i) = i*10.0;
    }
    ts_vector_t tsv;
    for(size_t i =0; i<n_scn;++i) {
        apoint_ts sum;
        for(size_t j=0;j<n_obj;++j) {
            size_t p = i*(2*n_obj) + (2*j);
            size_t c = p+1;
            apoint_ts eff=3.0*(stsv[p]-stsv[c]);
            if(j==0)
                sum = eff;
            else
                sum = sum + eff;
        }
        tsv.emplace_back(sum);
    }
    tsv = 1000.0*(tsv.average(gta));


    //-- establish compute binding time
    auto t0 = timing::now();
    size_t bind_count{0};
    for(auto&sts:tsv) {
        auto ts_refs=sts.find_ts_bind_info();
        for(auto& bi:ts_refs) {
            bi.ts.bind(rtsv[bi.reference]);
            bind_count++;
        }
    }
    for(auto&sts:tsv)
        sts.do_bind();

    auto t1 =timing::now();
    //-- establish timing for apoint_ts eval.
    auto xmr = deflate_ts_vector<point_ts<time_axis::generic_dt>>(tsv);
    auto t2 = timing::now();

    FAST_CHECK_EQ(xmr.size(),n_scn);
    FAST_CHECK_EQ(bind_count,n_ts);
    cout<<"bind phase n_ts= "<<n_ts<<", n="<<n<<", time="<<double(elapsed_us(t0,t1))/1000.0<<"ms\n";
    cout<<"eval phase n_ts= "<<n_ts<<", n="<<n<<", time="<<double(elapsed_us(t1,t2))/1000.0<<"ms ->"
    << double(n*n_ts*(2+1))/(elapsed_us(t1,t2)/1e6)/1e6<<" mops/s \n";

}

TEST_CASE("dtss_container_wrapping") {
    auto tmpdir = (fs::temp_directory_path()/fs::unique_path());

    shyft::core::calendar utc;

    using shyft::time_axis::generic_dt;
    using shyft::time_series::point_ts;
    using shyft::time_series::ts_point_fx;
    // -----
    using shyft::dtss::ts_db;
    using cwrp_t = shyft::dtss::container_wrapper<ts_db>;

    SUBCASE("dispatch save through container") {
        std::string ts_name{ "test" };

        generic_dt ta{ utc.time(2002, 2, 2), calendar::HOUR, 24 };
        point_ts<generic_dt> ts{ ta, 15., shyft::time_series::ts_point_fx::POINT_INSTANT_VALUE };

        cwrp_t container{ std::make_unique<ts_db>(tmpdir.string()) };

        container.save(ts_name, ts);
        FAST_REQUIRE_UNARY( fs::is_regular_file(tmpdir/ts_name) );
    }

    SUBCASE("dispatch read through container") {
        std::string ts_name{ "test" };

        generic_dt ta{ utc.time(2002, 2, 2), calendar::HOUR, 24 };
        point_ts<generic_dt> ts{ ta, 15., shyft::time_series::ts_point_fx::POINT_INSTANT_VALUE };

        cwrp_t container{ std::make_unique<ts_db>(tmpdir.string()) };

        container.save(ts_name, ts);
        FAST_CHECK_UNARY( fs::is_regular_file(tmpdir/ts_name) );

        auto ts_new = container.read(ts_name, ta.total_period());
        TS_ASSERT_EQUALS( ts_new, ts );
    }

    SUBCASE("dispatch remove through container") {
        std::string ts_name{ "test" };

        generic_dt ta{ utc.time(2002, 2, 2), calendar::HOUR, 24 };
        point_ts<generic_dt> ts{ ta, 15., shyft::time_series::ts_point_fx::POINT_INSTANT_VALUE };

        cwrp_t container{ std::make_unique<ts_db>(tmpdir.string()) };

        container.save(ts_name, ts);
        FAST_CHECK_UNARY( fs::is_regular_file(tmpdir/ts_name) );

        container.remove(ts_name);
        FAST_REQUIRE_UNARY_FALSE( fs::exists(tmpdir/ts_name) );
    }

    SUBCASE("dispatch get_ts_info through container") {
        std::string ts_name{ "test" };

        generic_dt ta{ utc.time(2002, 2, 2), calendar::HOUR, 24 };
        point_ts<generic_dt> ts{ ta, 15., shyft::time_series::ts_point_fx::POINT_INSTANT_VALUE };

        cwrp_t container{ std::make_unique<ts_db>(tmpdir.string()) };

        container.save(ts_name, ts);
        FAST_CHECK_UNARY( fs::is_regular_file(tmpdir/ts_name) );

        auto info = container.get_ts_info(ts_name);
        FAST_REQUIRE_EQ( info.name, ts_name );
    }

    SUBCASE("dispatch find through container") {
        std::string ts_name{ "test" };

        generic_dt ta{ utc.time(2002, 2, 2), calendar::HOUR, 24 };
        point_ts<generic_dt> ts{ ta, 15., shyft::time_series::ts_point_fx::POINT_INSTANT_VALUE };

        cwrp_t container{ std::make_unique<ts_db>(tmpdir.string()) };

        container.save(ts_name, ts);
        FAST_CHECK_UNARY( fs::is_regular_file(tmpdir/ts_name) );

        auto info = container.find(ts_name);
        FAST_REQUIRE_EQ( info.size(), 1 );
        FAST_REQUIRE_EQ( info[0].name, ts_name );
    }
}


namespace {

struct query_test_dtss_container {

    using ts_info = shyft::dtss::ts_info;
    using gta_t = shyft::time_axis::generic_dt;
    using gts_t = shyft::time_series::point_ts<gta_t>;
    using queries_t = std::map<std::string, std::string>;

    std::string root;

    query_test_dtss_container() = default;
    ~query_test_dtss_container() = default;

    query_test_dtss_container(const std::string & root_dir) : root{ root_dir } {}

    query_test_dtss_container(const query_test_dtss_container &) = default;
    query_test_dtss_container & operator=(const query_test_dtss_container &) = default;

    query_test_dtss_container(query_test_dtss_container &&) = default;
    query_test_dtss_container & operator=(query_test_dtss_container &&) = default;

    /*  Container API
     * =============== */

    void save(const std::string & fn, const gts_t & ts, bool overwrite = true, const queries_t & queries = queries_t{}) const {

        FAST_CHECK_UNARY( fn.find("?") == fn.npos );  // there should be no query part in the filename

        auto q = queries.find("my_query");
        FAST_REQUIRE_NE( q, queries.cend() );
        FAST_CHECK_EQ( q->first, std::string{"my_query"} );
        FAST_CHECK_EQ( q->second, std::string{"some_value"} );

        FAST_CHECK_EQ( queries.find("removed"), queries.cend() );
    }

    gts_t read(const std::string & fn, core::utcperiod p, const queries_t & queries = queries_t{}) const {

        FAST_CHECK_UNARY( fn.find("?") == fn.npos );  // there should be no query part in the filename

        auto q = queries.find("my_query");
        FAST_REQUIRE_NE( q, queries.cend() );
        FAST_CHECK_EQ( q->first, std::string{"my_query"} );
        FAST_CHECK_EQ( q->second, std::string{"some_value"} );

        FAST_CHECK_EQ( queries.find("removed"), queries.cend() );
        return gts_t{};
    }

    void remove(const std::string & fn, const queries_t & queries = queries_t{}) const {

        FAST_CHECK_UNARY( fn.find("?") == fn.npos );  // there should be no query part in the filename

        auto q = queries.find("my_query");
        FAST_REQUIRE_NE( q, queries.cend() );
        FAST_CHECK_EQ( q->first, std::string{"my_query"} );
        FAST_CHECK_EQ( q->second, std::string{"some_value"} );

        FAST_CHECK_EQ( queries.find("removed"), queries.cend() );
    }

    ts_info get_ts_info(const std::string & fn, const queries_t & queries = queries_t{}) const {

        FAST_CHECK_UNARY( fn.find("?") == fn.npos );  // there should be no query part in the filename

        auto q = queries.find("my_query");
        FAST_REQUIRE_NE( q, queries.cend() );
        FAST_CHECK_EQ( q->first, std::string{"my_query"} );
        FAST_CHECK_EQ( q->second, std::string{"some_value"} );

        FAST_CHECK_EQ( queries.find("removed"), queries.cend() );

        return ts_info{};
    }

    std::vector<ts_info> find(const std::string & match, const queries_t & queries = queries_t{}) const {
        auto q = queries.find("my_query");
        FAST_REQUIRE_NE( q, queries.cend() );
        FAST_CHECK_EQ( q->first, std::string{"my_query"} );
        FAST_CHECK_EQ( q->second, std::string{"some_value"} );

        FAST_CHECK_EQ( queries.find("removed"), queries.cend() );

        return std::vector<ts_info>{};
    };
};

struct query_test_dtss_dispatcher {

    using queries_t = std::map<std::string, std::string>;
    using container_wrapper_t = dtss::container_wrapper<query_test_dtss_container>;

    inline static const std::string container_query{ "container" };
    inline static const std::array<std::string, 2> remove_queries{{ "removed", container_query }};

    static void create_container(
        const std::string & container_name,
        const std::string & container_type,
        const std::string & root_path,
        shyft::dtss::server<query_test_dtss_dispatcher> & dtss_server
    ) {
        if ( container_type.empty() || container_type == "test" ) {
            dtss_server.container[std::string{"TEST_"} + container_name] = container_wrapper_t{ std::make_unique<query_test_dtss_container>(root_path) };
        } else {
            throw std::runtime_error{ std::string{"Cannot construct unknown container type: "} + container_type };
        }
    }

    static container_wrapper_t & get_container(
        const std::string & container_name, const std::string & container_query,
        dtss::server<query_test_dtss_dispatcher> & dtss_server
    ) {
        decltype(dtss_server.container)::iterator f;
        if ( container_query.empty() || container_query == "test" ) {
            f = dtss_server.container.find(std::string{"TEST_"} + container_name);
        } else {
            throw std::runtime_error{ std::string{"Cannot construct unknown container type: "} + container_query };
        }

        if( f == std::end(dtss_server.container) )
            throw runtime_error(std::string{"Failed to find shyft container: "} + container_name);

        return f->second;
    }
};

}

TEST_CASE("dtss_server_query_to_containers") {

    shyft::core::calendar utc;

    using shyft::time_axis::generic_dt;
    using shyft::time_series::dd::apoint_ts;
    using shyft::time_series::ts_point_fx;
    using ts_vector_t = shyft::time_series::dd::ats_vector;
    using test_server_t = shyft::dtss::server<query_test_dtss_dispatcher>;
    using queries_t = std::map<std::string, std::string>;

    const std::string container{"container"};

    test_server_t srv{};
    srv.add_container(container, "foo/bar/00");  // same as "first"
    const int port =srv.start_server();

    shyft::dtss::client cli{ std::string{"localhost:"} + std::to_string(port) };

    /* NOTE: `remove` and `get_ts_info` are not called from the server...
     */

    SUBCASE("Save using url's with queries") {
        core::utctime t = utc.time(2016, 1, 1);
        core::utctimespan dt = core::deltahours(1);
        size_t n = 10;
        generic_dt ta{ t, dt, n };

        ts_vector_t vec{};
        for( std::size_t i = 0; i < 10; ++i ) {
            vec.emplace_back(
                shyft::dtss::shyft_url(container, to_string(i), queries_t{ {"my_query", "some_value"}, {"removed", "value"} }),
                apoint_ts{ ta, i*10.0, ts_point_fx::POINT_INSTANT_VALUE }
            );
        }

        cli.store_ts(vec, false, false);
    }

    SUBCASE("Read using url's with queries") {
        utcperiod period{ utc.time(2016, 1, 1), utc.time(2017, 1, 1) };

        ts_vector_t vec{};
        for( std::size_t i = 0; i < 10; ++i ) {
            vec.emplace_back(
                apoint_ts{ shyft::dtss::shyft_url(container, to_string(i), queries_t{ {"my_query", "some_value"}, {"removed", "value"} }) }
            );
        }

        auto res = cli.evaluate(vec, period, false, false);
    }

    SUBCASE("Find using url's with queries") {
        cli.find(shyft::dtss::shyft_url(container, "/path/to/something", queries_t{ {"my_query", "some_value"}, {"removed", "value"} }));
    }
}


namespace {

template < int Idx >
struct multicontainer_test_dtss_container {

    using ts_info = shyft::dtss::ts_info;
    using gta_t = shyft::time_axis::generic_dt;
    using gts_t = shyft::time_series::point_ts<gta_t>;
    using queries_t = std::map<std::string, std::string>;

    std::string root;

    multicontainer_test_dtss_container() = default;
    ~multicontainer_test_dtss_container() = default;

    multicontainer_test_dtss_container(const std::string & root_dir) : root{ root_dir } {}

    multicontainer_test_dtss_container(const multicontainer_test_dtss_container &) = default;
    multicontainer_test_dtss_container & operator=(const multicontainer_test_dtss_container &) = default;

    multicontainer_test_dtss_container(multicontainer_test_dtss_container &&) = default;
    multicontainer_test_dtss_container & operator=(multicontainer_test_dtss_container &&) = default;

    /*  Container API
    * =============== */

    void save(const std::string & fn, const gts_t & ts, bool overwrite = true, const queries_t & queries = queries_t{}) const {

        FAST_CHECK_UNARY( fn.find("?") == fn.npos );  // there should be no query part in the filename

        auto q = queries.find("value");
        FAST_REQUIRE_NE( q, queries.cend() );
        FAST_CHECK_EQ( q->first, std::string{"value"} );
        FAST_CHECK_EQ( q->second, std::to_string(Idx) );
    }

    gts_t read(const std::string & fn, core::utcperiod p, const queries_t & queries = queries_t{}) const {

        FAST_CHECK_UNARY( fn.find("?") == fn.npos );  // there should be no query part in the filename

        auto q = queries.find("value");
        FAST_REQUIRE_NE( q, queries.cend() );
        FAST_CHECK_EQ( q->first, std::string{"value"} );
        FAST_CHECK_EQ( q->second, std::to_string(Idx) );
        return gts_t{};
    }

    void remove(const std::string & fn, const queries_t & queries = queries_t{}) const {

        FAST_CHECK_UNARY( fn.find("?") == fn.npos );  // there should be no query part in the filename

        auto q = queries.find("value");
        FAST_REQUIRE_NE( q, queries.cend() );
        FAST_CHECK_EQ( q->first, std::string{"value"} );
        FAST_CHECK_EQ( q->second, std::to_string(Idx) );
    }

    ts_info get_ts_info(const std::string & fn, const queries_t & queries = queries_t{}) const {

        FAST_CHECK_UNARY( fn.find("?") == fn.npos );  // there should be no query part in the filename

        auto q = queries.find("value");
        FAST_REQUIRE_NE( q, queries.cend() );
        FAST_CHECK_EQ( q->first, std::string{"value"} );
        FAST_CHECK_EQ( q->second, std::to_string(Idx) );

        return ts_info{};
    }

    std::vector<ts_info> find(const std::string & match, const queries_t & queries = queries_t{}) const {
        auto q = queries.find("value");
        FAST_REQUIRE_NE( q, queries.cend() );
        FAST_CHECK_EQ( q->first, std::string{"value"} );
        FAST_CHECK_EQ( q->second, std::to_string(Idx) );

        return std::vector<ts_info>{};
    };
};

struct multicontainer_test_dtss_dispatcher {

    using queries_t = std::map<std::string, std::string>;
    using container_wrapper_t = dtss::container_wrapper<
        multicontainer_test_dtss_container<0>, multicontainer_test_dtss_container<1> >;

    inline static const std::string container_query{ "container" };
    inline static const std::array<std::string, 2> remove_queries{{ "removed", container_query }};

    static void create_container(
        const std::string & container_name,
        const std::string & container_type,
        const std::string & root_path,
        shyft::dtss::server<multicontainer_test_dtss_dispatcher> & dtss_server
    ) {
        if ( container_type.empty() || container_type == "first" ) {
            dtss_server.container[std::string{"FIRST_"} + container_name] = container_wrapper_t{ std::make_unique<multicontainer_test_dtss_container<0>>(root_path) };
        } else if ( container_type == "second" ) {
            dtss_server.container[std::string{"SECOND_"} + container_name] = container_wrapper_t{ std::make_unique<multicontainer_test_dtss_container<1>>(root_path) };
        } else {
            throw std::runtime_error{ std::string{"Cannot construct unknown container type: "} + container_type };
        }
    }

    static container_wrapper_t & get_container(
        const std::string & container_name, const std::string & container_query,
        dtss::server<multicontainer_test_dtss_dispatcher> & dtss_server
    ) {
        decltype(dtss_server.container)::iterator f;
        if ( container_query.empty() || container_query == "first" ) {
            f = dtss_server.container.find(std::string{"FIRST_"} + container_name);
        } else if ( container_query == "second" ) {
            f = dtss_server.container.find(std::string{"SECOND_"} + container_name);
        } else {
            throw std::runtime_error{ std::string{"Cannot construct unknown container type: "} + container_query };
        }

        if( f == std::end(dtss_server.container) )
            throw runtime_error(std::string{"Failed to find shyft container: "} + container_name);

        return f->second;
    }
};

}

TEST_CASE("dtss_server_with_multiple_containers") {

    shyft::core::calendar utc;

    using shyft::time_axis::generic_dt;
    using shyft::time_series::dd::apoint_ts;
    using shyft::time_series::ts_point_fx;
    using ts_vector_t = shyft::time_series::dd::ats_vector;
    using test_server_t = shyft::dtss::server<multicontainer_test_dtss_dispatcher>;
    using queries_t = std::map<std::string, std::string>;

    test_server_t srv{};
    srv.add_container("container", "foo/bar/first", "first");
    srv.add_container("container", "foo/bar/second", "second");
    const int port=srv.start_server();

    shyft::dtss::client cli{ std::string{"localhost:"} + std::to_string(port) };

    /* NOTE: `remove` and `get_ts_info` are not called from the server...
    */

    SUBCASE("Save dispatched to the correct container") {
        core::utctime t = utc.time(2016, 1, 1);
        core::utctimespan dt = core::deltahours(1);
size_t n = 10;
generic_dt ta{ t, dt, n };

ts_vector_t vec{};
vec.emplace_back(
    shyft::dtss::shyft_url("container", "foo/bar/baz", queries_t{ {"value", "0"}, {"container", "first"} }),
    apoint_ts{ ta, 0 * 10.0, ts_point_fx::POINT_INSTANT_VALUE });
vec.emplace_back(
    shyft::dtss::shyft_url("container", "foo/bar/baz", queries_t{ {"value", "1"}, {"container", "second"} }),
    apoint_ts{ ta, 1 * 10.0, ts_point_fx::POINT_INSTANT_VALUE });

cli.store_ts(vec, false, false);
    }

    SUBCASE("Read dispatched to the correct container") {
        utcperiod period{ utc.time(2016, 1, 1), utc.time(2017, 1, 1) };

        ts_vector_t vec{};
        vec.emplace_back(shyft::dtss::shyft_url("container", "foo/bar/baz", queries_t{ {"value", "0"}, {"container", "first"} }));
        vec.emplace_back(shyft::dtss::shyft_url("container", "foo/bar/baz", queries_t{ {"value", "1"}, {"container", "second"} }));

        auto res = cli.evaluate(vec, period, false, false);
    }

    SUBCASE("Find dispatched to the correct container") {
        cli.find(shyft::dtss::shyft_url("container", "foo/bar/baz", queries_t{ {"value", "0"}, {"container", "first"} }));
        cli.find(shyft::dtss::shyft_url("container", "foo/bar/baz", queries_t{ {"value", "1"}, {"container", "second"} }));
    }
}

}
