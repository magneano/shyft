/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#include <map>
#include <algorithm>
#include <memory>
#include <utility>
#include <functional>
#include <cstring>
#include <regex>
#include <future>
#include <utility>
#include <chrono>

#include "dtss_client.h"
#include "dtss_url.h"
#include "dtss_msg.h"

#include "core_serialization.h"
#include "core_archive.h"
#include "expression_serialization.h"

#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>

namespace shyft {
namespace dtss {
using std::runtime_error;
using std::vector;
using std::string;
using std::int64_t;
using std::exception;
using std::future;
using std::min;
using std::chrono::milliseconds;
using std::this_thread::sleep_for;

using shyft::core::core_iarchive;
using shyft::core::core_oarchive;
using shyft::time_series::dd::apoint_ts;
using shyft::time_series::dd::aref_ts;
using shyft::time_series::dd::gta_t;
using shyft::time_series::dd::expression_compressor;
using shyft::time_series::statistics_property;
using dlib::socket_error;

void srv_connection::open(int timeout_ms) {
	io->open(host_port, max(timeout_ms, this->timeout_ms));
    is_open=true;
}
void srv_connection::close(int timeout_ms) {
    is_open=false; // if next line throws, consider closed anyway
	io->close(max(timeout_ms, this->timeout_ms));
}
void srv_connection::reopen(int timeout_ms) {
	io->open(host_port, max(timeout_ms, this->timeout_ms));
    is_open=true;
}

/** 
 * The socket connection created from the client stays connected from the first use until explicitely 
 * closed by client.
 * The life time of the underlying socket *after* the close is ~ 120 seconds(windows, similar linux).
 * If the dtss-server is restarted we would enjoy that this layer 'auto-repair' with the dtss-server
 * if possible.
 * This routine ensures that this can happen for any io function sequence F.
 * note that there is a requirement that f(sc) can be invoced several times 
 * without (unwanted) side-effects.
 */
template <class F>
void do_io_with_repair_and_retry(srv_connection&sc, F&&f) {
	for (int retry = 0; retry < 3; ++retry) {
		try {
			f(sc);
			return;
		} catch (const dlib::socket_error&) {
			sc.reopen();
		}
	}
	throw std::runtime_error("Failed to establish connection with " + sc.host_port);
}

/** helper-class to enable 'autoconnect', lazy connect. */
struct scoped_connect {
    client& c;
    scoped_connect (client& c):c(c){
        bool rethrow=false;
        runtime_error rt_re("");
        for(auto&sc:c.srv_con) {
            if (!sc.is_open ) { // auto-connect, and put some effort into succeeding
                bool attempt_more=false;
                int retry_count=5;
                do {
                    try {
                        sc.open(); // either we succeed and finishes, or we get an exception
                        attempt_more=false;
                    } catch(const socket_error&se) {// capture socket error that might go away if we stay alive a little bit more
                        if(--retry_count>0 && strstr(se.what(),"unable to connect") ) {
                            attempt_more=true;
                            sleep_for(milliseconds(100));
                        } else {
                            rt_re=runtime_error(se.what());
                            rethrow=true;
                            attempt_more=false;
                        }
                    } catch(const exception&re) { // just give up this type of error,
                            rt_re=runtime_error(re.what());
                            rethrow=true;
                            attempt_more=false;
                    }
                } while(attempt_more);
            }
        }
        if(rethrow)
            throw rt_re;
    }
    ~scoped_connect() noexcept(false) {
		// we don't disconnect, rather we try to fix broken connects.
    }
    scoped_connect (const scoped_connect&) = delete;
    scoped_connect ( scoped_connect&&) = delete;
    scoped_connect& operator=(const scoped_connect&) = delete;
    scoped_connect()=delete;
    scoped_connect& operator=(scoped_connect&&)=delete;
};

client::client ( const string& host_port, bool auto_connect, int timeout_ms )
    :auto_connect(auto_connect)
{
    srv_con.push_back(srv_connection{make_unique<dlib::iosockstream>(),host_port,timeout_ms});
}

client::client(const vector<string>& host_ports,bool auto_connect,int timeout_ms):auto_connect(auto_connect) {
    if(host_ports.size()==0)
        throw runtime_error("host_ports must contain at least one element");
    for(const auto &hp:host_ports) {
        srv_con.push_back(srv_connection{make_unique<dlib::iosockstream>(),hp,timeout_ms});
    }
}

void client::reopen(int timeout_ms) {
    for(auto&sc:srv_con)
        sc.reopen(timeout_ms);
}

void client::close(int timeout_ms) {
    bool rethrow = false;
    runtime_error rt_re("");
    for(auto&sc:srv_con) {
        try {
            sc.close(timeout_ms);
        } catch (const exception &re) { // ensure we try to close all
            rt_re=runtime_error(re.what());
            rethrow=true;
        }
    }
    if(rethrow)
        throw rt_re;
}

vector<apoint_ts>
client::percentiles(ts_vector_t const& tsv, utcperiod p, gta_t const&ta, const vector<int64_t>& percentile_spec,bool use_ts_cached_read,bool update_ts_cache) {
    if (tsv.size() == 0)
        throw runtime_error("percentiles requires a source ts-vector with more than 0 time-series");
    if (percentile_spec.size() == 0)
        throw std::runtime_error("percentile function require more than 0 percentiles specified");
    if (!p.valid())
        throw std::runtime_error("percentiles require a valid period-specification");
    if (ta.size() == 0)
        throw std::runtime_error("percentile function require a time-axis with more than 0 steps");


    if(srv_con.size()==1 || tsv.size()< srv_con.size()) {
        scoped_connect ac(*this);
		ts_vector_t r;
		do_io_with_repair_and_retry(srv_con[0],
			[this,&r,&ta,&tsv,&percentile_spec,&use_ts_cached_read,&update_ts_cache,&p](srv_connection&sc) {
				dlib::iosockstream& io = *sc.io;
				msg::write_type(compress_expressions ? message_type::EVALUATE_EXPRESSION_PERCENTILES : message_type::EVALUATE_TS_VECTOR_PERCENTILES, io);
				core_oarchive oa(io, core_arch_flags);
				oa << p;
				if (compress_expressions) {
					oa << expression_compressor::compress(tsv);
				} else {
					oa << tsv;
				}
				oa << ta << percentile_spec << use_ts_cached_read << update_ts_cache;
				auto response_type = msg::read_type(io);
				if (response_type == message_type::SERVER_EXCEPTION) {
					auto re = msg::read_exception(io);
					throw re;
				} else if (response_type == message_type::EVALUATE_TS_VECTOR_PERCENTILES) {
					core_iarchive ia(io, core_arch_flags);
					ia >> r;
				} else {
					throw std::runtime_error(std::string("Got unexpected response:") + std::to_string((int)response_type));
				}
			}
		);
		return r;
    } else {
        vector<int> p_spec;
        bool can_do_server_side_average=true; // in case we are searching for min-max extreme, we can not do server-side average
        for(size_t i=0;i<percentile_spec.size();++i) {
            p_spec.push_back(int(percentile_spec[i]));
            if(percentile_spec[i] == statistics_property::MIN_EXTREME || percentile_spec[i]==statistics_property::MAX_EXTREME)
                can_do_server_side_average=false;
        }
        auto atsv = evaluate(can_do_server_side_average?tsv.average(ta):tsv,p,use_ts_cached_read,update_ts_cache); // get the result we can do percentiles on
        return shyft::time_series::dd::percentiles(atsv, ta, p_spec);
    }
}

std::vector<apoint_ts>
client::evaluate(ts_vector_t const& tsv, utcperiod p,bool use_ts_cached_read,bool update_ts_cache) {
    if (tsv.size() == 0)
        throw std::runtime_error("evaluate requires a source ts-vector with more than 0 time-series");
    if (!p.valid())
        throw std::runtime_error("percentiles require a valid period-specification");
    scoped_connect ac(*this);
    // local lambda to ensure one definition of communication with the server
    auto eval_io = [this] (srv_connection &sc,const ts_vector_t& tsv,const utcperiod& p,bool use_ts_cached_read,bool update_ts_cache) ->ts_vector_t {
		ts_vector_t r;
		do_io_with_repair_and_retry(sc,
			[this,&r,&tsv,&p,&use_ts_cached_read,&update_ts_cache](srv_connection &sc) {
				dlib::iosockstream& io = *sc.io;
				msg::write_type(compress_expressions ? message_type::EVALUATE_EXPRESSION : message_type::EVALUATE_TS_VECTOR, io); {
					core_oarchive oa(io, core_arch_flags);
					oa << p;
					if (compress_expressions) { // notice that we stream out all in once here
						// .. just in case the destruction of the compressed expr take time (it could..)
						oa << expression_compressor::compress(tsv) << use_ts_cached_read << update_ts_cache;
					} else {
						oa << tsv << use_ts_cached_read << update_ts_cache;
					}
				}
				auto response_type = msg::read_type(io);
				if (response_type == message_type::SERVER_EXCEPTION) {
					auto re = msg::read_exception(io);
					throw re;
				} else if (response_type == message_type::EVALUATE_TS_VECTOR) {
					core_iarchive ia(io, core_arch_flags);
					ia >> r;
				} else {
					throw std::runtime_error(std::string("Got unexpected response:") + std::to_string((int)response_type));
				}
			}
		);
		return r;
    };

    if(srv_con.size()==1 || tsv.size() == 1) { // one server, or just one ts, do it easy
        return eval_io(srv_con[0],tsv,p,use_ts_cached_read,update_ts_cache);
    } else {
        ts_vector_t rt(tsv.size()); // make place for the result contributions from threads
        // lamda to eval partition on server
        auto eval_partition= [&rt,&tsv,&eval_io,p,use_ts_cached_read,update_ts_cache]
            (srv_connection& sc,size_t i0, size_t n) {
                ts_vector_t ptsv;ptsv.reserve(tsv.size());
                for(size_t i=i0;i<i0+n;++i) ptsv.push_back(tsv[i]);
                auto pt = eval_io(sc,ptsv,p,use_ts_cached_read,update_ts_cache);
                for(size_t i=0;i<pt.size();++i)
                    rt[i0+i]=pt[i];// notice how se stash data into common result variable, safe because each thread using it's own preallocated space
        };
        size_t n_ts_pr_server = tsv.size()/srv_con.size(); 
        size_t remainder = tsv.size() - n_ts_pr_server*srv_con.size();
        vector<future<void>> calcs;
        size_t b=0;
        for(size_t i=0;i<srv_con.size() && b<tsv.size();++i) {
            size_t e = b+n_ts_pr_server;
            if(remainder>0) {
                e+=1;// spread remainder equally between servers
                --remainder;
            } else if (e>tsv.size()) {
                e=tsv.size();
            }
                
            size_t n = e-b;
            calcs.push_back(std::async(
                            std::launch::async,
                            [this,i,b,n,&eval_partition] () {
                                eval_partition(srv_con[i],b,n);
                            }
                           )
                  );
            b=e;
        }
        for (auto &f : calcs)
            f.get();

        return rt;
    }
}

void
client::store_ts(const ts_vector_t &tsv, bool overwrite_on_write, bool cache_on_write) {
    if (tsv.size() == 0)
        return; //trivial and considered valid case
                // verify that each member of tsv is a gpoint_ts
    for (auto const &ats : tsv) {
        auto rts = dynamic_cast<aref_ts*>(ats.ts.get());
        if (!rts) throw std::runtime_error(std::string("attempt to store a null ts"));
        if (rts->needs_bind()) throw std::runtime_error(std::string("attempt to store unbound ts:") + rts->id);
    }
    scoped_connect ac(*this);
	do_io_with_repair_and_retry(srv_con[0],
		[this,&tsv,overwrite_on_write,cache_on_write](srv_connection&sc) {
			dlib::iosockstream& io = *(sc.io);
			msg::write_type(message_type::STORE_TS, io);
			{
				core_oarchive oa(io, core_arch_flags);
				oa << tsv << overwrite_on_write << cache_on_write;
			}
			auto response_type = msg::read_type(io);
			if (response_type == message_type::SERVER_EXCEPTION) {
				auto re = msg::read_exception(io);
				throw re;
			} else if (response_type == message_type::STORE_TS) {
				return;
			}
			throw std::runtime_error(std::string("Got unexpected response:") + std::to_string((int)response_type));
		}
	);
}

void
client::merge_store_ts(const ts_vector_t &tsv,bool cache_on_write) {
    if (tsv.size() == 0)
        return; //trivial and considered valid case
                // verify that each member of tsv is a gpoint_ts
    for (auto const &ats : tsv) {
        auto rts = dynamic_cast<aref_ts*>(ats.ts.get());
        if (!rts) throw std::runtime_error(std::string("attempt to store a null ts"));
        if (rts->needs_bind()) throw std::runtime_error(std::string("attempt to store unbound ts:") + rts->id);
    }
    scoped_connect ac(*this);
	do_io_with_repair_and_retry(srv_con[0],
		[this, &tsv, cache_on_write](srv_connection&sc) {
			dlib::iosockstream& io = *(sc.io);
			msg::write_type(message_type::MERGE_STORE_TS, io);
			{
				core_oarchive oa(io, core_arch_flags);
				oa << tsv << cache_on_write;
			}
			auto response_type = msg::read_type(io);
			if (response_type == message_type::SERVER_EXCEPTION) {
				auto re = msg::read_exception(io);
				throw re;
			} else if (response_type == message_type::MERGE_STORE_TS) {
				return;
			}
			throw std::runtime_error(std::string("Got unexpected response:") + std::to_string((int)response_type));
		}
	);
}

ts_info_vector_t
client::find(const std::string& search_expression) {
    scoped_connect ac(*this);
	ts_info_vector_t r;
	do_io_with_repair_and_retry(srv_con[0],
		[this,&r, &search_expression](srv_connection&sc) {
			auto& io = *(sc.io);
			msg::write_type(message_type::FIND_TS, io); {
				msg::write_string(search_expression, io);
			}
			auto response_type = msg::read_type(io);
			if (response_type == message_type::SERVER_EXCEPTION) {
				auto re = msg::read_exception(io);
				throw re;
			} else if (response_type == message_type::FIND_TS) {
					core_iarchive ia(io, core_arch_flags);
					ia >> r;
			} else {
				throw std::runtime_error(std::string("Got unexpected response:") + std::to_string((int)response_type));
			}
		}
	);
	return r;
}

ts_info
client::get_ts_info(const std::string & ts_url) {
    scoped_connect ac(*this);
    ts_info tsi;
    do_io_with_repair_and_retry(srv_con[0],
        [this, &tsi, &ts_url](srv_connection & sc) {
            auto& io = *(sc.io);
            msg::write_type(message_type::GET_TS_INFO, io); {
                msg::write_string(ts_url, io);
            }
            auto response_type = msg::read_type(io);
            if (response_type == message_type::SERVER_EXCEPTION) {
                auto re = msg::read_exception(io);
                throw re;
            } else if (response_type == message_type::GET_TS_INFO) {
                core_iarchive ia(io, core_arch_flags);
                ia >> tsi;
            } else {
                throw std::runtime_error(std::string("Got unexpected response:") + std::to_string((int)response_type));
            }
        }
    );
    return tsi;
}

void
client::remove(const string & ts_url) {
	scoped_connect ac(*this);
	do_io_with_repair_and_retry(srv_con[0],
		[this, &ts_url](srv_connection&sc) {
			auto& io = *(sc.io);
			msg::write_type(message_type::REMOVE_TS, io);{
				msg::write_string(ts_url, io);
			}
			auto response_type = msg::read_type(io);
			if (response_type == message_type::SERVER_EXCEPTION) {
				auto re = msg::read_exception(io);
				throw re;
			} else if (response_type == message_type::REMOVE_TS) {
				return;
			}
			throw std::runtime_error(std::string("Got unexpected response:") + std::to_string((int)response_type));
		}
	);
}

void
client::cache_flush() {
    scoped_connect ac(*this);
    for(auto& sc:srv_con) {
		do_io_with_repair_and_retry(sc,
			[this](srv_connection&sc) {
				auto& io = *(sc.io);
				msg::write_type(message_type::CACHE_FLUSH, io);
				auto response_type = msg::read_type(io);
				if (response_type == message_type::SERVER_EXCEPTION) {
					auto re = msg::read_exception(io);
					throw re;
				} else if (response_type != message_type::CACHE_FLUSH) {
					throw std::runtime_error(std::string("Got unexpected response:") + std::to_string((int)response_type));
				}
			}
		);
	}
}



cache_stats
client::get_cache_stats() {
    scoped_connect ac(*this);
    cache_stats s;
    for(auto& sc:srv_con) {
		do_io_with_repair_and_retry(sc, [&s](srv_connection&sc) {
			auto& io = *(sc.io);
			msg::write_type(message_type::CACHE_STATS, io);
			auto response_type = msg::read_type(io);
			if (response_type == message_type::CACHE_STATS) {
				cache_stats r;
				core_iarchive oa(io, core_arch_flags);
				oa >> r;
				s = s + r;
			} else if (response_type == message_type::SERVER_EXCEPTION) {
				auto re = msg::read_exception(io);
				throw re;
			} else {
				throw std::runtime_error(std::string("Got unexpected response:") + std::to_string((int)response_type));
			}
		}
		);
    }
    return s;
}

}
}
