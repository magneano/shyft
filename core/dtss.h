/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#pragma once

#include <cstdint>
#include <cstdio>
#include <string>
#include <string_view>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <memory>
#include <utility>
#include <functional>
#include <cstring>
#include <regex>
#include <variant>

#include "core_serialization.h"
#include "expression_serialization.h"
#include "core_archive.h"

#include <boost/functional/hash.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <dlib/server.h>
#include <dlib/iosockstream.h>
#include <dlib/logger.h>
#include <dlib/misc_api.h>

#include "core/time_series_dd.h"
#include "time_series_info.h"
#include "utctime_utilities.h"
#include "dtss_cache.h"
#include "dtss_url.h"
#include "dtss_msg.h"
#include "dtss_db.h"
#include "dtss_krls.h"


namespace shyft {
namespace dtss {

using shyft::core::utctime;
using shyft::core::utcperiod;
using shyft::core::utctimespan;
using shyft::core::no_utctime;
using shyft::core::calendar;
using shyft::core::deltahours;

using gta_t = shyft::time_axis::generic_dt;
using gts_t = shyft::time_series::point_ts<gta_t>;

using shyft::time_series::dd::apoint_ts;
using shyft::time_series::dd::gpoint_ts;
using shyft::time_series::dd::gts_t;
using shyft::time_series::dd::aref_ts;


// ========================================


using ts_vector_t = shyft::time_series::dd::ats_vector;
using ts_info_vector_t = std::vector<ts_info>;
using id_vector_t = std::vector<std::string>;
using read_call_back_t = std::function<ts_vector_t(const id_vector_t& ts_ids, utcperiod p)>;
using store_call_back_t = std::function<void(const ts_vector_t&)>;
using find_call_back_t = std::function<ts_info_vector_t(std::string search_expression)>;

using std::unique_ptr;
using std::make_unique;


// ========================================


/**
 */
template < class ... CIMPLs >  /* Container Implementations */
struct container_wrapper {

    using gta_t = shyft::time_axis::generic_dt;
    using gts_t = shyft::time_series::point_ts<gta_t>;
    using queries_t = std::map<std::string, std::string>;
    // -----
    using container_adt = std::variant<std::unique_ptr<CIMPLs>... >;
    // -----
    container_adt _container;

    container_wrapper() = default;
    ~container_wrapper() = default;
    // -----
    template < class CIMPL >
    container_wrapper(const std::unique_ptr<CIMPL> & c) : _container{ c } { }
    template < class CIMPL >
    container_wrapper(std::unique_ptr<CIMPL> && c) : _container{ std::forward<std::unique_ptr<CIMPL>>(c) } { }
    // -----
    container_wrapper(const container_wrapper &) = default;
    container_wrapper & operator=(const container_wrapper &) = default;
    // -----
    container_wrapper(container_wrapper &&) = default;
    container_wrapper & operator=(container_wrapper &&) = default;

    /* Container API
     * ------------- */

public:
    /** Save a time-series to the container. */
    void save(
        const std::string & tsid,
        const gts_t & ts,
        bool overwrite = true,
        const queries_t & queries = queries_t{}
    ) {
        std::visit([&](auto && var) {
            var->save(tsid, ts, overwrite, queries);
        }, _container);
    }
    /** Read a period from a time-series from the container. */
    gts_t read(
        const std::string & tsid,
        core::utcperiod period,
        const queries_t & queries = queries_t{}
    ) {
        return std::visit([&](auto && var) -> gts_t {
            return var->read(tsid, period, queries);
        }, _container);
    }
    /** Remove a time-series from the container. */
    void remove(
        const std::string & tsid,
        const queries_t & queries = queries_t{}
    ) {
        std::visit([&](auto && var) {
            return var->remove(tsid, queries);
        }, _container);
    }
    /** Get minimal info about a time-series stored in the container. */
    ts_info get_ts_info(
        const std::string & tsid,
        const queries_t & queries = queries_t{}
    ) {
        return std::visit([&](auto && var) -> ts_info {
            return var->get_ts_info(tsid, queries);
        }, _container);
    }
    /** Find minimal information about all time-series stored in the container matching a regex pattern. */
    std::vector<ts_info> find(
        const std::string & pattern,
        const queries_t & queries = queries_t{}
    ) {
        return std::visit([&](auto && var) -> std::vector<ts_info> {
            return var->find(pattern, queries);
        }, _container);
    }
};


/** @brief A dtss server with time-series server-side functions
 *
 * The dtss server listens on a port, receives messages, interpret them
 * and ship the response back to the client.
 *
 * Callbacks are provided for extending/delegating find/read_ts/store_ts,
 * as well as internal implementation of storing time-series
 * using plain binary files stored in containers(directory).
 *
 * Time-series are named with url's, and all request involving 'shyft://'
 * like
 *   shyft://<container>/<local_ts_name>
 * resolves to the internal implementation.
 *
 * @tparam CD ContainerDispatcher class that helps implementing static dispatch to special containers, like the KRLS for training expressions
 *
 */
template < class CD >
struct server : dlib::server_iostream {

    friend CD;  // the dispatcher have full access to the server

    using ts_cache_t = cache<apoint_ts_frag, apoint_ts>;
    using cwrp_t = typename CD::container_wrapper_t;

    // callbacks for extensions
    read_call_back_t bind_ts_cb;    ///< called to read non shyft:// unbound ts
    find_call_back_t find_ts_cb;    ///< called for all non shyft:// find operations
    store_call_back_t store_ts_cb;  ///< called for all non shyft:// store operations

    // shyft-internal implementation
    mutex c_mx;///< container mutex
    std::unordered_map<std::string, cwrp_t> container;  ///< mapping of internal shyft <container>
    ts_cache_t ts_cache{ 1000000 };  // default 1 mill ts in cache
    bool cache_all_reads{ false };
    bool can_remove{ false };

    // constructors
    server()=default;
    server(server&&)=delete;
    server(const server&) =delete;
    server& operator=(const server&)=delete;
    server& operator=(server&&)=delete;

    template < class CB >
    explicit server(CB && cb)
        : bind_ts_cb{ std::forward<CB>(cb) }  {
    }

    template < class RCB, class FCB>
    server(RCB && rcb, FCB && fcb)
        : bind_ts_cb{ std::forward<RCB>(rcb) }, find_ts_cb{ std::forward<FCB>(fcb) } {
    }

    template < class RCB, class FCB, class SCB >
    server(RCB && rcb, FCB && fcb, SCB && scb)
        : bind_ts_cb{ std::forward<RCB>(rcb) }, find_ts_cb{ std::forward<FCB>(fcb) },
          store_ts_cb{ std::forward<SCB>(scb) } {
    }

    ~server() =default;

    //-- container management
    void add_container(
        const std::string & container_name, const std::string & root_dir,
        std::string container_type = std::string{}
    ) {
        unique_lock<mutex> sl(c_mx);
        CD::create_container(container_name, container_type, root_dir, *this);
    }

    cwrp_t & internal(const std::string & container_name, const std::string & container_query = std::string{}) {
        return CD::get_container(container_name, container_query, *this);
    }

    /** start the server in background, return the listening port used in case it was set unspecified */
    int start_server() {
        if(get_listening_port()==0) {
            start_async();
            while(is_running()&& get_listening_port()==0) //because dlib do not guarantee that listening port is set
                std::this_thread::sleep_for(std::chrono::milliseconds(10)); // upon return, so we have to wait until it's done
        } else {
            start_async();
        }
        return get_listening_port();
    }
    //-- expose cache functions

    void add_to_cache(id_vector_t&ids, ts_vector_t& tss) { ts_cache.add(ids,tss);}
    void remove_from_cache(id_vector_t &ids) { ts_cache.remove(ids);}
    cache_stats get_cache_stats() { return ts_cache.get_cache_stats();}
    void clear_cache_stats() { ts_cache.clear_cache_stats();}
    void flush_cache() { return ts_cache.flush();}
    void set_cache_size(std::size_t max_size) { ts_cache.set_capacity(max_size);}
    void set_auto_cache(bool active) { cache_all_reads=active;}
    std::size_t get_cache_size() const {return ts_cache.get_capacity();}

    void set_can_remove(bool can_remove) { this->can_remove = can_remove; }

    ts_info_vector_t do_find_ts(const std::string& search_expression);
    ts_info do_get_ts_info(const std::string & ts_url);

    std::string extract_url(const apoint_ts&ats) const {
        auto rts = dynamic_pointer_cast<aref_ts>(ats.ts);
        if(rts)
            return rts->id;
        throw runtime_error("dtss store.extract_url:supplied type must be of type ref_ts");
    }

    void do_cache_update_on_write(const ts_vector_t&tsv,bool overwrite_on_write);

    void do_store_ts(const ts_vector_t & tsv, bool overwrite_on_write, bool cache_on_write);

    void do_merge_store_ts(const ts_vector_t & tsv, bool cache_on_write);
    /** @brief Read the time-series from providers for specified period
    *
    * @param ts_ids identifiers, url form, where shyft://.. is specially filtered
    * @param p the period to read
    * @param use_ts_cached_read allow reading results from already existing cached results
    * @param update_ts_cache when reading, also update the ts-cache with the results
    * \return read ts-vector in the order of the ts_ids
    */
    ts_vector_t do_read(const id_vector_t& ts_ids,utcperiod p,bool use_ts_cached_read,bool update_ts_cache);
    void do_remove_ts(const std::string & ts_url);
    void do_bind_ts(utcperiod bind_period, ts_vector_t& atsv,bool use_ts_cached_read,bool update_ts_cache);
    ts_vector_t do_evaluate_ts_vector(utcperiod bind_period, ts_vector_t& atsv,bool use_ts_cached_read,bool update_ts_cache,utcperiod clip_period);
    ts_vector_t do_evaluate_percentiles(utcperiod bind_period, ts_vector_t& atsv, gta_t const&ta,std::vector<int64_t> const& percentile_spec,bool use_ts_cached_read,bool update_ts_cache);

    // ref. dlib, all connection calls are directed here
    void on_connect(
        std::istream & in,
        std::ostream & out,
        const std::string & foreign_ip,
        const std::string & local_ip,
        unsigned short foreign_port,
        unsigned short local_port,
        dlib::uint64 connection_id
    );
};


// --------------------------------------------------------------------------------


struct standard_dtss_dispatcher {

    using queries_t = std::map<std::string, std::string>;

    /** Alias for the wrapping of containers. */
    using container_wrapper_t = container_wrapper<ts_db, krls_pred_db>;

    /** Query key used to specify container types. */
    inline static const std::string container_query{ "container" };

    /** Sequence of query keys to remove before passing the query map to the containers. */
    inline static const std::array<std::string, 1> remove_queries{{ container_query }};

    /** @brief Create a new container into the server.
     * @param  container_name  Name of the new container. Used to access the container.
     * @param  container_type  Identifier of the container type to add.
     * @param  root_path  Root path to where the container data is stored.
     * @param  dtss_server  Dtss server instance the container is added to.
     */
    static void create_container(
        const std::string & container_name,
        const std::string & container_type,
        const std::string & root_path,
        dtss::server<standard_dtss_dispatcher> & dtss_server
    ) {
        std::string server_container_name;
        if ( container_type.empty() || container_type == "ts_db" ) {
            server_container_name = container_name;
            dtss_server.container[server_container_name] = container_wrapper_t{ std::make_unique<ts_db>( root_path ) };
        } else if ( container_type == "krls" ) {
            server_container_name = std::string{"KRLS_"} + container_name;
            dtss_server.container[server_container_name] = container_wrapper_t {
                std::make_unique<krls_pred_db>(
                    root_path,
                    [&dtss_server](const std::string & tsid, utcperiod period, bool use_ts_cached_read, bool update_ts_cache) -> ts_vector_t {
                        id_vector_t id_vec{ tsid };
                        return dtss_server.do_read(id_vec, period, use_ts_cached_read, update_ts_cache);
                    }
                )
            };
        } else {
            throw std::runtime_error{ std::string{"Cannot construct unknown container type: "} + container_type };
        }

        if ( dtss_server.container[container_name]._container.valueless_by_exception() ) {
            throw std::runtime_error{"Container was constructed in an invalid state"};
        }
    }

    static container_wrapper_t & get_container(
        const std::string & container_name, const std::string & container_query,
        dtss::server<standard_dtss_dispatcher> & dtss_server
    ) {
        decltype(dtss_server.container)::iterator f;
        if ( container_query.empty() || container_query == "ts_db" ) {
            f = dtss_server.container.find(container_name);
        } else if ( container_query == "krls" ) {
            f = dtss_server.container.find(std::string{"KRLS_"} + container_name);
        }

        if( f == std::end(dtss_server.container) )
            throw std::runtime_error(std::string("Failed to find shyft container: ")+container_name);

        return f->second;
    }
};


// --------------------------------------------------------------------------------


namespace detail {  // namespace for implementation details not considered public api

// formerly this was only implemented and used in the dtss.cpp translation unit

struct utcperiod_hasher {
    std::size_t operator()(const utcperiod&k) const {
        using boost::hash_value;
        using boost::hash_combine;
        std::size_t seed=0;
        hash_combine(seed,hash_value(k.start.count()));
        hash_combine(seed,hash_value(k.end.count()));
        return seed;
    }
};

}

// NOTE: The following methods may be implemented in-class, remove inline in that case.

template < class CD >
inline ts_info_vector_t server<CD>::do_find_ts(const std::string& search_expression) {
    // 1. filter shyft://<container>/
    auto pattern = extract_shyft_url_container(search_expression);
    if ( pattern.size() > 0 ) {
        // assume it is a shyft url -> look for query flags
        auto queries = extract_shyft_url_query_parameters(search_expression);
        auto container_query_it = queries.find(CD::container_query);
        if ( ! queries.empty() && container_query_it != queries.end() ) {
            auto container_query = container_query_it->second;
            filter_shyft_url_parsed_queries(queries, CD::remove_queries);
            return internal(pattern, container_query).find(extract_shyft_url_path(search_expression), queries);
        } else {
            filter_shyft_url_parsed_queries(queries, CD::remove_queries);
            return internal(pattern).find(extract_shyft_url_path(search_expression), queries);
        }
    } else if (find_ts_cb) {
        return find_ts_cb(search_expression);
    } else {
        return ts_info_vector_t();
    }
}


template < class CD >
inline ts_info server<CD>::do_get_ts_info(const std::string & ts_name) {
    // 1. filter shyft://<container>/
    auto pattern = extract_shyft_url_container(ts_name);
    if ( pattern.size() > 0 ) {
        // assume it is a shyft url -> look for query flags
        auto queries = extract_shyft_url_query_parameters(ts_name);
        auto container_query_it = queries.find(CD::container_query);
        if ( ! queries.empty() && container_query_it != queries.end() ) {
            auto container_query = container_query_it->second;
             filter_shyft_url_parsed_queries(queries, CD::remove_queries);
            return internal(pattern, container_query).get_ts_info(extract_shyft_url_path(ts_name), queries);
        } else {
            filter_shyft_url_parsed_queries(queries, CD::remove_queries);
            return internal(pattern).get_ts_info(extract_shyft_url_path(ts_name), queries);
        }
    } else {
        return ts_info{};
    }
}


/** if overwrite on write, then flush the cache prior to writing */
template < class CD >
inline void server<CD>::do_cache_update_on_write(const ts_vector_t&tsv,bool overwrite_on_write) {
    vector<string> ts_ids;ts_ids.reserve(tsv.size()); 
    ts_vector_t tss;
    for (std::size_t i = 0; i < tsv.size(); ++i) {
        auto rts = dynamic_pointer_cast<aref_ts>(tsv[i].ts);
        ts_ids.emplace_back(rts->id);
        tss.emplace_back(apoint_ts(rts->rep));
    }
    ts_cache.add(ts_ids,tss,overwrite_on_write);
}


template < class CD >
inline void server<CD>::do_store_ts(const ts_vector_t & tsv, bool overwrite_on_write, bool cache_on_write) {
    if(tsv.size()==0) return;
    // 1. filter out all shyft://<container>/<ts-path> elements
    //    and route these to the internal storage controller (threaded)
    //    std::map<std::string, ts_db> shyft_internal;
    //
    std::vector<std::size_t> other;
    other.reserve(tsv.size());
    for( std::size_t i = 0; i < tsv.size(); ++i ) {
        auto rts = dynamic_pointer_cast<aref_ts>(tsv[i].ts);
        if ( ! rts )
            throw std::runtime_error("dtss store: require ts with url-references");

        auto c = extract_shyft_url_container(rts->id);
        if( c.size() > 0 ) {
            auto queries = extract_shyft_url_query_parameters(rts->id);
            auto container_query_it = queries.find(CD::container_query);
            if ( ! queries.empty() && container_query_it != queries.end() ) {
                auto container_query = container_query_it->second;
                filter_shyft_url_parsed_queries(queries, CD::remove_queries);
                internal(c, container_query).save(
                    extract_shyft_url_path(rts->id),  // path
                    rts->core_ts(),      // ts to save
                    overwrite_on_write,  // should do overwrite instead of merge
                    queries              // query key/values from url
                );
            } else {
                filter_shyft_url_parsed_queries(queries, CD::remove_queries);
                internal(c).save(
                    extract_shyft_url_path(rts->id),  // path
                    rts->core_ts(),      // ts to save
                    overwrite_on_write,  // should do overwrite instead of merge
                    queries              // query key/values from url
                );
            }
            if ( cache_on_write ) { // ok, this ends up in a copy, and lock for each item(can be optimized if many)
                if(overwrite_on_write)
                    ts_cache.remove(rts->id);//invalidate previous defs. if any
                ts_cache.add(rts->id, apoint_ts(rts->rep));
            }
        } else {
            other.push_back(i); // keep idx of those we have not saved
        }
    }

    // 2. for all non shyft:// forward those to the
    //    store_ts_cb
    if(store_ts_cb && other.size()) {
        if(other.size()==tsv.size()) { //avoid copy/move if possible
            store_ts_cb(tsv);
            if (cache_on_write) do_cache_update_on_write(tsv,overwrite_on_write);
        } else { // have to do a copy to new vector
            ts_vector_t r;
            for(auto i:other) r.push_back(tsv[i]);
            store_ts_cb(r);
            if (cache_on_write) do_cache_update_on_write(r,overwrite_on_write);
        }
    }
}


template < class CD >
inline void server<CD>::do_merge_store_ts(const ts_vector_t& tsv, bool cache_on_write) {
    if ( tsv.size() == 0 )
        return;

    //
    // 0. check & prepare the read time-series in tsv for the specified period of each ts
    //    (we optimize a little bit grouping on common period, and reading in batches with equal periods)
    //
    id_vector_t ts_ids; ts_ids.reserve(tsv.size());
    std::unordered_map<utcperiod, id_vector_t, detail::utcperiod_hasher> read_map;
    std::unordered_map<std::string, apoint_ts> id_map;

    for ( std::size_t i=0; i < tsv.size(); ++i ) {
        auto rts = dynamic_pointer_cast<aref_ts>(tsv[i].ts);
        if ( ! rts )
            throw std::runtime_error("dtss store merge: require ts with url-references");
        // sanity check
        if ( id_map.find(rts->id) != end(id_map) )
            throw std::runtime_error("dtss store merge requires distinct set of ids, first duplicate found:" + rts->id);
        id_map[rts->id] = apoint_ts(rts->rep);
        // then just build up map[period] = list of time-series to read
        auto rp = rts->rep->total_period();
        if ( read_map.find(rp) != end(read_map) ) {
            read_map[rp].push_back(rts->id);
        } else {
            read_map[rp] = id_vector_t{ rts->id };
        }
    }

    //
    // 1. do the read-merge for each common period, append to final minimal write list
    //
    ts_vector_t tsv_store; tsv_store.reserve(tsv.size());
    for ( auto rr = read_map.begin(); rr != read_map.end(); ++rr ) {
        auto read_ts = do_read(rr->second, rr->first, cache_on_write, cache_on_write);
        // read_ts is in the order of the ts-id-list rr->second
        for ( std::size_t i = 0; i < read_ts.size(); ++i ) {
            auto ts_id = rr->second[i];
            read_ts[i].merge_points(id_map[ts_id]);
            tsv_store.push_back(apoint_ts(ts_id, read_ts[i]));
        }
    }

    //
    // 2. finally write the merged result back to whatever store is there
    //
    do_store_ts(tsv_store, false, cache_on_write);

}


template < class CD >
inline ts_vector_t server<CD>::do_read(const id_vector_t & ts_ids, utcperiod p, bool use_ts_cached_read, bool update_ts_cache) {
    if( ts_ids.size() == 0 )
        return ts_vector_t{};

    // should cache?
    bool cache_read_results = update_ts_cache || cache_all_reads;

    // 0. filter out ts's we can get from cache, given we are allowed to use cache
    std::unordered_map<std::string, apoint_ts> cc;  // cached series
    if( use_ts_cached_read )
        cc = ts_cache.get(ts_ids, p);

    ts_vector_t results(ts_ids.size());
    std::vector<std::size_t> external_idxs;
    if ( cc.size() == ts_ids.size() ) {
        // if we got all ts's from cache -> map in the results
        for( std::size_t i = 0; i < ts_ids.size(); ++i )
            results[i] = cc[ts_ids[i]];
    } else {
        // 1. filter out shyft://
        //    if all shyft: return internal read
        external_idxs.reserve(ts_ids.size()); // only reserve space when needed
        for ( std::size_t i = 0; i < ts_ids.size(); ++i ) {
            if ( cc.find(ts_ids[i]) == cc.end() ) {
                auto c = extract_shyft_url_container(ts_ids[i]);
                if ( c.size() > 0 ) {
                    // check for queries in shyft:// url's
                    auto queries = extract_shyft_url_query_parameters(ts_ids[i]);
                    auto container_query_it = queries.find(CD::container_query);
                    if ( ! queries.empty() && container_query_it != queries.end() ) {
                        auto container_query = container_query_it->second;
                        filter_shyft_url_parsed_queries(queries, CD::remove_queries);
                        results[i] = apoint_ts(make_shared<gpoint_ts>(internal(c, container_query).read(
                            extract_shyft_url_path(ts_ids[i]), p, queries)));
                    } else {
                        filter_shyft_url_parsed_queries(queries, CD::remove_queries);
                        results[i] = apoint_ts(make_shared<gpoint_ts>(internal(c).read(
                            extract_shyft_url_path(ts_ids[i]), p, queries)));
                    }
                    // caching?
                    if ( cache_read_results )
                        ts_cache.add(ts_ids[i], results[i]);
                } else
                    external_idxs.push_back(i);
            } else {
                results[i] = cc[ts_ids[i]];
            }
        }
    }

    // 2. if other/more than shyft get all those
    if( external_idxs.size() > 0 ) {
        if( ! bind_ts_cb )
            throw std::runtime_error("dtss: read-request to external ts, without external handler");

        // only externaly handled series => return only external result
        if( external_idxs.size() == ts_ids.size() ) {
            auto rts = bind_ts_cb(ts_ids,p);
            if( cache_read_results )
                ts_cache.add(ts_ids,rts);

            return rts;
        }

        // collect & handle external references
        vector<string> external_ts_ids; external_ts_ids.reserve(external_idxs.size());
        for( auto i : external_idxs )
            external_ts_ids.push_back(ts_ids[i]);
        auto ext_resolved = bind_ts_cb(external_ts_ids, p);

        // caching?
        if( cache_read_results )
            ts_cache.add(external_ts_ids, ext_resolved);

        // merge external results into output results
        for(std::size_t i=0;i<ext_resolved.size();++i)
            results[external_idxs[i]] = ext_resolved[i];
    }

    return results;
}


template < class CD >
inline void server<CD>::do_bind_ts(utcperiod bind_period, ts_vector_t& atsv, bool use_ts_cached_read, bool update_ts_cache) {

    using shyft::time_series::dd::ts_bind_info;

    std::unordered_map<std::string, std::vector<ts_bind_info>> ts_bind_map;
    std::vector<std::string> ts_id_list;

    // step 1: bind not yet bound time-series (ts with only symbol, needs to be resolved using bind_cb)
    for ( auto & ats : atsv ) {
        auto ts_refs = ats.find_ts_bind_info();
        for ( const auto & bi : ts_refs ) {
            if ( ts_bind_map.find(bi.reference) == ts_bind_map.end() ) { // maintain unique set
                ts_id_list.push_back(bi.reference);
                ts_bind_map[bi.reference] = std::vector<ts_bind_info>();
            }
            ts_bind_map[bi.reference].push_back(bi);
        }
    }

    // step 2: (optional) bind_ts callback should resolve symbol time-series with content
    if ( ts_bind_map.size() > 0 ) {
        auto bts = do_read(ts_id_list, bind_period, use_ts_cached_read, update_ts_cache);
        if ( bts.size() != ts_id_list.size() )
            throw std::runtime_error(std::string{"failed to bind all of "} + std::to_string(bts.size()) + std::string{" ts"});

        for ( std::size_t i = 0; i < ts_id_list.size(); ++i ) {
            for ( auto & bi : ts_bind_map[ts_id_list[i]] )
                bi.ts.bind(bts[i]);
        }
    }

    // step 3: after the symbolic ts are read and bound, we iterate over the
    //         expression tree and calls .do_bind() so that
    //         the new information is taken into account and the expression tree are
    //         ready for evaluate with everything const so threading is safe.
    for ( auto & ats : atsv )
        ats.do_bind();
}


template < class CD >
inline ts_vector_t server<CD>::do_evaluate_ts_vector(utcperiod bind_period, ts_vector_t& atsv,bool use_ts_cached_read,bool update_ts_cache,utcperiod clip_period) {
    do_bind_ts(bind_period, atsv, use_ts_cached_read, update_ts_cache);
    if(clip_period.valid()) 
        return clip_to_period(ts_vector_t{ shyft::time_series::dd::deflate_ts_vector<apoint_ts>(atsv) },clip_period);
    else
        return ts_vector_t{ shyft::time_series::dd::deflate_ts_vector<apoint_ts>(atsv) };
}


template < class CD >
inline ts_vector_t server<CD>::do_evaluate_percentiles(utcperiod bind_period, ts_vector_t& atsv, gta_t const&ta, std::vector<int64_t> const& percentile_spec,bool use_ts_cached_read,bool update_ts_cache) {
    do_bind_ts(bind_period, atsv,use_ts_cached_read,update_ts_cache);
    std::vector<int64_t> p_spec;for(const auto p:percentile_spec) p_spec.push_back(int(p));// convert
    return percentiles(atsv, ta, p_spec);// we can assume the result is trivial to serialize
}


template < class CD >
inline void server<CD>::do_remove_ts(const std::string & ts_url) {
    if ( ! can_remove ) {
        throw std::runtime_error("dtss::server: server does not support removing");
    }

    // 1. filter shyft://<container>/
    auto pattern = extract_shyft_url_container(ts_url);
    if ( pattern.size() > 0 ) {
        // assume it is a shyft url -> look for query flags
        auto queries = extract_shyft_url_query_parameters(ts_url);
        auto container_query_it = queries.find(CD::container_query);
        auto shyft_ts_url = extract_shyft_url_path(ts_url);
        if ( ! queries.empty() && container_query_it != queries.end() ) {
            auto container_query = container_query_it->second;
            filter_shyft_url_parsed_queries(queries, CD::remove_queries);
            internal(pattern, container_query).remove(shyft_ts_url, queries);
        } else {
            filter_shyft_url_parsed_queries(queries, CD::remove_queries);
            internal(pattern).remove(shyft_ts_url, queries);
        }
        ts_cache.remove(shyft_ts_url);// remove it from cache as well!
    } else {
        throw std::runtime_error("dtss::server: server does not allow removing for non shyft-url type data");
    }
}


template < class CD >
inline void server<CD>::on_connect(
    std::istream & in,
    std::ostream & out,
    const std::string & foreign_ip,
    const std::string & local_ip,
    unsigned short foreign_port,
    unsigned short local_port,
    dlib::uint64 connection_id
) {

    using shyft::core::core_iarchive;
    using shyft::core::core_oarchive;
    using shyft::time_series::dd::expression_decompressor;
    using shyft::time_series::dd::compressed_ts_expression;

    while (in.peek() != EOF) {
        auto msg_type= msg::read_type(in);
        try { // scoping the binary-archive could be ok, since it forces destruction time (considerable) to taken immediately, reduce memory foot-print early
              //  at the cost of early& fast response. I leave the commented scopes in there for now, and aim for fastest response-time
            switch (msg_type) { // currently switch, later maybe table[msg_type]=msg_handler
            case message_type::EVALUATE_TS_VECTOR:
            case message_type::EVALUATE_TS_VECTOR_CLIP:
            case message_type::EVALUATE_EXPRESSION:
            case message_type::EVALUATE_EXPRESSION_CLIP:{
                utcperiod bind_period;bool use_ts_cached_read,update_ts_cache;utcperiod clip_period;
                ts_vector_t rtsv;
                core_iarchive ia(in,core_arch_flags);
                ia>>bind_period;
                if(msg_type==message_type::EVALUATE_EXPRESSION || msg_type==message_type::EVALUATE_EXPRESSION_CLIP) {
                    compressed_ts_expression c_expr;
                    ia>>c_expr;
                    rtsv=expression_decompressor::decompress(c_expr);
                } else {
                    ia>>rtsv;
                }
                ia>>use_ts_cached_read>>update_ts_cache;
                if(msg_type==message_type::EVALUATE_TS_VECTOR_CLIP || msg_type==message_type::EVALUATE_EXPRESSION_CLIP) {
                    ia>>clip_period;
                }
                auto result=do_evaluate_ts_vector(bind_period, rtsv,use_ts_cached_read,update_ts_cache,clip_period);//first get result
                msg::write_type(message_type::EVALUATE_TS_VECTOR,out);// then send
                core_oarchive oa(out,core_arch_flags);
                oa<<result;

            } break;
            case message_type::EVALUATE_EXPRESSION_PERCENTILES:
            case message_type::EVALUATE_TS_VECTOR_PERCENTILES: {
                utcperiod bind_period;bool use_ts_cached_read,update_ts_cache;
                ts_vector_t rtsv;
                std::vector<int64_t> percentile_spec;
                gta_t ta;
                core_iarchive ia(in,core_arch_flags);
                ia >> bind_period;
                if(msg_type==message_type::EVALUATE_EXPRESSION_PERCENTILES) {
                    compressed_ts_expression c_expr;
                    ia>>c_expr;
                    rtsv=expression_decompressor::decompress(c_expr);
                } else {
                    ia>>rtsv;
                }

                ia>>ta>>percentile_spec>>use_ts_cached_read>>update_ts_cache;

                auto result = do_evaluate_percentiles(bind_period, rtsv,ta,percentile_spec,use_ts_cached_read,update_ts_cache);//{
                msg::write_type(message_type::EVALUATE_TS_VECTOR_PERCENTILES, out);
                core_oarchive oa(out,core_arch_flags);
                oa << result;
            } break;
            case message_type::FIND_TS: {
                std::string search_expression; //{
                search_expression = msg::read_string(in);// >> search_expression;
                auto find_result = do_find_ts(search_expression);
                msg::write_type(message_type::FIND_TS, out);
                core_oarchive oa(out,core_arch_flags);
                oa << find_result;
            } break;
            case message_type::GET_TS_INFO: {
                std::string ts_url;
                ts_url = msg::read_string(in);
                auto result = do_get_ts_info(ts_url);
                msg::write_type(message_type::GET_TS_INFO, out);
                core_oarchive oa(out, core_arch_flags);
                oa << result;
            } break;
            case message_type::STORE_TS: {
                ts_vector_t rtsv;
                bool overwrite_on_write{ true };
                bool cache_on_write{ false };
                core_iarchive ia(in,core_arch_flags);
                ia >> rtsv >> overwrite_on_write >> cache_on_write;
                do_store_ts(rtsv, overwrite_on_write, cache_on_write);
                msg::write_type(message_type::STORE_TS, out);
            } break;
            case message_type::MERGE_STORE_TS: {
                ts_vector_t rtsv;
                bool cache_on_write{ false };
                core_iarchive ia(in,core_arch_flags);
                ia >> rtsv >> cache_on_write;
                do_merge_store_ts(rtsv, cache_on_write);
                msg::write_type(message_type::MERGE_STORE_TS, out);
            } break;
            case message_type::REMOVE_TS: {
                std::string ts_url = msg::read_string(in);
                do_remove_ts(ts_url);
                msg::write_type(message_type::REMOVE_TS, out);
            } break;
            case message_type::CACHE_FLUSH: {
                flush_cache();
                clear_cache_stats();
                msg::write_type(message_type::CACHE_FLUSH,out);
            } break;
            case message_type::CACHE_STATS: {
                auto cs = get_cache_stats();
                msg::write_type(message_type::CACHE_STATS,out);
                core_oarchive oa(out,core_arch_flags);
                oa<<cs;
            } break;
            default:
                throw std::runtime_error(std::string("Server got unknown message type:") + std::to_string((int)msg_type));
            }
        } catch (std::exception const& e) {
            msg::send_exception(e,out);
        }
    }
}

} // shyft::dtss
} // shyft
