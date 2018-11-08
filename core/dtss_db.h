/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#pragma once

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <memory>
#include <utility>
#include <functional>
#include <cstring>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <regex>

#include <fcntl.h>

#include "time_series_dd.h"
#include "dtss_mutex.h"

#include "time_series_info.h"
#include "utctime_utilities.h"

namespace shyft {
namespace dtss {

using std::move;
using std::vector;
using std::map;
using std::shared_ptr;
using std::make_shared;
using std::unique_ptr;
using std::make_unique;
using std::runtime_error;
using std::string;
using std::size_t;
using std::fopen;
using std::fclose;
using std::fwrite;
using std::fread;
using std::fseek;
using std::FILE;

using shyft::time_series::dd::apoint_ts;
using shyft::time_series::dd::gpoint_ts;
using shyft::time_series::dd::gts_t;
using shyft::time_series::dd::aref_ts;
using shyft::core::utctime;
using shyft::core::utcperiod;
using shyft::core::utctimespan;
using shyft::core::no_utctime;
using shyft::core::max_utctime;
using shyft::core::min_utctime;
using shyft::core::calendar;
using shyft::core::deltahours;
using shyft::core::to_seconds64;
using shyft::core::seconds;
using shyft::time_series::ts_point_fx;
using gta_t = shyft::time_axis::generic_dt;
using ta_generic_type=shyft::time_axis::generic_dt::generic_type;
using gts_t = shyft::time_series::point_ts<gta_t>;


/** \brief Storage layer header-record.
 *
 * The storage layer header-record contains enough
 * information to
 *  a) identify the format version (and thus the rest of the layout)
 *  b) provide  minimal quick summary information that could ease search
 *  c) vital features of the ts that is invariant
 *
 * \note This is a packed record, and should portable through small-endian cpu-arch.
 *
 * Format specification:
 * The binary format of the file is then defined as :
 *   <ts.db.file>   -> <header><time-axis><values>
 *       <header>   -> ()'TS1'|'TS2') <point_fx> <ta_type> <n> <data_period>
 *                            note:
 *                               if 'TS1', then all time-values (int64_t) are in seconds,
 *                               if 'TS2,' all time-values(int64_t) are in micro seconds
 *     <point_fx>   -> ts_point_fx:uint8_t
 *      <ta_type>   -> time_axis::generic_dt::generic_type:uint8_t
 *            <n>   -> uint32_t // number points
 *  <data_period>   -> int64_t int64_t // from .. until
 *
 * <time-axis>      ->
 *   if ta_type == fixed_dt:
 *        <start>   -> int64_t
 *        <delta_t> -> int64_t
 *
 *   if ta_type == calendar_dt:
 *        <start>   -> int64_t
 *        <delta_t> -> int64_t
 *        <tz_sz>   -> uint32_t // the size of tz-info string bytes following
 *        <tz_name> -> uint8_t[<tz_sz>] // length given by  tz_sz above
 *
 *   if ta_type == point_dt:
 *        <t_end>   -> int64_t // the end of the last interval, aka t_end
 *        <t>       -> int64_t[<n>] // <n> from the header
 *
 * <values>         -> double[<n>] // <n> from the header
 *
 */
#pragma pack(push,1)
struct ts_db_header {
    char signature[4] = { 'T','S','1','\0' }; ///< signature header with version #
    time_series::ts_point_fx point_fx = time_series::POINT_AVERAGE_VALUE; ///< point_fx feature
    time_axis::generic_dt::generic_type ta_type = time_axis::generic_dt::FIXED; ///< time-axis type
    uint32_t n = 0; ///< number of points in the time-series (time-axis and values)
    utcperiod data_period; ///< [from..until> period range

    ts_db_header() = default;
    ts_db_header(time_series::ts_point_fx point_fx, time_axis::generic_dt::generic_type ta_type, uint32_t n, utcperiod data_period, char v='1')
        : point_fx(point_fx), ta_type(ta_type), n(n), data_period(data_period) {
     signature[2]=v;
    }
    bool is_seconds() const noexcept {return signature[2]=='1';}
};
#pragma pack(pop)


/** native-time io, using direct write/read utctime, that is micro seconds,
  * zero overhead.
  */
struct native_time_io {
        static char version() {return '2';}
        static void write(FILE *fh,const utctime &t) {
            if( fwrite(&t,sizeof(utctime),1,fh)!=1)
                throw runtime_error("dtss_store: failed to write timepoint to disk");
        }
        static void write(FILE *fh,const vector<utctime> &t) {
            if( fwrite(t.data(),sizeof(utctime),t.size(),fh)!=t.size())
                throw runtime_error("dtss_store: failed to write timepoints to disk");
        }
         static void read(FILE *fh,utctime&t ) {
            if( fread(&t,sizeof(utctime),1,fh)!=1)
                throw runtime_error("dtss_store: failed to read timepoint from from disk");
        }
        static  void read(FILE *fh, vector<utctime> &t) {
            if( fread(t.data(),sizeof(utctime),t.size(),fh)!=t.size())
                throw runtime_error("dtss_store: failed to read timepoints  from disk");
        }
};

/** seconds-based -time io, needs conversion/to/from us/seconds
  * computational  overhead.
  */
struct seconds_time_io {
        static char version() {return '1';}
        static void write(FILE *fh,const utctime &t_us) {
            int64_t t=to_seconds64(t_us);
            if( fwrite(&t,sizeof(utctime),1,fh)!=1)
                throw runtime_error("dtss_store: failed to write timepoint to disk");
        }
        static void write(FILE *fh,const vector<utctime> &t_us) {
            vector<int64_t> t;t.reserve(t_us.size());// have to alloc and convert
            for(const auto&tx:t_us) t.emplace_back(to_seconds64(tx));
            if( fwrite(t.data(),sizeof(utctime),t.size(),fh)!=t.size())
                throw runtime_error("dtss_store: failed to write timepoints to disk");
        }
         static void read(FILE *fh,utctime& t_us) {
            int64_t t;
            if( fread(&t,sizeof(utctime),1,fh)!=1)
                throw runtime_error("dtss_store: failed to read timeponit from from disk");
            t_us=seconds(t);
        }
        static  void read(FILE *fh, vector<utctime> &t_us) {
            vector<int64_t>t;t.resize(t_us.size());
            if( fread(t.data(),sizeof(utctime),t.size(),fh)!=t.size())
                throw runtime_error("dtss_store: failed to read timepoints from disk");
            for(size_t i=0;i<t.size();++i) t_us[i]=seconds(t[i]);
        }
};

/** \brief A simple file-io based internal time-series storage for the dtss.
 *
 * Utilizing standard c++ libraries to store time-series
 * to regular files, that resides in directory containers.
 * Features are limited to simple write/replace, read, search and delete.
 *
 * The simple idea is just to store one time-series(fragment) pr. file.
 *
 *
 * Using a simple client side url naming:
 *
 *  shyft://<container>/<container-relative-path>
 *
 * inside the ts_db at the server-side there is a map
 *   <container> -> root_dir
 * and thus the fullname of the ts is
 *   <container>.root_dir/<container-relative-path>
 *
 *  e.g.:
 *              (proto)  (container ) (   path within container     )
 *  client url: 'shyft://measurements/hydmet_station/1/temperature_1'
 *
 *  server side:
 *    ts_db_container['measurements']=ts_db('/srv/shyft/ts_db/measurements')
 *
 *   which then would resolve into the full-path for the stored ts-file:
 *
 *     '/srv/shyft/ts_db/measurements/hydmet_station/1/temperature_1'
 *
 * \note that ts-urls that do not match internal 'shyft' protocol are dispatched
 *      to the external setup callback (if any). This way we support both
 *      internally managed as well as externally mapped ts-db
 *
 */
struct ts_db {

    using queries_t = map<string, string>;

    string root_dir; ///< root_dir points to the top of the container

  private:
      /** helper class needed for win compensating code */
    struct close_write_handle {
        bool win_thread_close = false;
        mutable ts_db * parent = nullptr;

		close_write_handle() noexcept {};// minimum fix for clang ref. https://stackoverflow.com/questions/43819314/default-member-initializer-needed-within-definition-of-enclosing-class-outside
        close_write_handle(bool wtc) noexcept : win_thread_close{ wtc } {};
        close_write_handle(const close_write_handle &) noexcept = default;

        void operator()(FILE * fh) const {
#ifdef _WIN32WORKAROUND
            if (win_thread_close && parent) {
                parent->fclose_me(fh);
            } else {
                fclose(fh); // takes forever in windows, by design
            }
#else
            fclose(fh);
#endif
        }
    };

    map<string, shared_ptr<calendar>> calendars;

    //--section dealing with windows and (postponing slow) closing files
#ifdef _WIN32WORKAROUND
    mutable mutex fclose_mx;
    mutable vector<std::future<void>> fclose_windows;

    void fclose_me(FILE *fh) {
		lock_guard<decltype(fclose_mx)> sl(fclose_mx);
        fclose_windows.emplace_back(std::async(std::launch::async, [fh]() { fclose(fh); }));
    }

    void wait_for_close_fh() const noexcept {
        try {
            lock_guard<decltype(fclose_mx)> scope_lock(fclose_mx);
            for (auto& fc : fclose_windows)
                fc.get();
            fclose_windows.clear();
        } catch (...) {

        }
    }
#else
    void wait_for_close_fh() const noexcept {}
#endif

	file_lock_manager f_mx;

public:
    bool time_format_micro_seconds=true;///< for testing only, force save in old seconds format, to ensure we can read it.
    ts_db() = default;

    /** constructs a ts_db with specified container root */
    explicit ts_db(const string& root_dir) :root_dir(root_dir) {
        if (!fs::is_directory(root_dir)) {
            if (!fs::exists(root_dir)) {
                if (!fs::create_directories(root_dir)) {
                    throw runtime_error(string("ts_db: failed to create root directory :") + root_dir);
                }
            } else {
                throw runtime_error(string("ts_db: designated root directory is not a directory:") + root_dir);
            }
        }
        make_calendar_lookups();
    }

    ~ts_db() {
        wait_for_close_fh();
    }

	// not supported:
	ts_db(const ts_db&)=delete;
	ts_db(ts_db&&)=delete;
	ts_db & operator=(const ts_db&)=delete;
	ts_db & operator=( ts_db&&) =delete;

    /** \brief Save a time-series to a file, *overwrite* any existing file with contents.
     *
     * \note Windows NTFS is incredibly slow on close, so as a workaround using
     *       handles and separate thread for the close-job for now.
     *
     * \param fn  Pathname to save the time-series at.
     * \param ts  Time-series to save.
     * \param win_thread_close  Only meaningfull on the Windows platform.
     *                          Use a deatached background thread to close the file.
     *                          Defaults to true.
     */
    void save(const string& fn, const gts_t& ts, bool overwrite = true, const queries_t & queries = queries_t{}, bool win_thread_close = true) {
        wait_for_close_fh();

        string ffp = make_full_path(fn, true);
        writer_file_lock lck(f_mx,ffp);
        unique_ptr<FILE, close_write_handle> fh;  // zero-initializes deleter
        fh.get_deleter().win_thread_close = win_thread_close;
        fh.get_deleter().parent = const_cast<ts_db*>(this);
        ts_db_header old_header;

        bool do_merge = false;
        if (!overwrite && save_path_exists(fn)) {
            fh.reset(fopen(ffp.c_str(), "r+b"));
            old_header = read_header(fh.get());
            if (ts.total_period().contains(old_header.data_period)) {
                // old data is completly contained in the new => start the file anew
                //  - reopen, as there is no simple way to truncate an open file...
                //fseek(fh.get(), 0, SEEK_SET);
                wait_for_close_fh();
                fh.reset(fopen(ffp.c_str(), "wb"));
            } else {
                do_merge = true;
            }
        } else {
            fh.reset(fopen(ffp.c_str(), "wb"));
        }
        if (!do_merge) {
            if(time_format_micro_seconds)
                write_ts<native_time_io>(fh.get(), ts);
            else
                write_ts<seconds_time_io>(fh.get(), ts);
        } else {
            if(old_header.is_seconds())
                merge_ts<seconds_time_io>(fh.get(), old_header, ts);
            else 
                merge_ts<native_time_io>(fh.get(), old_header, ts);
        }
    }

    /** read a ts from specified file */
    gts_t read(const string& fn, utcperiod p, const queries_t & queries = queries_t{}) {
        wait_for_close_fh();
        string ffp = make_full_path(fn);
        reader_file_lock lck(f_mx, ffp);

        unique_ptr<FILE, decltype(&fclose)> fh{ fopen(ffp.c_str(), "rb"), &fclose };
        if(!fh.get()) {
            throw runtime_error(string("shyft-read time-series internal: Could not open file ")+ffp );
        }
        return read_ts(fh.get(), p);
    }

    /** removes a ts from the container */
    void remove(const string& fn, const queries_t & queries = queries_t{}) {
        wait_for_close_fh();
        auto fp = make_full_path(fn);
        writer_file_lock lck(f_mx, fp);

        for (size_t retry = 0; retry < 10; ++retry) {
            try {
                fs::remove(fp);
                return;
            } catch (...) { // windows usually fails, due to delayed file-close/file-release so we retry 10 x 0.3 seconds
                std::this_thread::sleep_for(std::chrono::duration<int, std::milli>(300));
            }
        }
        throw runtime_error("failed to remove file '" + fp + "' after 10 repeated attempts lasting for 3 seconds");
    }

    /** get minimal ts-information from specified fn */
    ts_info get_ts_info(const string& fn, const queries_t & queries = queries_t{}) {
        wait_for_close_fh();
        auto ffp = make_full_path(fn);
        reader_file_lock lck(f_mx, ffp);

        if ( save_path_exists(fn) ) {
            unique_ptr<FILE, decltype(&fclose)> fh{ fopen(ffp.c_str(), "rb"), &fclose };
            auto h = read_header(fh.get());
            ts_info i;
            i.name = fn;
            i.point_fx = h.point_fx;
            i.modified = utctime(seconds(fs::last_write_time(ffp)));
            i.data_period = h.data_period;
            // consider time-axis type info, dt.. as well
            return i;
        } else {
            throw runtime_error(string{"ts_db: no ts named: "}+fn);
        }
    }

    /** find all ts_info s that matches the specified re match string
     *
     * e.g.: match= 'hydmet_station/.*_id/temperature'
     *    would find all time-series /hydmet_station/xxx_id/temperature
     */
    vector<ts_info> find(const string& match, const queries_t & queries = queries_t{}) {
        wait_for_close_fh();
        fs::path root(root_dir);
        vector<ts_info> r;
        std::regex r_match(match, std::regex_constants::ECMAScript | std::regex_constants::icase);
		for (auto& x : fs::recursive_directory_iterator(root)) {
            if (fs::is_regular(x.path())) {
                string fn = x.path().lexically_relative(root).generic_string(); // x.path() except root-part
                if (std::regex_search(fn, r_match)) {
                    r.push_back(get_ts_info(fn, queries)); // TODO: maybe multi-core this into a job-queue
                }
            } else if (fs::is_directory(x.path())) {
                // TODO: elide recursion into the x, calling
                // if x.path is not part of match
                //   x.no_push();
            }
        }
        return r;
    }

private:
    shared_ptr<calendar> lookup_calendar(const string& tz) const {
        auto it = calendars.find(tz);
        if (it == calendars.end())
            return make_shared<calendar>(tz);
        return it->second;
    }

    void make_calendar_lookups() {
        for (int hour = -11; hour < 12; hour++) {
            auto c = make_shared<calendar>(deltahours(hour));
            calendars[c->tz_info->name()] = c;
        }
    }

    bool save_path_exists(const string & fn) const {
        fs::path fn_path{ fn }, root_path{ root_dir };
        if (fn_path.is_relative()) {
            fn_path = root_path / fn_path;
        } else {
            // questionable: should we allow outside container specs?
            // return false;
        }
        return fs::is_regular_file(fn_path);
    }
    string make_full_path(const string& fn, bool create_paths = false) const {
        fs::path fn_path{ fn }, root_path{ root_dir };
        // determine path type
        if (fn_path.is_relative()) {
            fn_path = root_path / fn_path;
        } else {  // fn_path.is_absolute()
                  // questionable: should we allow outside container specs?
                  //  - if determined to be fully allowed: remove this branch or throw
        }
        // not a directory and create missing path
        if (fs::is_directory(fn_path)) {
            throw runtime_error(fn_path.string() + " is a directory. Should be a file.");
        } else if (!fs::exists(fn_path) && create_paths) {
            fs::path rp = fn_path.parent_path();
            if (rp.compare(root_path) > 0) {  // if fn contains sub-directory, we have to check that it exits
                if (!fs::is_directory(rp)) {
                    fs::create_directories(rp);
                }
            }
        }
        // -----
        return fn_path.string();
    }

    template<class T>
    ts_db_header mk_header(ts_point_fx pfx, ta_generic_type gt, size_t n, utcperiod p) const {
         if(T::version()=='1') {
             p.start /= utctime::period::den;
             p.end  /= utctime::period::den;
        }
         return ts_db_header{ pfx,gt,uint32_t(n),p,T::version()};
    }
    
    template<class T>
    ts_db_header mk_header(const gts_t& ts)  const  {
        return mk_header<T>(ts.point_interpretation(),ts.time_axis().gt,ts.size(),ts.total_period());
    }

    // ----------

    inline void write(FILE * fh, const void* d, size_t sz) const {
        if (fwrite(d, sizeof(char), sz, fh) != sz)
            throw runtime_error("dtss_store: failed to write do disk");
    }
    template<class T>
    void write_header(FILE * fh, const gts_t & ats) const {
        ts_db_header h = mk_header<T>(ats);
        write(fh, static_cast<const void*>(&h), sizeof(h));
    }
    
    template <class T>
    void write_time_axis(FILE * fh, const gta_t& ta) const {
        switch (ta.gt) {
        case time_axis::generic_dt::FIXED: {
            T::write(fh, ta.f.t);
            T::write(fh, ta.f.dt);
        } break;
        case time_axis::generic_dt::CALENDAR: {
            T::write(fh, ta.c.t);
            T::write(fh, ta.c.dt);
            string tz = ta.c.cal->tz_info->name();
            uint32_t sz = tz.size();
            write(fh, static_cast<const void*>(&sz), sizeof(uint32_t));
            write(fh, static_cast<const void*>(tz.c_str()), sz);
        } break;
        case time_axis::generic_dt::POINT: {
            T::write(fh, ta.p.t_end);
            T::write(fh, ta.p.t);
        } break;
        }
    }
    void write_values(FILE * fh, const vector<double>& v) const {
        write(fh, static_cast<const void*>(v.data()), sizeof(double)*v.size());
    }
    template<class T>
    void write_ts(FILE * fh, const gts_t& ats) const {
        write_header<T>(fh, ats);
        write_time_axis<T>(fh, ats.ta);
        write_values(fh, ats.v);
    }

    // ----------

    template<class T>
    void do_merge(FILE * fh, const ts_db_header & old_header, const time_axis::generic_dt & old_ta, const gts_t & new_ts) const {

        // assume the two time-axes have the same type and are aligned

        const auto old_p = old_header.data_period,
            new_p = new_ts.ta.total_period();

        switch (old_header.ta_type) {
        case time_axis::generic_dt::FIXED: {
            utctime t0, tn;          // first/last time point in the merged data
            size_t keep_old_v = 0,    // number of old values to keep
                nan_v_before = 0,  // number of NaN values to insert between new and old
                nan_v_after = 0;   // number of NaN values to insert between old and new
            vector<double> old_tail(0, 0.);
            // determine first time-axis
            if (old_p.start <= new_p.start) {  // old start before new
                t0 = old_p.start;

                // is there a gap between?
                if (old_p.end < new_p.start) {
                    nan_v_after = (new_p.start - old_p.end) / new_ts.ta.f.dt;
                }

                // determine number of old values to keep
                keep_old_v = (new_p.start - old_p.start) / new_ts.ta.f.dt - nan_v_after;
            } else {  // new start before old
                t0 = new_p.start;

                // is there a gap between?
                if (new_p.end < old_p.start) {
                    nan_v_before = (old_p.start - new_p.end) / new_ts.ta.f.dt;
                }

                // is there is a tail of old values
                if (new_p.end < old_p.end) {
                    size_t count = (old_p.end - new_p.end) / new_ts.ta.f.dt - nan_v_before;
                    old_tail.resize(count);
                    fseek(fh,
                        sizeof(ts_db_header)  // header
                        + 2 * sizeof(int64_t)  // + fixed_dt time-axis
                        + (old_header.n - count) * sizeof(double),  // + values except count last
                        SEEK_SET);
                    read(fh, static_cast<void*>(old_tail.data()), count * sizeof(double));
                }
            }

            // determine last time-axis
            if (old_p.end <= new_p.end) {  // old end before new
                tn = new_p.end;
            } else {
                tn = old_p.end;
            }

            // write header
            auto new_header=mk_header<T>(
                new_ts.fx_policy,
                old_header.ta_type, 
                static_cast<uint32_t>((tn - t0) / new_ts.ta.f.dt), 
                utcperiod{ t0, tn }
            );
            // -----
            fseek(fh, 0, SEEK_SET);  // seek to begining
            write(fh, static_cast<const void*>(&new_header), sizeof(ts_db_header));

            // write time-axis
            T::write(fh, t0);

            // write values
            //  - seek past old values to keep
            fseek(fh, sizeof(int64_t) + keep_old_v * sizeof(double), SEEK_CUR);
            //  - if gap after new -> write NaN
            if (nan_v_after > 0) {
                vector<double> tmp(nan_v_after, shyft::nan);
                write(fh, static_cast<const void*>(tmp.data()), nan_v_after * sizeof(double));
            }
            //  - write new values
            write(fh, static_cast<const void*>(new_ts.v.data()), new_ts.v.size() * sizeof(double));
            //  - if gap before old -> write NaN
            if (nan_v_before > 0) {
                vector<double> tmp(nan_v_before, shyft::nan);
                write(fh, static_cast<const void*>(tmp.data()), nan_v_before * sizeof(double));
            }
            //  - if old tail values -> write them back
            if (old_tail.size() > 0) {
                write(fh, static_cast<const void*>(old_tail.data()), old_tail.size() * sizeof(double));
            }
        } break;
        case time_axis::generic_dt::CALENDAR: {
            utctime t0, tn;          // first/last time point in the merged data
            size_t keep_old_v = 0,    // number of old values to keep
                nan_v_before = 0,  // number of NaN values to insert between new and old
                nan_v_after = 0;   // number of NaN values to insert between old and new
            vector<double> old_tail(0, 0.);
            // determine first time-axis
            if (old_p.start <= new_p.start) {  // old start before new
                t0 = old_p.start;

                // is there a gap between?
                if (old_p.end < new_p.start) {
                    nan_v_after = new_ts.ta.c.cal->diff_units(old_p.end, new_p.start, new_ts.ta.c.dt);
                }

                // determine number of old values to keep
                keep_old_v = new_ts.ta.c.cal->diff_units(old_p.start, new_p.start, new_ts.ta.c.dt) - nan_v_after;
            } else {  // new start before old
                t0 = new_p.start;

                // is there a gap between?
                if (new_p.end < old_p.start) {
                    nan_v_before = new_ts.ta.c.cal->diff_units(new_p.end, old_p.start, new_ts.ta.c.dt);
                }

                // is there is a tail of old values
                if (new_p.end < old_p.end) {
                    size_t count = new_ts.ta.c.cal->diff_units(new_p.end, old_p.end, new_ts.ta.c.dt) - nan_v_before;
                    old_tail.resize(count);
                    fseek(fh,
                        sizeof(ts_db_header)  // header
                        + 2 * sizeof(int64_t),  // + first part of calendar_dt time-axis
                        SEEK_SET);
                    uint32_t tz_sz{};
                    read(fh, static_cast<void *>(&tz_sz), sizeof(tz_sz));  // read size of tz_name
                    fseek(fh,
                        tz_sz * sizeof(uint8_t)  // + second part of calendar_dt time-axis
                        + (old_header.n - count) * sizeof(double),  // values to overwrite
                        SEEK_CUR);
                    read(fh, static_cast<void*>(old_tail.data()), count * sizeof(double));
                }
            }

            // determine last time-axis
            if (old_p.end <= new_p.end) {  // old end before new
                tn = new_p.end;
            } else {
                tn = old_p.end;
            }

            // write header
            auto  new_header= mk_header<T>(
                new_ts.fx_policy, old_header.ta_type,
                static_cast<uint32_t>(new_ts.ta.c.cal->diff_units(t0, tn, new_ts.ta.c.dt)),
                utcperiod{ t0, tn }
            );
            // -----
            fseek(fh, 0, SEEK_SET);  // seek to begining
            write(fh, static_cast<const void*>(&new_header), sizeof(ts_db_header));

            // update time-axis
            T::write(fh, t0);
            fseek(fh, sizeof(int64_t), SEEK_CUR);
            {
                uint32_t tz_sz{};
                read(fh, static_cast<void *>(&tz_sz), sizeof(tz_sz));  // read size of calendar str
                fseek(fh, tz_sz * sizeof(uint8_t), SEEK_CUR);  // seek past tz_name
            }

            // write values
            //  - seek past old values to keep
            fseek(fh, keep_old_v * sizeof(double), SEEK_CUR);
            //  - if gap after new -> write NaN
            if (nan_v_after > 0) {
                vector<double> tmp(nan_v_after, shyft::nan);
                write(fh, static_cast<const void*>(tmp.data()), nan_v_after * sizeof(double));
            }
            //  - write new values
            write(fh, static_cast<const void*>(new_ts.v.data()), new_ts.v.size() * sizeof(double));
            //  - if gap before old -> write NaN
            if (nan_v_before > 0) {
                vector<double> tmp(nan_v_before, shyft::nan);
                write(fh, static_cast<const void*>(tmp.data()), nan_v_before * sizeof(double));
            }
            //  - if old tail values -> write them back
            if (old_tail.size() > 0) {
                write(fh, static_cast<const void*>(old_tail.data()), old_tail.size() * sizeof(double));
            }
        } break;
        case time_axis::generic_dt::POINT: {
            vector<utctime> merged_t;
            merged_t.reserve(old_ta.size() + new_ts.size() + 1);  // guaranteed large enough
                                                                  // -----
            vector<double> merged_v(0u, 0.);
            merged_v.reserve(old_ta.size() + new_ts.size() + 1);  // guaranteed large enough
                                                                  // -----
            const size_t old_v_offset = sizeof(ts_db_header) + (old_header.n + 1) * sizeof(int64_t);

            // new start AFTER  =>  start at old
            if (old_p.start < new_p.start) {

                // get iterator into old time-axis at first where p.start >= old.start
                //  - [old_ta.begin(), old_end) is the part of old we want to keep
                auto old_end = std::lower_bound(
                    old_ta.p.t.cbegin(), old_ta.p.t.cend(),
                    new_ts.total_period().start);

                // store portion of old time-points to keep
                merged_t.insert(merged_t.end(), old_ta.p.t.cbegin(), old_end);
                // store portion of old values to keep
                const size_t to_insert = std::distance(old_ta.p.t.cbegin(), old_end);
                auto it = merged_v.insert(merged_v.end(), to_insert, 0.);
                fseek(fh, old_v_offset, SEEK_SET);
                read(fh, static_cast<void*>(&(*it)), to_insert * sizeof(double));

                // if NEW start truly after OLD include the end point
                if (old_end == old_ta.p.t.cend() && new_p.start > old_ta.p.t_end) {
                    merged_t.emplace_back(old_ta.p.t_end);  // include the end point
                    merged_v.emplace_back(shyft::nan);      // nan for the gap
                }
            }

            // read new into merge_ts
            merged_t.insert(merged_t.end(), new_ts.ta.p.t.cbegin(), new_ts.ta.p.t.cend());
            merged_t.emplace_back(new_ts.ta.p.t_end);
            merged_v.insert(merged_v.end(), new_ts.v.cbegin(), new_ts.v.cend());
            // if new end BEFORE start of old  =>  insert NaN
            if (new_p.end < old_p.start) {
                merged_v.emplace_back(shyft::nan);
            }

            // new end BEFORE end of old  =>  read more from old
            if (new_p.end < old_p.end) {

                // determine first period in old NOT CONTAINING new.end
                auto old_begin = std::upper_bound(
                    old_ta.p.t.cbegin(), old_ta.p.t.cend(),
                    new_ts.ta.p.t_end);

                // store any trailing old time-points
                merged_t.insert(merged_t.end(), old_begin, old_ta.p.t.cend());
                merged_t.emplace_back(old_ta.p.t_end);
                // store portion of old values to keep
                size_t to_insert = std::distance(old_begin, old_ta.p.t.cend());
                // new end INSIDE of AT start of old  =>  insert value from old where new end
                if (new_p.end >= old_p.start) {
                    to_insert += 1;
                }
                auto it = merged_v.insert(merged_v.end(), to_insert, 0.);
                fseek(fh, old_v_offset + (old_header.n - to_insert) * sizeof(double), SEEK_SET);
                read(fh, static_cast<void*>(&(*it)), to_insert * sizeof(double));
            }

            // write header
            auto new_header=mk_header<T>(
                new_ts.fx_policy, old_header.ta_type,
                static_cast<uint32_t>(merged_t.size() - 1),
                utcperiod{ merged_t.at(0), merged_t.at(merged_t.size() - 1) }
            );
            // -----
            fseek(fh, 0, SEEK_SET);  // seek to begining
            write(fh, static_cast<const void *>(&new_header), sizeof(ts_db_header));

            // write time-axis
            T::write(fh, merged_t.at(merged_t.size() - 1));
            merged_t.pop_back();
            T::write(fh, merged_t);
            // write values
            write(fh, static_cast<const void *>(merged_v.data()), merged_v.size() * sizeof(double));
        } break;
        }
    }
    void check_ta_alignment(FILE * fh, const ts_db_header & old_header, const time_axis::generic_dt & old_ta, const gts_t & ats) const {
        if (ats.ta.gt != old_header.ta_type) {
            throw runtime_error("dtss_store: cannot merge with different ta type");
        } else {
            // parse specific ta data to determine compatibility
            switch (old_header.ta_type) {
            case time_axis::generic_dt::FIXED:
            {
                // refuse to merge into old second base storage. Old series should be converted
                if ( old_header.is_seconds() && ats.ta.f.t.count() % core::seconds{1}.count() != 0 ) {
                    throw runtime_error("dtss_store: cannot merge microseconds to old seconds based storage ts-file");
                } else if (old_ta.f.dt != ats.ta.f.dt || (old_ta.f.t - ats.ta.f.t).count() % old_ta.f.dt.count() != 0) {
                    throw runtime_error("dtss_store: cannot merge unaligned fixed_dt");
                }
            } break;
            case time_axis::generic_dt::CALENDAR:
            {
                if (ats.ta.c.cal->tz_info->tz.tz_name == old_ta.c.cal->tz_info->tz.tz_name) {
                    utctimespan remainder;
                    ats.ta.c.cal->diff_units(ats.ta.c.t, seconds{0}, seconds{1}, remainder);
                    if ( old_header.is_seconds() && remainder.count() != 0 ) {
                        throw runtime_error("dtss_store: cannot merge microseconds to old seconds based storage ts-file");
                    } else {
                        ats.ta.c.cal->diff_units(old_ta.c.t, ats.ta.c.t, old_ta.c.dt, remainder);
                        if (old_ta.c.dt != ats.ta.c.dt || remainder != utctimespan{0}) {
                            throw runtime_error("dtss_store: cannot merge unaligned calendar_dt");
                        }
                    }
                } else {
                    throw runtime_error("dtss_store: cannot merge calendar_dt with different calendars");
                }
            } break;
            case time_axis::generic_dt::POINT:
                if(old_header.is_seconds()) { // check if us resolution in new time-axis.. throw if ..
                    for(const auto& t:ats.ta.p.t) {
                        if( t!= seconds(to_seconds64(t)))
                            throw runtime_error("dtss_store: can not merge us resolution to old seconds based ts-file");
                    }
                    if(ats.ta.p.t_end!=no_utctime && (ats.ta.p.t_end!=seconds(to_seconds64(ats.ta.p.t_end))))
                        throw runtime_error("dtss_store: can not merge us resolution to old seconds based ts-file");
                }
            break;
            }
        }
    }
    template <class T>
    void merge_ts(FILE * fh, const ts_db_header & old_header, const gts_t & ats) const {
        // read time-axis
        size_t ignored{};
        time_axis::generic_dt old_ta = read_time_axis<T>(fh, old_header, old_header.data_period, ignored);
        check_ta_alignment(fh, old_header, old_ta, ats);
        do_merge<T>(fh, old_header, old_ta, ats);
    }

    // ----------

    inline void read(FILE * fh, void* d, size_t sz) const {
        size_t rsz = fread(d, sizeof(char), sz, fh);
		if (rsz != sz) {
			string fn{"?"};
			throw runtime_error("dtss_store: failed to read '" + fn + "'from disk expected size="+std::to_string(sz)+"!="+std::to_string(rsz));
		}
    }
    ts_db_header read_header(FILE * fh) const {
        ts_db_header h;
        fseek(fh, 0, SEEK_SET);
        read(fh, static_cast<void *>(&h), sizeof(ts_db_header));
        if(h.is_seconds()) {
            h.data_period.start *= utctime::period::den;
            h.data_period.end *= utctime::period::den;
        }
        return h;
    }
    template <class T>
    gta_t read_time_axis(FILE * fh, const ts_db_header& h, const utcperiod p, size_t& skip_n) const {

        // seek to beginning of time-axis
        fseek(fh, sizeof(ts_db_header), SEEK_SET);

        gta_t ta;
        ta.gt = h.ta_type;

        utctime t_start = p.start;
        utctime t_end = p.end;
        if (t_start == no_utctime)
            t_start = min_utctime;
        if (t_end == no_utctime)
            t_end = max_utctime;

        // no overlap?
        if (h.data_period.end <= t_start || h.data_period.start >= t_end) {
            return ta;
        }

        skip_n = 0;
        utctime t0 = no_utctime;
        switch (h.ta_type) {
        case time_axis::generic_dt::FIXED: {
            // read start & step
            T::read(fh,t0);
            T::read(fh, ta.f.dt);
            // handle various overlaping periods
            if (t_start <= h.data_period.start && t_end >= h.data_period.end) {
                // fully around or exact
                ta.f.t = t0;
                ta.f.n = h.n;
            } else {
                size_t drop_n = 0;
                if (t_start > h.data_period.start)  // start inside
                    skip_n = (t_start - h.data_period.start) / ta.f.dt;
                if (t_end < h.data_period.end)  // end inside
                    drop_n = (h.data_period.end - t_end) / ta.f.dt;
                // -----
                ta.f.t = t0 + ta.f.dt*skip_n;
                ta.f.n = h.n - skip_n - drop_n;
            }
        } break;
        case time_axis::generic_dt::CALENDAR: {
            // read start & step
            T::read(fh, t0);
            T::read(fh, ta.c.dt);
            // read tz_info
            uint32_t sz{ 0 };
            read(fh, static_cast<void *>(&sz), sizeof(uint32_t));
            string tz(sz, '\0');
            {
                unique_ptr<char[]> tmp_ptr = make_unique<char[]>(sz);
                read(fh, static_cast<void*>(tmp_ptr.get()), sz);
                tz.replace(0, sz, tmp_ptr.get(), sz);
            }
            ta.c.cal = lookup_calendar(tz);
            // handle various overlaping periods
            if (t_start <= h.data_period.start && t_end >= h.data_period.end) {
                // fully around or exact
                ta.c.t = t0;
                ta.c.n = h.n;
            } else {
                size_t drop_n = 0;
                if (t_start > h.data_period.start)  // start inside
                    skip_n = ta.c.cal->diff_units(h.data_period.start, t_start, ta.c.dt);
                if (t_end < h.data_period.end)  // end inside
                    drop_n = ta.c.cal->diff_units(t_end, h.data_period.end, ta.c.dt);
                // -----
                ta.c.t = ta.c.cal->add(t0, ta.c.dt, skip_n);
                ta.c.n = h.n - skip_n - drop_n;
            }
        } break;
        case time_axis::generic_dt::POINT: {
            if (t_start <= h.data_period.start && t_end >= h.data_period.end) {
                // fully around or exact
                ta.p.t.resize(h.n);
                T::read(fh,ta.p.t_end);
                T::read(fh, ta.p.t);
            } else {
                utctime f_time{seconds{0}};
                vector<utctime> tmp(h.n, seconds{0});
                T::read(fh, f_time);
                T::read(fh, tmp);
                // -----
                auto it_b = tmp.begin();
                if (t_start > h.data_period.start) {
                    it_b = std::upper_bound(tmp.begin(), tmp.end(), t_start, std::less<utctime>());
                    if (it_b != tmp.begin()) {
                        std::advance(it_b, -1);
                    }
                }
                // -----
                auto it_e = tmp.end();
                if (t_end < h.data_period.end) {
                    it_e = std::upper_bound(it_b, tmp.end(), t_end, std::less<utctime>());
                    if (it_e != tmp.end()) {
                        f_time = *it_e;
                    }
                }
                // -----
                skip_n = std::distance(tmp.begin(), it_b);
                ta.p.t.reserve(std::distance(it_b, it_e));
                ta.p.t.assign(it_b, it_e);
                ta.p.t_end = f_time;
            }
        } break;
        }
        return ta;
    }
    
    vector<double> read_values(FILE * fh, const ts_db_header& h, const gta_t& ta, const size_t skip_n) const {

        // seek to beginning of values
        fseek(fh, sizeof(ts_db_header), SEEK_SET);
        switch (h.ta_type) {
        case time_axis::generic_dt::FIXED: {
            fseek(fh, 2 * sizeof(int64_t), SEEK_CUR);
        } break;
        case time_axis::generic_dt::CALENDAR: {
            fseek(fh, 2 * sizeof(int64_t), SEEK_CUR);
            uint32_t sz{};
            read(fh, static_cast<void*>(&sz), sizeof(uint32_t));
            fseek(fh, sz * sizeof(uint8_t), SEEK_CUR);
        } break;
        case time_axis::generic_dt::POINT: {
            fseek(fh, (h.n + 1) * sizeof(int64_t), SEEK_CUR);
        } break;
        }

        const size_t points_n = ta.size();
        vector<double> val(points_n, 0.);
        fseek(fh, sizeof(double)*skip_n, SEEK_CUR);
        read(fh, static_cast<void *>(val.data()), sizeof(double)*points_n);
        return val;
    }

    gts_t read_ts(FILE * fh, const utcperiod p) const {
        size_t skip_n = 0u;
        ts_db_header h = read_header(fh);
        gta_t ta = h.is_seconds()?read_time_axis<seconds_time_io>(fh, h, p, skip_n):read_time_axis<native_time_io>(fh,h,p,skip_n);
        vector<double> v = read_values(fh, h, ta, skip_n);
        return gts_t{ move(ta),move(v),h.point_fx };
    }

};

}
}
