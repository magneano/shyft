/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#include <dlib/statistics.h>
#include <memory>
#include "time_series_dd.h"
#include "time_series_merge.h"
#include "time_series_qm.h"
#include "time_series_point_merge.h"

namespace shyft{
    namespace time_series {
        /*
         Some comments:
         The api-type ts, time-axis gives some challenges when it comes to performance on large amount of data.
         The computation tends to be memory bound, - and the .values() approach kind of does not help in that matter.
         On the other side, a zillion virtual-dispatches down to the terminal-values (when evaluating expressions)
         also have higher cost (tested), although one could argue that in an ideal case, it should be more core-bound
         doing just the math on fewer memory-refs.

         As a common approach here, we say that if the set of time-series do have common time-axis (same-resolution/range),
         we try to perform the operations as fast as vector vs vector operations.

         Tests comparing with armadillo, indicates that we succeeds with this approach.

         The system is still fast, when time-axis are not aligned, but it takes the penalty of aligning the accesses (which is needed!).
        */

        /** ts_src template helps using the average_value/accumulate value template functions
        *
        *  During the evaluation of apoint_ts, .average /.integral functions, the ts at hand
        *  can be anything from a complex expression, up to a terminal-value (a concrete point ts).
        *  To avoid making a copy of terminal-value ts, we check if the argument is a terminal-value type of ts
        *  and instead of copy, just reference the terminal-ts values (if the time-axis are aligned).
        *  This saves us a copy of the terminal-value ts (which is significant part of computation),
        */
        template<class TA>
        struct ts_src {
            const TA& ta;
            const std::vector<double>& v;
            ts_src(const TA& ta,const std::vector<double>&v):ta{ta},v{v} {}
            point get(size_t i ) const {return point{ta.time(i),v[i]};}
            size_t size() const {return v.size();}
        };
        /* specialization of the hint-based searches,
         * SiH: really questionable if we should use open_range_index here.
         */
        template<>
        inline size_t hint_based_search<ts_src<time_axis::fixed_dt>>(const ts_src<time_axis::fixed_dt>& source, const utcperiod& p, size_t i) {
            return source.ta.index_of(p.start);
        }
        template<>
        inline size_t hint_based_search<ts_src<time_axis::calendar_dt>>(const ts_src<time_axis::calendar_dt>& source, const utcperiod& p, size_t i) {
            return source.ta.index_of(p.start);
        }
        template<>
        inline size_t hint_based_search<ts_src<time_axis::point_dt>>(const ts_src<time_axis::point_dt>& source, const utcperiod& p, size_t i) {
            return source.ta.index_of(p.start,i);
        }
        template<>
        inline size_t hint_based_search<ts_src<time_axis::generic_dt>>(const ts_src<time_axis::generic_dt>& source, const utcperiod& p, size_t i) {
            return source.ta.index_of(p.start,i);
        }

        namespace dd {

        static inline double do_op(double a, iop_t op, double b) {
            switch (op) {
            case iop_t::OP_ADD:return a + b;
            case iop_t::OP_SUB:return a - b;
            case iop_t::OP_DIV:return a / b;
            case iop_t::OP_MUL:return a * b;
            case iop_t::OP_MAX:return std::max(a, b);
            case iop_t::OP_MIN:return std::min(a, b);
            case iop_t::OP_POW:return std::pow(a, b);
            case iop_t::OP_NONE:break;// just fall to exception
            }
            throw std::runtime_error("unsupported shyft::api::iop_t");
        }
        // add operators and functions to the apoint_ts class, of all variants that we want to expose
        apoint_ts average(const apoint_ts& ts, const gta_t& ta/*fx-type */) { return apoint_ts(std::make_shared<average_ts>(ta, ts)); }
        apoint_ts average(apoint_ts&& ts, const gta_t& ta) { return apoint_ts(std::make_shared<average_ts>(ta, std::move(ts))); }
        apoint_ts integral(const apoint_ts& ts, const gta_t& ta/*fx-type */) { return apoint_ts(std::make_shared<integral_ts>(ta, ts)); }
        apoint_ts integral(apoint_ts&& ts, const gta_t& ta) { return apoint_ts(std::make_shared<integral_ts>(ta, std::move(ts))); }

        apoint_ts accumulate(const apoint_ts& ts, const gta_t& ta/*fx-type */) { return apoint_ts(std::make_shared<accumulate_ts>(ta, ts)); }
        apoint_ts accumulate(apoint_ts&& ts, const gta_t& ta) { return apoint_ts(std::make_shared<accumulate_ts>(ta, std::move(ts))); }

        apoint_ts create_periodic_pattern_ts(const vector<double>& pattern, utctimespan dt, utctime pattern_t0, const gta_t& ta) { return apoint_ts(make_shared<periodic_ts>(pattern, dt, pattern_t0, ta)); }

        apoint_ts operator+(const apoint_ts& lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_ts       >(lhs, iop_t::OP_ADD, rhs)); }
        apoint_ts operator+(const apoint_ts& lhs, double           rhs) { return apoint_ts(std::make_shared<abin_op_ts_scalar>(lhs, iop_t::OP_ADD, rhs)); }
        apoint_ts operator+(double           lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_scalar_ts>(lhs, iop_t::OP_ADD, rhs)); }

        apoint_ts operator-(const apoint_ts& lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_ts       >(lhs, iop_t::OP_SUB, rhs)); }
        apoint_ts operator-(const apoint_ts& lhs, double           rhs) { return apoint_ts(std::make_shared<abin_op_ts_scalar>(lhs, iop_t::OP_SUB, rhs)); }
        apoint_ts operator-(double           lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_scalar_ts>(lhs, iop_t::OP_SUB, rhs)); }
        apoint_ts operator-(const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_scalar_ts>(-1.0, iop_t::OP_MUL, rhs)); }

        apoint_ts operator/(const apoint_ts& lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_ts       >(lhs, iop_t::OP_DIV, rhs)); }
        apoint_ts operator/(const apoint_ts& lhs, double           rhs) { return apoint_ts(std::make_shared<abin_op_ts_scalar>(lhs, iop_t::OP_DIV, rhs)); }
        apoint_ts operator/(double           lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_scalar_ts>(lhs, iop_t::OP_DIV, rhs)); }

        apoint_ts operator*(const apoint_ts& lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_ts       >(lhs, iop_t::OP_MUL, rhs)); }
        apoint_ts operator*(const apoint_ts& lhs, double           rhs) { return apoint_ts(std::make_shared<abin_op_ts_scalar>(lhs, iop_t::OP_MUL, rhs)); }
        apoint_ts operator*(double           lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_scalar_ts>(lhs, iop_t::OP_MUL, rhs)); }


        apoint_ts max(const apoint_ts& lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_ts       >(lhs, iop_t::OP_MAX, rhs)); }
        apoint_ts max(const apoint_ts& lhs, double           rhs) { return apoint_ts(std::make_shared<abin_op_ts_scalar>(lhs, iop_t::OP_MAX, rhs)); }
        apoint_ts max(double           lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_scalar_ts>(lhs, iop_t::OP_MAX, rhs)); }

        apoint_ts min(const apoint_ts& lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_ts>(lhs, iop_t::OP_MIN, rhs)); }
        apoint_ts min(const apoint_ts& lhs, double           rhs) { return apoint_ts(std::make_shared<abin_op_ts_scalar>(lhs, iop_t::OP_MIN, rhs)); }
        apoint_ts min(double           lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_scalar_ts>(lhs, iop_t::OP_MIN, rhs)); }

        apoint_ts pow(const apoint_ts& lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_ts>(lhs, iop_t::OP_POW, rhs)); }
        apoint_ts pow(const apoint_ts& lhs, double           rhs) { return apoint_ts(std::make_shared<abin_op_ts_scalar>(lhs, iop_t::OP_POW, rhs)); }
        apoint_ts pow(double           lhs, const apoint_ts& rhs) { return apoint_ts(std::make_shared<abin_op_scalar_ts>(lhs, iop_t::OP_POW, rhs)); }

        
        double abin_op_ts::value_at(utctime t) const {
            if (!time_axis().total_period().contains(t))
                return nan;
            return do_op(lhs(t), op, rhs(t));// this might cost a 2xbin-search if not the underlying ts have smart incremental search (at the cost of thread safety)
        }

        static vector<double> ts_op_ts_values(const vector<double>& l, iop_t op, const vector<double>& r) {
            vector<double> x; x.reserve(l.size());
            switch (op) {
            case OP_ADD:for (size_t i = 0; i < r.size(); ++i) x.emplace_back(l[i] + r[i]); return x;
            case OP_SUB:for (size_t i = 0; i < r.size(); ++i) x.emplace_back(l[i] - r[i]); return x;
            case OP_MUL:for (size_t i = 0; i < r.size(); ++i) x.emplace_back(l[i] * r[i]); return x;
            case OP_DIV:for (size_t i = 0; i < r.size(); ++i) x.emplace_back(l[i] / r[i]); return x;
            case OP_MAX:for (size_t i = 0; i < r.size(); ++i) x.emplace_back(std::max(l[i], r[i])); return x;
            case OP_MIN:for (size_t i = 0; i < r.size(); ++i) x.emplace_back(std::min(l[i], r[i])); return x;
            case OP_POW:for (size_t i = 0; i < r.size(); ++i) x.emplace_back(std::pow(l[i], r[i])); return x;
            default: break;
            }
            throw runtime_error("Unsupported operation " + to_string(int(op)));
        }

        static void lhs_in_place_ts_op_ts_values(vector<double>& l, iop_t op, const vector<double>& r) {
            switch (op) {
            case OP_ADD:for (size_t i = 0; i < r.size(); ++i) l[i] += r[i]; return;
            case OP_SUB:for (size_t i = 0; i < r.size(); ++i) l[i] -= r[i]; return;
            case OP_MUL:for (size_t i = 0; i < r.size(); ++i) l[i] *= r[i]; return;
            case OP_DIV:for (size_t i = 0; i < r.size(); ++i) l[i] /= r[i]; return;
            case OP_MAX:for (size_t i = 0; i < r.size(); ++i) l[i] = std::max(l[i], r[i]); return;
            case OP_MIN:for (size_t i = 0; i < r.size(); ++i) l[i] = std::min(l[i], r[i]); return;
            case OP_POW:for (size_t i = 0; i < r.size(); ++i) l[i] = std::pow(l[i], r[i]); return;
            default: break;
            }
            throw runtime_error("Unsupported operation " + to_string(int(op)));
        }
        static void rhs_in_place_ts_op_ts_values(const vector<double>& l, iop_t op, vector<double>& r) {
            switch (op) {
            case OP_ADD:for (size_t i = 0; i < r.size(); ++i) r[i] += l[i]; return;
            case OP_SUB:for (size_t i = 0; i < r.size(); ++i) r[i] = l[i] - r[i]; return;
            case OP_MUL:for (size_t i = 0; i < r.size(); ++i) r[i] *= l[i]; return;
            case OP_DIV:for (size_t i = 0; i < r.size(); ++i) r[i] = l[i] / r[i]; return;
            case OP_MAX:for (size_t i = 0; i < r.size(); ++i) r[i] = std::max(l[i], r[i]); return;
            case OP_MIN:for (size_t i = 0; i < r.size(); ++i) r[i] = std::min(l[i], r[i]); return;
            case OP_POW:for (size_t i = 0; i < r.size(); ++i) r[i] = std::pow(l[i], r[i]); return;
            default: break;
            }
            throw runtime_error("Unsupported operation " + to_string(int(op)));
        }

        static inline const vector<double>* terminal_values(const shared_ptr<ipoint_ts>& ts) {
            if (dynamic_pointer_cast<aref_ts>(ts))
                return &dynamic_pointer_cast<aref_ts>(ts)->core_ts().v;
            if (dynamic_pointer_cast<gpoint_ts>(ts))
                return &dynamic_pointer_cast<gpoint_ts>(ts)->core_ts().v;
            return nullptr;
        }

        static inline const vector<double>* terminal_values(const apoint_ts&ats) {
            return terminal_values(ats.ts);

        }

        /** implementation of apoint_ts op apoint_ts
        *
        *  Try to optimize the evaluation by considering
        *  equal time-axis.
        *  Then also check  rhs and lhs.
        *  If one is a terminal-value (concrete point ts)
        *  then just reference the values, otherwise get the computed values.
        *
        * If time-axis not aligned, just compute value-by-value.
        *
        */
        std::vector<double> abin_op_ts::values() const {
            if (lhs.time_axis() == rhs.time_axis()) {
                const vector<double>* lhs_v{ terminal_values(lhs) };
                const vector<double>* rhs_v{ terminal_values(rhs) };
                if (lhs_v && rhs_v) {
                    return ts_op_ts_values(*lhs_v, op, *rhs_v);
                } else if (lhs_v) {
                    auto r{ rhs.values() };
                    rhs_in_place_ts_op_ts_values(*lhs_v, op, r);
                    return r;
                } else if (rhs_v) {
                    auto l{ lhs.values() };
                    lhs_in_place_ts_op_ts_values(l, op, *rhs_v);
                    return l;
                } else {
                    auto l{ lhs.values() };
                    auto r{ rhs.values() };
                    lhs_in_place_ts_op_ts_values(l, op, r);
                    return l;
                }
            } else {
                std::vector<double> r; r.reserve(time_axis().size());
                for (size_t i = 0; i < time_axis().size(); ++i) {
                    r.push_back(value(i));//TODO: improve speed using accessors with ix-hint for lhs/rhs stepwise traversal
                }
                return r;
            }
        }


        /** make_interval template
        *
        * Help creating implementation of interval type of results (average/integral value)
        * at highest possible speed
        */
        template < class FX, class TA, class TS, class ATA >
        static std::vector<double> make_interval(FX&&fx, const TA &ta, const TS& ts, const ATA& avg_ta) {
            auto pfx = ts->point_interpretation();
            const bool linear_interpretation = pfx == ts_point_fx::POINT_INSTANT_VALUE;
            const vector<double>* src_v{ terminal_values(ts) };//fetch terminal value
            if (src_v) {// if values are already in place, use reference
                const ts_src<TA> rts{ ta,*src_v };
                std::vector<double> r; r.reserve(avg_ta.size());
                size_t ix_hint = ta.index_of(avg_ta.time(0));
                for (size_t i = 0; i < avg_ta.size(); ++i)
                    r.emplace_back(fx(rts, avg_ta.period(i), ix_hint, linear_interpretation));
                return r;
            } else {
                auto tsv{ ts->values() };// pull out values from underlying expression.
                const ts_src<TA> rts{ ta,tsv };// note that ->values() deflates any underlying expression, this could be sub optimal in some cases.
                std::vector<double> r; r.reserve(avg_ta.size());
                size_t ix_hint = ta.index_of(avg_ta.time(0));
                for (size_t i = 0; i < avg_ta.size(); ++i)
                    r.emplace_back(fx(rts, avg_ta.period(i), ix_hint, linear_interpretation));
                return r;
            }
        }

        /** Implementation of the average_ts::values
        *
        *  Important here is to ensure that accessing the time-axis is
        *  direct into the core-type implementation of the generic_dt.
        *  (avoid switch every single lookup during integration)
        */
        std::vector<double> average_ts::values() const {
            if (ts->time_axis() == ta && ts->point_interpretation() == ts_point_fx::POINT_AVERAGE_VALUE) {
                return ts->values(); // elide trivial average_ts cases
            }
            switch (ts->time_axis().gt) { // pull out the real&fast time-axis before doing computations here:
            case time_axis::generic_dt::FIXED:    return make_interval(average_value<ts_src<time_axis::fixed_dt>>, ts->time_axis().f, ts, ta);
            case time_axis::generic_dt::CALENDAR: return make_interval(average_value<ts_src<time_axis::calendar_dt>>, ts->time_axis().c, ts, ta);
            case time_axis::generic_dt::POINT:    return make_interval(average_value<ts_src<time_axis::point_dt>>, ts->time_axis().p, ts, ta);
            }
            return make_interval(average_value<ts_src<time_axis::generic_dt>>, ts->time_axis(), ts, ta);
        }

        template<class S>
        inline double integral_value(const S& source, const utcperiod&p, size_t& last_idx, bool linear) {
            utctimespan tsum{ 0 };
            return accumulate_value<S>(source, p, last_idx, tsum, linear);
        }

        std::vector<double> integral_ts::values() const {
            switch (ts->time_axis().gt) { // pull out the real&fast time-axis before doing computations here:
            case time_axis::generic_dt::FIXED:    return make_interval(integral_value<ts_src<time_axis::fixed_dt>>, ts->time_axis().f, ts, ta);
            case time_axis::generic_dt::CALENDAR: return make_interval(integral_value<ts_src<time_axis::calendar_dt>>, ts->time_axis().c, ts, ta);
            case time_axis::generic_dt::POINT:    return make_interval(integral_value<ts_src<time_axis::point_dt>>, ts->time_axis().p, ts, ta);
            }
            return make_interval(integral_value<ts_src<time_axis::generic_dt>>, ts->time_axis(), ts, ta);
        }
        /** helpers */
        static inline utctimespan fixed_timestep(const time_axis::fixed_dt& ta) {
            return ta.dt;
        }
        static inline utctimespan fixed_timestep(const time_axis::point_dt& ta) {
            return utctimespan{0};
        }
        static inline utctimespan fixed_timestep(const time_axis::calendar_dt& ta) {
            return ta.dt<time_axis::calendar_dt::dt_tz_semantics?ta.dt:utctimespan{0};
        }
        static inline utctimespan fixed_timestep(const time_axis::generic_dt& ta) {
            if(ta.gt==time_axis::generic_dt::FIXED) {
                return ta.f.dt;
            } else if(ta.gt==time_axis::generic_dt::CALENDAR && ta.c.dt<time_axis::calendar_dt::dt_tz_semantics){
                return ta.c.dt;
            }
            return utctimespan{0};
        }

        /**We need a temporary slice of a time-axis */
        template <class TA>
        struct ta_slice {
            const TA& ta;
            size_t i0;
            size_t n;
            ta_slice(const TA&ta,size_t i0,size_t n):ta(ta),i0(i0),n(n){}
            size_t size() const {return n;}
            utcperiod period(size_t i) const {return ta.period(i0+i);}
            utctime time(size_t i) const {return ta.time(i0+i);}
        };

        template <class TA>
        static inline utctimespan fixed_timestep(const ta_slice<TA>& ta) {
            return fixed_timestep(ta.ta);
        }

        /** */
        template <class TA>
        void derivative_fx(const TA& ta,vector<double>& v,derivative_method dm) {
            if(v.size()<2) {
                if(v.size())
                    v.front()=isfinite(v.front())?0.0:nan;// flat or nan
            } else {
                auto dt=fixed_timestep(ta);                 // if fixed interval(really) ->highspeed
                if(dt>utctimespan{0}) {
                    switch(dm) {
                        case derivative_method::default_diff:
                        case derivative_method::center_diff:{// maybe average_diff is better name
                        double v0=v[0];
                        v[0]= isfinite(v0)?(isfinite(v[1])? (v[1]-v0)/to_seconds(2*dt) : 0.0) : nan;// first value special case
                        for(size_t i=1;i+1<v.size();++i) {
                            double v1=v[i];
                            if(isfinite(v1)) {
                                if(isfinite(v0)) {
                                    if(isfinite(v[i+1])) {
                                        v[i] = (v[i+1] - v0)/to_seconds(2*dt);
                                    } else {
                                        v[i] = (v1 - v0)/to_seconds(2*dt);
                                    }
                                } else {
                                    if(isfinite(v[i+1])) {
                                        v[i] = (v[i+1] - v1)/to_seconds(2*dt);
                                    } else {
                                        v[i] = 0.0;
                                    }
                                }
                            } else {
                                v[i]=nan;
                            }
                            v0=v1;
                        }
                        v.back()=isfinite(v.back())?(isfinite(v0)?(v.back() - v0)/to_seconds(2*dt):0.0):nan;// last value special case

                        }break;
                        case derivative_method::forward_diff:{// center forward diff
                        for(size_t i=0;i+1<v.size();++i) {
                            v[i]= isfinite(v[i])?(isfinite(v[i+1])?(v[i+1]-v[i])/to_seconds(dt):0.0):nan;
                        }
                        v.back()=isfinite(v.back())?0.0:nan;
                        } break;
                        case derivative_method::backward_diff:{// center backward diff
                        for(size_t i=v.size()-1;i>0;--i) {
                            v[i]= isfinite(v[i])?(isfinite(v[i-1])?(v[i]-v[i-1])/to_seconds(dt):0.0):nan;
                        }
                        v.front()=isfinite(v.front())?0.0:nan;
                        } break;
                    }
                } else { // variable intervals
                    switch(dm) {
                        case derivative_method::default_diff:
                        case derivative_method::center_diff:{// maybe average_diff is better name
                        double v0=v[0];
                        auto p0=ta.period(0);
                        auto p1=ta.period(1);
                        v[0]= isfinite(v0)?(isfinite(v[1])?(v[1]-v[0])/to_seconds(p1.end - p0.start):0.0):nan;// first value special case
                        for(size_t i=1;i+1<v.size();++i) {
                            p1=ta.period(i);
                            auto v1=v[i];
                            auto p2=ta.period(i+1);
                            auto v2=v[i+1];
                            if(isfinite(v1)){
                                if(isfinite(v0)) {// got left valid value
                                    if(isfinite(v2)) { // got right valid value
                                        v[i] = 2*(v2-v0)/to_seconds((p2.start-p0.start) +(p2.end-p0.end));
                                    } else { // right is nan, finish flat last half
                                        v[i] = (v1-v0)/to_seconds(p1.end - p0.start);
                                    }
                                } else if(isfinite(v2)) {// valid right side, nan on the left
                                    v[i] = (v2-v1)/to_seconds(p2.end - p1.start);
                                } else { // both sides nan, flat
                                    v[i] = 0.0;
                                }
                            } else {
                                v[i] = nan;
                            }
                            v0=v1;
                            p0=p1;
                        }
                        // we got p0,v0,
                        p1=ta.period(v.size()-1);
                        auto v1=v.back();
                        if(isfinite(v1)) {
                            if(isfinite(v0)) {
                                v.back() = (v1-v0)/to_seconds(p1.end-p0.start);
                            } else {
                                v.back() = 0.0;
                            }
                        } else {
                            v.back()=nan;
                        }

                        }break;
                        case derivative_method::forward_diff:{// center forward diff
                            auto p0 = ta.period(0);
                            auto v0 = v[0];
                            for(size_t i=0;i<v.size()-1;++i) {
                                auto p1=ta.period(i+1);
                                auto v1=v[i+1];
                                if(isfinite(v0)) {
                                    if(isfinite(v1)) {
                                        v[i] = 2*(v1 - v0)/ to_seconds((p1.start-p0.start) + (p1.end-p0.end));
                                    } else {
                                        v[i] = 0.0; // flat out
                                    }
                                } else {
                                    v[i]=nan;
                                }
                                p0=p1;
                                v0=v1;
                            }
                            v.back()=isfinite(v.back())?0.0:nan;// flat out last step.
                        } break;
                        case derivative_method::backward_diff:{// center backward diff
                            auto p0 = ta.period(0);
                            auto v0 = v[0];
                            v.front()=isfinite(v.front())?0.0:nan;
                            for(size_t i=1;i<v.size();++i) {
                                auto p1=ta.period(i);
                                auto v1=v[i];
                                if(isfinite(v1)) {
                                    if(isfinite(v0)) {
                                        v[i] = 2*(v1 - v0)/ to_seconds((p1.start-p0.start) + (p1.end-p0.end));
                                    } else {
                                        v[i] = 0.0; // flat out
                                    }
                                } else {
                                    v[i]=nan;
                                }
                                p0=p1;
                                v0=v1;
                            }
                        } break;
                    }
                }
            }
        }

        std::vector<double> derivative_ts::values() const {
            if(!ts) throw runtime_error("derivative of null ts attempted");
            auto v = ts->values();
            if(ts->point_interpretation()==POINT_INSTANT_VALUE) {
                for(size_t i=0;i+1<v.size();++i) {
                    v[i]= (v[i+1]-v[i])/to_seconds(ts->time(i+1)-ts->time(i));
                }
                if(v.size())
                    v[v.size()-1] = nan;
            } else {  // implement dm for stair-case
                derivative_fx(ts->time_axis(),v,dm);
            }
            return v;
        }
        double derivative_ts::value(size_t i) const {
            if(!ts) throw runtime_error("derivative of null ts attempted");
            if(ts->point_interpretation()==POINT_INSTANT_VALUE) {
                if(i+1<ts->size()) {
                    return (ts->value(i+1)-ts->value(i))/to_seconds(ts->time_axis().period(i).timespan());
                } else {
                    return nan;
                }
            } else {
                vector<double> v;v.reserve(3);
                size_t i0=i;
                if(i>0){ v.push_back(ts->value(i-1));i0=i-1;}
                v.push_back(ts->value(i));
                if(i+1<ts->size()) v.push_back(ts->value(i+1));
                //const auto & ta = ts->time_axis();
                ta_slice<decltype(ts->time_axis())> ta(ts->time_axis(),i0,v.size());
                derivative_fx(ta,v,dm);
                return v[i>0?1:0];
            }
            return nan;
        }

        double derivative_ts::value_at(utctime t) const {
            if(!ts) throw runtime_error("derivative of null ts attempted");
            size_t ix=ts->index_of(t); // easy, just figure out the ix of the period
            if(ix == string::npos)
                return nan;
            return value(ix);// and forward it to the value method.
        }

        apoint_ts apoint_ts::derivative(derivative_method dm) const {
            return apoint_ts(make_shared<derivative_ts>(ts,dm));
        }
        // implement popular ct for apoint_ts to make it easy to expose & use
        apoint_ts::apoint_ts(const time_axis::generic_dt& ta, double fill_value, ts_point_fx point_fx)
            :ts(std::make_shared<gpoint_ts>(ta, fill_value, point_fx)) {
        }
        apoint_ts::apoint_ts(const time_axis::generic_dt& ta, const std::vector<double>& values, ts_point_fx point_fx)
            : ts(std::make_shared<gpoint_ts>(ta, values, point_fx)) {
        }
        apoint_ts::apoint_ts(const time_axis::generic_dt& ta, std::vector<double>&& values, ts_point_fx point_fx)
            : ts(std::make_shared<gpoint_ts>(ta, std::move(values), point_fx)) {
        }

        apoint_ts::apoint_ts(std::string ref_ts_id)
            : ts(std::make_shared<aref_ts>(ref_ts_id)) {
        }
        apoint_ts::apoint_ts(std::string ref_ts_id, const apoint_ts&bts)
            : ts(std::make_shared<aref_ts>(ref_ts_id)) {
            bind(bts);// bind the symbolic ts directly
        }

        void apoint_ts::bind(const apoint_ts& bts) {
            if (!dynamic_pointer_cast<aref_ts>(ts))
                throw runtime_error("this time-series is not bindable");
            if (dynamic_pointer_cast<gpoint_ts>(bts.ts)) {
                dynamic_pointer_cast<aref_ts>(ts)->rep = dynamic_pointer_cast<gpoint_ts>(bts.ts);
            } else if (!bts.needs_bind()) {
                dynamic_pointer_cast<aref_ts>(ts)->rep = make_shared<gpoint_ts>(bts.time_axis(), bts.values(), bts.point_interpretation());
            } else {
                throw runtime_error("the supplied argument time-series must be a point ts or something that directly resolves to one");
            }

        }
        string apoint_ts::id() const {
            if (!dynamic_pointer_cast<aref_ts>(ts))
                return string{};
            return dynamic_pointer_cast<aref_ts>(ts)->id;
        }

        // and python needs these:
        apoint_ts::apoint_ts(const time_axis::fixed_dt& ta, double fill_value, ts_point_fx point_fx)
            :apoint_ts(time_axis::generic_dt(ta), fill_value, point_fx) {
        }
        apoint_ts::apoint_ts(const time_axis::fixed_dt& ta, const std::vector<double>& values, ts_point_fx point_fx)
            : apoint_ts(time_axis::generic_dt(ta), values, point_fx) {
        }

        // and python needs these:
        apoint_ts::apoint_ts(const time_axis::point_dt& ta, double fill_value, ts_point_fx point_fx)
            : apoint_ts(time_axis::generic_dt(ta), fill_value, point_fx) {
        }
        apoint_ts::apoint_ts(const time_axis::point_dt& ta, const std::vector<double>& values, ts_point_fx point_fx)
            : apoint_ts(time_axis::generic_dt(ta), values, point_fx) {
        }
        apoint_ts::apoint_ts(const rts_t &rts) :
            apoint_ts(time_axis::generic_dt(rts.ta), rts.v, rts.point_interpretation()) {
        }

        apoint_ts::apoint_ts(const vector<double>& pattern, utctimespan dt, const time_axis::generic_dt& ta) :
            apoint_ts(make_shared<periodic_ts>(pattern, dt, ta)) {
        }
        apoint_ts::apoint_ts(const vector<double>& pattern, utctimespan dt, utctime pattern_t0, const time_axis::generic_dt& ta) :
            apoint_ts(make_shared<periodic_ts>(pattern, dt, pattern_t0, ta)) {
        }

        apoint_ts::apoint_ts(time_axis::generic_dt&& ta, std::vector<double>&& values, ts_point_fx point_fx)
            : ts(std::make_shared<gpoint_ts>(std::move(ta), std::move(values), point_fx)) {
        }
        apoint_ts::apoint_ts(time_axis::generic_dt&& ta, double fill_value, ts_point_fx point_fx)
            : ts(std::make_shared<gpoint_ts>(std::move(ta), fill_value, point_fx)) {
        }
        apoint_ts apoint_ts::average(const gta_t &ta) const {
            return shyft::time_series::dd::average(*this, ta);
        }
        apoint_ts apoint_ts::integral(const gta_t &ta) const {
            return shyft::time_series::dd::integral(*this, ta);
        }
        apoint_ts apoint_ts::accumulate(const gta_t &ta) const {
            return shyft::time_series::dd::accumulate(*this, ta);
        }
        apoint_ts apoint_ts::time_shift(utctimespan dt) const {
            return shyft::time_series::dd::time_shift(*this, dt);
        }
        apoint_ts apoint_ts::extend(
            const apoint_ts & ts,
            extend_ts_split_policy split_policy, extend_ts_fill_policy fill_policy,
            utctime split_at, double fill_value
        ) const {
            return shyft::time_series::dd::extend(
                *this, ts,
                split_policy, fill_policy,
                split_at, fill_value
            );
        }

        /** recursive function to dig out bind_info */
        static void find_ts_bind_info(const std::shared_ptr<ipoint_ts>&its, std::vector<ts_bind_info>&r) {
            using namespace shyft;
            if (its == nullptr)
                return;
            if (dynamic_pointer_cast<const aref_ts>(its)) {
                auto rts = dynamic_pointer_cast<const aref_ts>(its);
                if (rts)
                    r.push_back(ts_bind_info(rts->id, apoint_ts(its)));
                else
                    ;// maybe throw ?
            } else if (dynamic_pointer_cast<const average_ts>(its)) {
                find_ts_bind_info(dynamic_cast<const average_ts*>(its.get())->ts, r);
            } else if (dynamic_cast<const integral_ts*>(its.get())) {
                find_ts_bind_info(dynamic_cast<const integral_ts*>(its.get())->ts, r);
            } else if (dynamic_cast<const accumulate_ts*>(its.get())) {
                find_ts_bind_info(dynamic_cast<const accumulate_ts*>(its.get())->ts, r);
            } else if (dynamic_cast<const time_shift_ts*>(its.get())) {
                find_ts_bind_info(dynamic_cast<const time_shift_ts*>(its.get())->ts, r);
            } else if (dynamic_cast<const abin_op_ts*>(its.get())) {
                auto bin_op = dynamic_cast<const abin_op_ts*>(its.get());
                find_ts_bind_info(bin_op->lhs.ts, r);
                find_ts_bind_info(bin_op->rhs.ts, r);
            } else if (dynamic_cast<const abin_op_scalar_ts*>(its.get())) {
                auto bin_op = dynamic_cast<const abin_op_scalar_ts*>(its.get());
                find_ts_bind_info(bin_op->rhs.ts, r);
            } else if (dynamic_cast<const abin_op_ts_scalar*>(its.get())) {
                auto bin_op = dynamic_cast<const abin_op_ts_scalar*>(its.get());
                find_ts_bind_info(bin_op->lhs.ts, r);
            } else if (dynamic_cast<const abs_ts*>(its.get())) {
                find_ts_bind_info(dynamic_cast<const abs_ts*>(its.get())->ts, r);
            } else if (dynamic_cast<const extend_ts*>(its.get())) {
                auto ext = dynamic_cast<const extend_ts*>(its.get());
                find_ts_bind_info(ext->lhs.ts, r);
                find_ts_bind_info(ext->rhs.ts, r);
            } else if (dynamic_cast<const ice_packing_ts*>(its.get())) {
                find_ts_bind_info(dynamic_cast<const ice_packing_ts*>(its.get())->ts.temp_ts.ts, r);
            } else if (dynamic_cast<const ice_packing_recession_ts*>(its.get())) {
                auto iprt = dynamic_cast<const ice_packing_recession_ts*>(its.get());
                find_ts_bind_info(iprt->flow_ts.ts, r);
                find_ts_bind_info(iprt->ice_packing_ts.ts, r);
            } else if (dynamic_cast<const rating_curve_ts*>(its.get())) {
                find_ts_bind_info(dynamic_cast<const rating_curve_ts*>(its.get())->ts.level_ts.ts, r);
            } else if (dynamic_cast<const krls_interpolation_ts*>(its.get())) {
                find_ts_bind_info(dynamic_cast<const krls_interpolation_ts*>(its.get())->ts.ts, r);
            } else if (dynamic_cast<const qac_ts*>(its.get())) {
                find_ts_bind_info(dynamic_cast<const qac_ts*>(its.get())->ts, r);
                find_ts_bind_info(dynamic_cast<const qac_ts*>(its.get())->cts, r);
            } else if (dynamic_cast<const inside_ts*>(its.get())) {
                find_ts_bind_info(dynamic_cast<const inside_ts*>(its.get())->ts, r);
            }
        }

        std::vector<ts_bind_info> apoint_ts::find_ts_bind_info() const {
            std::vector<ts_bind_info> r;
            shyft::time_series::dd::find_ts_bind_info(ts, r);
            return r;
        }

        ats_vector apoint_ts::partition_by(const calendar& cal, utctime t, utctimespan partition_interval, size_t n_partitions, utctime common_t0) const {
            // some very rudimentary argument checks:
            if (n_partitions < 1)
                throw std::runtime_error("n_partitions should be > 0");
            if (partition_interval.count() <= 0)
                throw std::runtime_error("partition_interval should be > 0, typically Calendar::YEAR|MONTH|WEEK|DAY");
            auto mk_raw_time_shift = [](const apoint_ts& ts, utctimespan dt)->apoint_ts {
                return apoint_ts(std::make_shared<time_shift_ts>(ts, dt));
            };
            auto r = shyft::time_series::partition_by<apoint_ts>(*this, cal, t, partition_interval, n_partitions, common_t0, mk_raw_time_shift);
            return ats_vector(r.begin(), r.end());
        }

        void apoint_ts::set(size_t i, double x) {
            gpoint_ts *gpts = dynamic_cast<gpoint_ts*>(ts.get());
            if (!gpts)
                throw std::runtime_error("apoint_ts::set(i,x) only allowed for ts of non-expression types");
            gpts->set(i, x);
        }
        void apoint_ts::fill(double x) {
            gpoint_ts *gpts = dynamic_cast<gpoint_ts*>(ts.get());
            if (!gpts)
                throw std::runtime_error("apoint_ts::fill(x) only allowed for ts of non-expression types");
            gpts->fill(x);
        }
        void apoint_ts::scale_by(double x) {
            gpoint_ts *gpts = dynamic_cast<gpoint_ts*>(ts.get());
            if (!gpts)
                throw std::runtime_error("apoint_ts::scale_by(x) only allowed for ts of non-expression types");
            gpts->scale_by(x);
        }

        bool apoint_ts::operator==(const apoint_ts& other) const {
            if (ts.get() == other.ts.get()) // equal by reference
                return true;
            // the time-consuming part, equal by value
            // SiH: I am not sure if this is ever useful, nevertheless it's doable
            //      given that it's ok with the zero-limit(not generally applicable!)

            const double zero_limit = 1e-9;// in hydrology, this is usually a small number
            if (ts->size() != other.ts->size())
                return false;
            for (size_t i = 0; i < ts->size(); ++i) {
                if (ts->time_axis().period(i) != other.ts->time_axis().period(i))
                    return false;
                if (fabs(ts->value(i) - other.ts->value(i)) > zero_limit)
                    return false;
            }
            return true;
        }

        apoint_ts apoint_ts::max(double a) const { return shyft::time_series::dd::max(*this, a); }
        apoint_ts apoint_ts::min(double a) const { return shyft::time_series::dd::min(*this, a); }

        apoint_ts apoint_ts::max(const apoint_ts& other) const { return shyft::time_series::dd::max(*this, other); }
        apoint_ts apoint_ts::min(const apoint_ts& other) const { return shyft::time_series::dd::min(*this, other); }

        apoint_ts apoint_ts::max(const apoint_ts &a, const apoint_ts&b) { return shyft::time_series::dd::max(a, b); }
        apoint_ts apoint_ts::min(const apoint_ts &a, const apoint_ts&b) { return shyft::time_series::dd::min(a, b); }

        apoint_ts apoint_ts::pow(const apoint_ts &a, const apoint_ts&b) { return shyft::time_series::dd::pow(a, b); }
        apoint_ts apoint_ts::pow(double a) const { return shyft::time_series::dd::pow(*this, a); }
        apoint_ts apoint_ts::pow(const apoint_ts& other) const { return shyft::time_series::dd::pow(*this, other); }

        apoint_ts apoint_ts::convolve_w(const std::vector<double> &w, shyft::time_series::convolve_policy conv_policy) const {
            return apoint_ts(std::make_shared<convolve_w_ts>(*this, w, conv_policy));
        }
        apoint_ts apoint_ts::slice(int i0, int n) const {
            gpoint_ts *gpts = dynamic_cast<gpoint_ts*>(ts.get());
            if (!gpts)
                throw std::runtime_error("apoint_ts::slice() only allowed for ts of non-expression types");
            return apoint_ts(make_shared<gpoint_ts>(gpts->slice(i0, n)));
        }

        apoint_ts apoint_ts::rating_curve(const rating_curve_parameters & rc_param) const {
            return apoint_ts(std::make_shared<rating_curve_ts>(*this, rc_param));
        }

        apoint_ts apoint_ts::ice_packing(const ice_packing_parameters & ip_param, ice_packing_temperature_policy ipt_policy) const {
            return apoint_ts(std::make_shared<ice_packing_ts>(*this, ip_param, ipt_policy));
        }

        apoint_ts apoint_ts::ice_packing_recession(
            const apoint_ts & ice_packing_ts,
            const ice_packing_recession_parameters & ipr_param
        ) const {
            return apoint_ts(std::make_shared<ice_packing_recession_ts>(
                *this, ice_packing_ts, ipr_param
            ));
        }

        apoint_ts apoint_ts::krls_interpolation(core::utctimespan dt, double rbf_gamma, double tol, std::size_t size) const {
            return apoint_ts(std::make_shared<krls_interpolation_ts>(*this, dt, rbf_gamma, tol, size));
        }
        prediction::krls_rbf_predictor apoint_ts::get_krls_predictor(core::utctimespan dt, double rbf_gamma, double tol, std::size_t size) const {
            if (needs_bind())
                throw std::runtime_error("cannot get predictor for unbound ts");
            shyft::prediction::krls_rbf_predictor predictor{ dt, rbf_gamma, tol, size };
            predictor.train(*this);
            return predictor;
        }
        apoint_ts apoint_ts::merge_points(const apoint_ts& o) {
            if(!o.ts)
                return *this;
            if(!ts) { // we are empty, create a new point ts
                auto a=make_shared<gpoint_ts>();
                ts_point_merge(a->rep,o);
                ts=a;
            } else {// we have existing ts, verify type:
                if(auto a=dynamic_pointer_cast<gpoint_ts>(ts)) {
                    ts_point_merge(a->rep,o);// merge
                } else if(auto a=dynamic_pointer_cast<aref_ts>(ts)) {
                    if(a->rep) { // its already bound
                        ts_point_merge(a->rep->rep,o);
                    } else {
                        auto r=make_shared<gpoint_ts>();
                        ts_point_merge(r->rep,o);
                        a->rep=r;
                    }
                } else {
                    throw runtime_error("self.merge_points_from:self ts must be a concrete point ts");
                }
            }
            return *this;
        }

        template<class TSV>
        static bool all_same_generic_time_axis_type(const TSV&tsv) {
            if (tsv.size() == 0) return true;
            auto gt = tsv[0].time_axis().gt;
            for (const auto&ts : tsv) {
                if (ts.time_axis().gt != gt)
                    return false;
            }
            return true;
        }

        std::vector<apoint_ts> percentiles(const std::vector<apoint_ts>& tsv1, const gta_t& ta, const intv_t& percentile_list) {
            std::vector<apoint_ts> r; r.reserve(percentile_list.size());
            auto tsvx = deflate_ts_vector<gts_t>(tsv1);
            // check of all tsvx.time_axis is of same type
            //  make that vector type (move values), and run percentile calc for that
            //  the objective is to avoid traffic over the 'switch' in generic_dt
            //  and provide direct access to representative time-axis instead
            if (all_same_generic_time_axis_type(tsvx)) {
                if (tsvx.size()) {
                    auto gt = tsvx[0].time_axis().gt;
                    switch (gt) {
                    case time_axis::generic_dt::FIXED: {
                        std::vector<point_ts<time_axis::fixed_dt>> tsv; tsv.reserve(tsvx.size());
                        for (auto&ts : tsvx) tsv.emplace_back(move(ts.ta.f), move(ts.v), ts.fx_policy);
                        auto rp = shyft::time_series::calculate_percentiles(ta, tsv, percentile_list);
                        for (auto&ts : rp) r.emplace_back(ta, std::move(ts.v), POINT_AVERAGE_VALUE);
                    } break;
                    case time_axis::generic_dt::CALENDAR: {
                        std::vector<point_ts<time_axis::calendar_dt>> tsv; tsv.reserve(tsvx.size());
                        for (auto&ts : tsvx) tsv.emplace_back(move(ts.ta.c), move(ts.v), ts.fx_policy);
                        auto rp = shyft::time_series::calculate_percentiles(ta, tsv, percentile_list);
                        for (auto&ts : rp) r.emplace_back(ta, std::move(ts.v), POINT_AVERAGE_VALUE);
                    } break;
                    case time_axis::generic_dt::POINT: {
                        std::vector<point_ts<time_axis::point_dt>> tsv; tsv.reserve(tsvx.size());
                        for (auto&ts : tsvx) tsv.emplace_back(move(ts.ta.p), move(ts.v), ts.fx_policy);
                        auto rp = shyft::time_series::calculate_percentiles(ta, tsv, percentile_list);
                        for (auto&ts : rp) r.emplace_back(ta, std::move(ts.v), POINT_AVERAGE_VALUE);
                    } break;
                    }
                } else {
                    for (size_t i = 0; i < percentile_list.size(); ++i)
                        r.emplace_back(ta, shyft::nan, POINT_AVERAGE_VALUE);
                }
            } else {
                auto rp = shyft::time_series::calculate_percentiles(ta, tsvx, percentile_list);
                for (auto&ts : rp) r.emplace_back(ta, std::move(ts.v), POINT_AVERAGE_VALUE);
            }
            return r;
        }

        std::vector<apoint_ts> percentiles(const std::vector<apoint_ts>& ts_list, const time_axis::fixed_dt& ta, const intv_t& percentile_list) {
            return percentiles(ts_list, time_axis::generic_dt(ta), percentile_list);
        }

        double abin_op_ts::value(size_t i) const {
            if (i == std::string::npos || i >= time_axis().size())
                return nan;
            return value_at(time_axis().time(i));
        }
        double abin_op_scalar_ts::value_at(utctime t) const {
            bind_check();
            return do_op(lhs, op, rhs(t));
        }
        double abin_op_scalar_ts::value(size_t i) const {
            bind_check();
            return do_op(lhs, op, rhs.value(i));
        }

        std::vector<double> abin_op_scalar_ts::values() const {
            bind_check();
            const vector<double> *rhs_v{ terminal_values(rhs) };
            if (rhs_v) {
                const auto& r_v = *rhs_v;
                std::vector<double> r; r.reserve(r_v.size());
                auto l = lhs;
                switch (op) {
                case OP_ADD:for (const auto&v : r_v) r.emplace_back(v + l); return r;
                case OP_SUB:for (const auto&v : r_v) r.emplace_back(l - v); return r;
                case OP_MUL:for (const auto&v : r_v) r.emplace_back(l*v); return r;
                case OP_DIV:for (const auto&v : r_v) r.emplace_back(l / v); return r;
                case OP_MAX:for (const auto&v : r_v) r.emplace_back(std::max(v, l)); return r;
                case OP_MIN:for (const auto&v : r_v) r.emplace_back(std::min(v, l)); return r;
                case OP_POW:for (const auto&v : r_v) r.emplace_back(std::pow(l, v)); return r;
                default: throw runtime_error("Unsupported operation " + to_string(int(op)));
                }
            } else {
                std::vector<double> r(rhs.values());
                auto l = lhs;
                switch (op) {
                case OP_ADD:for (size_t i = 0; i < r.size(); ++i) r[i] += l; return r;
                case OP_SUB:for (size_t i = 0; i < r.size(); ++i) r[i] = l - r[i]; return r;
                case OP_MUL:for (size_t i = 0; i < r.size(); ++i) r[i] *= l; return r;
                case OP_DIV:for (size_t i = 0; i < r.size(); ++i) r[i] = l / r[i]; return r;
                case OP_MAX:for (size_t i = 0; i < r.size(); ++i) r[i] = std::max(r[i], l); return r;
                case OP_MIN:for (size_t i = 0; i < r.size(); ++i) r[i] = std::min(r[i], l); return r;
                case OP_POW:for (size_t i = 0; i < r.size(); ++i) r[i] = std::pow(l,r[i]); return r;
                default: throw runtime_error("Unsupported operation " + to_string(int(op)));
                }
            }
        }

        double abin_op_ts_scalar::value_at(utctime t) const {
            bind_check();
            return do_op(lhs(t), op, rhs);
        }
        double abin_op_ts_scalar::value(size_t i) const {
            bind_check();
            return do_op(lhs.value(i), op, rhs);
        }
        std::vector<double> abin_op_ts_scalar::values() const {
            bind_check();
            const vector<double>* lhs_v{ terminal_values(lhs) };
            if (lhs_v) { // avoid a copy, but does not help much..
                std::vector<double> r; r.reserve(lhs_v->size());
                auto rv = rhs;
                switch (op) {
                case OP_ADD:for (const auto&lv : *lhs_v) r.emplace_back(lv + rv); return r;
                case OP_SUB:for (const auto&lv : *lhs_v) r.emplace_back(lv - rv); return r;
                case OP_MUL:for (const auto&lv : *lhs_v) r.emplace_back(lv*rv); return r;
                case OP_DIV:for (const auto&lv : *lhs_v) r.emplace_back(lv / rv); return r;
                case OP_MAX:for (const auto&lv : *lhs_v) r.emplace_back(std::max(lv, rv)); return r;
                case OP_MIN:for (const auto&lv : *lhs_v) r.emplace_back(std::min(lv, rv)); return r;
                case OP_POW:for (const auto&lv : *lhs_v) r.emplace_back(std::pow(lv, rv)); return r;
                default: throw runtime_error("Unsupported operation " + to_string(int(op)));
                }
            } else {
                std::vector<double> l(lhs.values());
                auto r = rhs;
                switch (op) {
                case OP_ADD:for (size_t i = 0; i < l.size(); ++i) l[i] += r; return l;
                case OP_SUB:for (size_t i = 0; i < l.size(); ++i) l[i] -= r; return l;
                case OP_MUL:for (size_t i = 0; i < l.size(); ++i) l[i] *= r; return l;
                case OP_DIV:for (size_t i = 0; i < l.size(); ++i) l[i] /= r; return l;
                case OP_MAX:for (size_t i = 0; i < l.size(); ++i) l[i] = std::max(l[i], r); return l;
                case OP_MIN:for (size_t i = 0; i < l.size(); ++i) l[i] = std::min(l[i], r); return l;
                case OP_POW:for (size_t i = 0; i < l.size(); ++i) l[i] = std::pow(l[i], r); return l;
                default: throw runtime_error("Unsupported operation " + to_string(int(op)));
                }
            }
        }

        apoint_ts time_shift(const apoint_ts& ts, utctimespan dt) {
            return apoint_ts(std::make_shared<time_shift_ts>(ts, dt));
        }

        apoint_ts apoint_ts::abs() const {
            return apoint_ts(std::make_shared<abs_ts>(ts));
        }

        apoint_ts apoint_ts::min_max_check_linear_fill(double min_x, double max_x, utctimespan max_dt) const {
            return apoint_ts(make_shared<qac_ts>(*this, qac_parameter::create_min_max_linear_interpolation_parameters(
                /* nan_qa */ true, /* infinity_qa */ true, min_x, max_x, max_dt
            )));
        }

        apoint_ts apoint_ts::min_max_check_ts_fill(double min_x, double max_x, utctimespan max_dt, const apoint_ts& cts) const {
            return apoint_ts(make_shared<qac_ts>(*this,qac_parameter::create_min_max_linear_interpolation_parameters(
                /* nan_qa */ true, /* infinity_qa */ true, min_x, max_x, max_dt
            ), cts));
        }

        apoint_ts apoint_ts::inside(double min_v,double max_v,double nan_v,double inside_v,double outside_v) const {
            return apoint_ts(make_shared<inside_ts>(*this, inside_parameter{ min_v,max_v,nan_v,inside_v,outside_v }));
        }

        apoint_ts apoint_ts::decode(int start_bit,int n_bits) const {
            if(start_bit <0 || start_bit >51)
                throw runtime_error("start_bit must be in range [0..51], was " + to_string(start_bit));
            if(n_bits < 1 || start_bit + n_bits > 51)
                throw runtime_error("n_bits must be > 0 and start_bit+n_bits <= 51: n_bits =" + to_string(n_bits) + ", start_bit="+to_string(start_bit));
            return apoint_ts(make_shared<decode_ts>(*this, bit_decoder{(unsigned int)start_bit,(unsigned int)n_bits }));
        }

        double nash_sutcliffe(const apoint_ts& observation_ts, const apoint_ts& model_ts, const gta_t &ta) {
            average_accessor<apoint_ts, gta_t> o(observation_ts, ta);
            average_accessor<apoint_ts, gta_t> m(model_ts, ta);
            return 1.0 - shyft::time_series::nash_sutcliffe_goal_function(o, m);
        }

        double kling_gupta(const apoint_ts & observation_ts, const apoint_ts &model_ts, const gta_t & ta, double s_r, double s_a, double s_b) {
            average_accessor<apoint_ts, gta_t> o(observation_ts, ta);
            average_accessor<apoint_ts, gta_t> m(model_ts, ta);
            return 1.0 - shyft::time_series::kling_gupta_goal_function<dlib::running_scalar_covariance<double>>(o, m, s_r, s_a, s_b);
        }
        // glacier_melt_ts as apoint_ts with it's internal being a glacier_melt_ts
        struct aglacier_melt_ts :ipoint_ts {
            glacier_melt_ts<std::shared_ptr<ipoint_ts>> gm;
            //-- default stuff, ct/copy etc goes here
            aglacier_melt_ts() = default;

            //-- useful ct goes here
            aglacier_melt_ts(const apoint_ts& temp, const apoint_ts& sca_m2, double glacier_area_m2, double dtf) :
                gm(temp.ts, sca_m2.ts, glacier_area_m2, dtf) {
            }
            //aglacier_melt_ts(apoint_ts&& ats, utctimespan dt):ts(std::move(ats.ts)),ta(time_axis::time_shift(ats.time_axis(),dt)),dt(dt) {}
            //aglacier_melt_ts(const std::shared_ptr<ipoint_ts> &ts, utctime dt ):ts(ts),ta(time_axis::time_shift(ts->time_axis(),dt)),dt(dt){}

            // implement ipoint_ts contract:
            virtual ts_point_fx point_interpretation() const { return gm.fx_policy; }
            virtual void set_point_interpretation(ts_point_fx point_interpretation) { gm.fx_policy = point_interpretation; }
            virtual const gta_t& time_axis() const { return gm.time_axis(); }
            virtual utcperiod total_period() const { return gm.time_axis().total_period(); }
            virtual size_t index_of(utctime t) const { return gm.time_axis().index_of(t); }
            virtual size_t size() const { return gm.time_axis().size(); }
            virtual utctime time(size_t i) const { return gm.time_axis().time(i); };
            virtual double value(size_t i) const { return gm.value(i); }
            virtual double value_at(utctime t) const { return gm(t); }
            virtual std::vector<double> values() const {
                std::vector<double> r; r.reserve(size());
                for (size_t i = 0; i < size(); ++i) r.push_back(value(i));
                return r;
            }
            virtual bool needs_bind() const { return gm.temperature->needs_bind() || gm.sca_m2->needs_bind(); }
            virtual void do_bind() { gm.temperature->do_bind(); gm.sca_m2->do_bind(); }

        };
        apoint_ts create_glacier_melt_ts_m3s(const apoint_ts & temp, const apoint_ts& sca_m2, double glacier_area_m2, double dtf) {
            return apoint_ts(make_shared<aglacier_melt_ts>(temp, sca_m2, glacier_area_m2, dtf));
        }



        std::vector<char> apoint_ts::serialize_to_bytes() const {
            auto ss = serialize();
            return std::vector<char>(std::begin(ss), std::end(ss));
        }
        apoint_ts apoint_ts::deserialize_from_bytes(const std::vector<char>&ss) {
            return deserialize(std::string(ss.data(), ss.size()));
        }

        //--ats_vector impl.
        // multiply operators
        ats_vector operator*(ats_vector const &a, double b) { ats_vector r; r.reserve(a.size()); for (auto const&ts : a) r.push_back(ts*b); return r; }
        ats_vector operator*(double a, ats_vector const &b) { return b * a; }
        ats_vector operator*(ats_vector const &a, ats_vector const& b) {
            if (a.size() != b.size()) throw runtime_error(string("ts-vector multiply require same sizes: lhs.size=") + std::to_string(a.size()) + string(",rhs.size=") + std::to_string(b.size()));
            ats_vector r; r.reserve(a.size()); for (size_t i = 0; i < a.size(); ++i) r.push_back(a[i] * b[i]);
            return r;
        }
        ats_vector operator*(ats_vector::value_type const &a, ats_vector const& b) { ats_vector r; r.reserve(b.size()); for (size_t i = 0; i < b.size(); ++i) r.push_back(a*b[i]); return r; }
        ats_vector operator*(ats_vector const& b, ats_vector::value_type const &a) { return a * b; }


        // divide operators
        ats_vector operator/(ats_vector const &a, double b) { return a * (1.0 / b); }
        ats_vector operator/(double a, ats_vector const &b) { ats_vector r; r.reserve(b.size()); for (auto const&ts : b) r.push_back(a / ts); return r; }
        ats_vector operator/(ats_vector const &a, ats_vector const& b) {
            if (a.size() != b.size()) throw runtime_error(string("ts-vector divide require same sizes: lhs.size=") + std::to_string(a.size()) + string(",rhs.size=") + std::to_string(b.size()));
            ats_vector r; r.reserve(a.size()); for (size_t i = 0; i < a.size(); ++i) r.push_back(a[i] / b[i]);
            return r;
        }
        ats_vector operator/(ats_vector::value_type const &a, ats_vector const& b) { ats_vector r; r.reserve(b.size()); for (size_t i = 0; i < b.size(); ++i) r.push_back(a / b[i]); return r; }
        ats_vector operator/(ats_vector const& b, ats_vector::value_type const &a) { ats_vector r; r.reserve(b.size()); for (size_t i = 0; i < b.size(); ++i) r.push_back(b[i] / a); return r; }

        // add operators
        ats_vector operator+(ats_vector const &a, double b) { ats_vector r; r.reserve(a.size()); for (auto const&ts : a) r.push_back(ts + b); return r; }
        ats_vector operator+(double a, ats_vector const &b) { return b + a; }
        ats_vector operator+(ats_vector const &a, ats_vector const& b) {
            // additional rule nice for reduce(add.. )
            // if one of tsv is 0, result is a tsv.
            if (a.size()==0 && b.size()!=0)
                return b;
            if (a.size()!=0 && b.size()==0)
                return a;
            if (a.size() != b.size()) throw runtime_error(string("ts-vector add require same sizes: lhs.size=") + std::to_string(a.size()) + string(",rhs.size=") + std::to_string(b.size()));
            ats_vector r; r.reserve(a.size()); for (size_t i = 0; i < a.size(); ++i) r.push_back(a[i] + b[i]);
            return r;
        }
        ats_vector operator+(ats_vector::value_type const &a, ats_vector const& b) { ats_vector r; r.reserve(b.size()); for (size_t i = 0; i < b.size(); ++i) r.push_back(a + b[i]); return r; }
        ats_vector operator+(ats_vector const& b, ats_vector::value_type const &a) { return a + b; }

        // sub operators
        ats_vector operator-(const ats_vector& a) { ats_vector r; r.reserve(a.size()); for (auto const&ts : a) r.push_back(-ts); return r; }

        ats_vector operator-(ats_vector const &a, double b) { ats_vector r; r.reserve(a.size()); for (auto const&ts : a) r.push_back(ts - b); return r; }
        ats_vector operator-(double a, ats_vector const &b) { ats_vector r; r.reserve(b.size()); for (auto const&ts : b) r.push_back(a - ts); return r; }
        ats_vector operator-(ats_vector const &a, ats_vector const& b) {
            if (a.size()==0 && b.size()!=0) {
                ats_vector r; r.reserve(b.size()); for (size_t i = 0; i < b.size(); ++i) r.push_back( - b[i]);
                return r;
            }
            if (a.size()!=0 && b.size()==0)
                return a;

            if (a.size() != b.size()) throw runtime_error(string("ts-vector sub require same sizes: lhs.size=") + std::to_string(a.size()) + string(",rhs.size=") + std::to_string(b.size()));
            ats_vector r; r.reserve(a.size()); for (size_t i = 0; i < a.size(); ++i) r.push_back(a[i] - b[i]);
            return r;
        }
        ats_vector operator-(ats_vector::value_type const &a, ats_vector const& b) { ats_vector r; r.reserve(b.size()); for (size_t i = 0; i < b.size(); ++i) r.push_back(a - b[i]); return r; }
        ats_vector operator-(ats_vector const& b, ats_vector::value_type const &a) { ats_vector r; r.reserve(b.size()); for (size_t i = 0; i < b.size(); ++i) r.push_back(b[i] - a); return r; }

        // max/min operators
        ats_vector ats_vector::min(ats_vector const& x) const {
            if (size() != x.size()) throw runtime_error(string("ts-vector min require same sizes: lhs.size=") + std::to_string(size()) + string(",rhs.size=") + std::to_string(x.size()));
            ats_vector r; r.reserve(size()); for (size_t i = 0; i < size(); ++i) r.push_back((*this)[i].min(x[i]));
            return r;
        }
        ats_vector ats_vector::max(ats_vector const& x) const {
            if (size() != x.size()) throw runtime_error(string("ts-vector max require same sizes: lhs.size=") + std::to_string(size()) + string(",rhs.size=") + std::to_string(x.size()));
            ats_vector r; r.reserve(size()); for (size_t i = 0; i < size(); ++i) r.push_back((*this)[i].max(x[i]));
            return r;
        }

        ats_vector ats_vector::pow(ats_vector const& x) const {
            if (size() != x.size()) throw runtime_error(string("ts-vector pow require same sizes: lhs.size=") + std::to_string(size()) + string(",rhs.size=") + std::to_string(x.size()));
            ats_vector r; r.reserve(size()); for (size_t i = 0; i < size(); ++i) r.push_back((*this)[i].pow(x[i]));
            return r;
        }

        ats_vector ats_vector::inside(double min_v,double max_v,double nan_v, double inside_v, double outside_v) const {
            ats_vector r; r.reserve(size()); for (size_t i = 0; i < size(); ++i) r.push_back((*this)[i].inside(min_v,max_v,nan_v,inside_v,outside_v));
            return r;
        }
        ats_vector min(ats_vector const &a, double b) { return a.min(b); }
        ats_vector min(double b, ats_vector const &a) { return a.min(b); }
        ats_vector min(ats_vector const &a, apoint_ts const& b) { return a.min(b); }
        ats_vector min(apoint_ts const &b, ats_vector const& a) { return a.min(b); }
        ats_vector min(ats_vector const &a, ats_vector const &b) { return a.min(b); }

        ats_vector max(ats_vector const &a, double b) { return a.max(b); }
        ats_vector max(double b, ats_vector const &a) { return a.max(b); }
        ats_vector max(ats_vector const &a, apoint_ts const & b) { return a.max(b); }
        ats_vector max(apoint_ts const &b, ats_vector const &a) { return a.max(b); }
        ats_vector max(ats_vector const &a, ats_vector const & b) { return a.max(b); }

        ats_vector pow(ats_vector const &a, double b) { return a.pow(b); }
        ats_vector pow(double b, ats_vector const &a) { 
            ats_vector r;
            r.reserve(a.size());
            for(const auto&ts:a) 
                r.push_back(pow(b,ts));
            return r; 
        }
        ats_vector pow(ats_vector const &a, apoint_ts const & b) { return a.pow(b); }
        ats_vector pow(apoint_ts const &b, ats_vector const &a) { 
            ats_vector r;
            r.reserve(a.size());
            for(const auto&ts:a) 
                r.push_back(pow(b,ts));
            return r;             
        }
        ats_vector pow(ats_vector const &a, ats_vector const & b) { return a.pow(b); }

        apoint_ts  ats_vector::forecast_merge(utctimespan lead_time, utctimespan fc_interval) const {
            //verify arguments
            if (lead_time.count() < 0)
                throw runtime_error("lead_time parameter should be 0 or a positive number giving number of seconds into each forecast to start the merge slice");
            if (fc_interval.count() <= 0)
                throw runtime_error("fc_interval parameter should be positive number giving number of seconds between first time point in each of the supplied forecast");
            for (size_t i = 1; i < size(); ++i) {
                if ((*this)[i - 1].total_period().start + fc_interval > (*this)[i].total_period().start) {
                    throw runtime_error(
                        string("The suplied forecast vector should be strictly ordered by increasing t0 by length at least fc_interval: requirement broken at index:")
                        + std::to_string(i)
                    );
                }
            }
            return time_series::forecast_merge<apoint_ts>(*this, lead_time, fc_interval);

        }
        double ats_vector::nash_sutcliffe(apoint_ts const &obs, utctimespan t0_offset, utctimespan dt, int n) const {
            if (n < 0)
                throw runtime_error("n, number of intervals, must be specified as > 0");
            if (dt <= utctimespan{0})
                throw runtime_error("dt, average interval, must be specified as > 0 s");
            if (t0_offset < utctimespan{0})
                throw runtime_error("lead_time,t0_offset,must be specified  >= 0 s");
            return time_series::nash_sutcliffe(*this, obs, t0_offset, dt, (size_t)n);
        }

        ats_vector ats_vector::average_slice(utctimespan t0_offset, utctimespan dt, int n) const {
            if (n < 0)
                throw runtime_error("n, number of intervals, must be specified as > 0");
            if (dt <= utctimespan{0})
                throw runtime_error("dt, average interval, must be specified as > 0 s");
            if (t0_offset < utctimespan{0})
                throw runtime_error("lead_time,t0_offset,must be specified  >= 0 s");
            ats_vector r;
            for (size_t i = 0; i < size(); ++i) {
                auto const& ts = (*this)[i];
                if (ts.size()) {
                    gta_t ta(ts.time_axis().time(0) + t0_offset, dt, n);
                    r.push_back((*this)[i].average(ta));
                } else {
                    r.push_back(ts);
                }
            }
            return r;
        }
        /** \see shyft::qm::quantile_map_forecast */
        ats_vector quantile_map_forecast(vector<ats_vector> const & forecast_sets,
            vector<double> const& set_weights,
            ats_vector const& historical_data,
            gta_t const&time_axis,
            utctime interpolation_start,
            utctime interpolation_end,
            bool interpolated_quantiles
        ) {
            // since this is scripting access, verify all parameters here
            if (forecast_sets.size() < 1)
                throw runtime_error("forecast_set must contain at least one forecast");
            if (historical_data.size() < 2)
                throw runtime_error("historical_data should have more than one time-series");
            if (forecast_sets.size() != set_weights.size())
                throw runtime_error(string("The size of weights (") + to_string(set_weights.size()) + string("), must match number of forecast-sets (") + to_string(forecast_sets.size()) + string(""));
            if (time_axis.size() == 0)
                throw runtime_error("time-axis should have at least one step");
            if (core::is_valid(interpolation_start)) {
                if (!time_axis.total_period().contains(interpolation_start)) {
                    calendar utc;
                    auto ts = utc.to_string(interpolation_start);
                    auto ps = utc.to_string(time_axis.total_period());
                    throw runtime_error("interpolation_start " + ts + " is not within time_axis period " + ps);
                }
                if (core::is_valid(interpolation_end) && !time_axis.total_period().contains(interpolation_end)) {
                    calendar utc;
                    auto ts = utc.to_string(interpolation_end);
                    auto ps = utc.to_string(time_axis.total_period());
                    throw runtime_error("interpolation_end " + ts + " is not within time_axis period " + ps);
                }
            }
            return qm::quantile_map_forecast<time_series::average_accessor<apoint_ts, gta_t> >(forecast_sets, set_weights, historical_data, time_axis, interpolation_start, interpolation_end, interpolated_quantiles);

        }

        apoint_ts extend(
            const apoint_ts & lhs_ts,
            const apoint_ts & rhs_ts,
            extend_ts_split_policy split_policy, extend_ts_fill_policy fill_policy,
            utctime split_at, double fill_value
        ) {
            return apoint_ts(std::make_shared<extend_ts>(
                lhs_ts, rhs_ts,
                split_policy, fill_policy,
                split_at, fill_value
                ));
        }

        std::vector<double> extend_ts::values() const {
            this->bind_check();

            const utctime split_at = this->get_split_at();
            const auto lhs_p = this->lhs.time_axis().total_period();
            const auto rhs_p = this->rhs.time_axis().total_period();

            // get values
            std::vector<double> lhs_values{}, rhs_values{};
            if (split_at >= lhs_p.start) lhs_values = this->lhs.values();
            if (split_at <= lhs_p.end)   rhs_values = this->rhs.values();

            // possibly to long, but not too short, all values default to nan
            std::vector<double> result;
            result.reserve(lhs.size() + rhs_values.size());

            auto res_oit = std::back_inserter(result);  // output iterator

            // values from the lhs
            if (split_at >= lhs_p.end) {  // use all of lhs
                res_oit = std::copy(lhs_values.begin(), lhs_values.end(), res_oit);
            } else if (split_at >= lhs_p.start) {  // split inside lhs
                size_t lhs_i = this->lhs.time_axis().index_of(split_at);
                auto lhs_end_it = lhs_values.begin();
                std::advance(lhs_end_it, lhs_i);

                res_oit = std::copy(lhs_values.begin(), lhs_end_it, res_oit);
            }

            // values from the rhs
            if (split_at <= rhs_p.start) {  // use all of rhs
                std::copy(rhs_values.begin(), rhs_values.end(), res_oit);
            } else if (split_at <= rhs_p.end) {  // split inside rhs
                size_t rhs_i = this->rhs.time_axis().index_of(split_at);
                auto rhs_start_it = rhs_values.begin();
                std::advance(rhs_start_it, rhs_i);

                std::copy(rhs_start_it, rhs_values.end(), res_oit);
            }

            return result;
        }

        double extend_ts::value_at(utctime t) const {
            //this->bind_check();  // done in time_axis()
            if (!this->time_axis().total_period().contains(t)) {
                return nan;
            }

            utctime split_at = this->get_split_at();

            if (t < split_at) {  // lhs
                if (this->lhs.time_axis().total_period().contains(t)) {
                    return this->lhs(t);
                } else {
                    // greater than lhs.end -> use policy
                    switch (this->ets_fill_p) {
                    default:
                    case EPF_NAN:  return nan;
                    case EPF_FILL: return this->fill_value;
                    case EPF_LAST: return this->lhs.value(lhs.size() - 1);
                    }
                }
            } else {  // rhs
                if (this->rhs.time_axis().total_period().contains(t)) {
                    return this->rhs(t);
                } else {
                    // less than rhs.start -> use policy
                    switch (this->ets_fill_p) {
                    default:
                    case EPF_NAN:  return nan;
                    case EPF_FILL: return this->fill_value;
                    case EPF_LAST: return this->lhs.value(lhs.size() - 1);
                    }
                }
            }
        }

        double extend_ts::value(size_t i) const {
            // this->bind_check();  // done in value_at()
            if (i == std::string::npos || i >= time_axis().size()) {
                return nan;
            }
            return value_at(time_axis().time(i));
        }

        // ================================================================================

        void qac_ts::_scan_repeating_values(const std::vector<double> & data) const {
            const std::size_t data_size = this->ts->size();

            // ready cache vector
            this->_repeating_cache.clear();
            this->_repeating_cache.reserve(data_size);

            // scan data
            double previous_value{};
            std::size_t repeat_counter{ 0u };
            std::size_t repeat_start_i{ 0u };
            for ( std::size_t i = 0u; i < data_size; ++i ) {
                this->_repeating_cache.emplace_back(0);

                double value = data[i];

                // after the first value
                if ( i > 0u ) {
                    if ( ::fabs(previous_value - value) <= this->p.repeat_tolerance ) {
                        // value close to previous
                        if ( repeat_counter == 0u ) {
                            repeat_counter = 1u;
                            repeat_start_i = i;
                        }
                        repeat_counter += 1u;
                    } else {
                        // value not close to previous -> check if this is a long enough repetition sequence
                        utctime t0 = this->time(repeat_start_i);
                        utctime t1 = this->time(i);  // not i-1 to count the full period of the last repeating value

                        // long enough period?
                        if ( t1 - t0 >= this->p.repeat_timespan ) {
                            // then replace the last few values 
                            for ( std::size_t j = repeat_start_i; j < i; ++j ) {
                                this->_repeating_cache[j] = 1;
                            }
                            repeat_counter = 0u;
                        }
                    }
                }
                previous_value = value;
            }

            // repeating at the end?
            if ( repeat_counter > 0 ) {
                utctime t0 = this->time(repeat_start_i);
                utctime t1 = this->total_period().end;  // end of series

                // long enough period?
                if ( t1 - t0 >= this->p.repeat_timespan ) {
                    // then replace the last few values 
                    for ( std::size_t j = repeat_start_i; j < data_size; ++j ) {
                        this->_repeating_cache[j] = 1;
                    }
                }
            }
        }

        // ----------------------------------------

        double qac_ts:: _fill_value(size_t i) const {
            // use a correction value ts if available
            if ( cts ) {
                return cts->value_at(ts->time(i)); // we do not check this value, assume ok!
            } else if ( this->p.can_do_linear_interpolation() ) {
                // do linear interpolation between previous--next *valid* point
                return this->_fill_linear_interpolation(i);
            } else {
                // fill a constant value if no other is set
                return this->_fill_constant();
            }
            
        }
        double qac_ts::_fill_constant() const {
            return this->p.constant_filler;
        }
        double qac_ts::_fill_linear_interpolation(std::size_t i) const {
            const auto t = ts->time(i);
            const std::size_t n = ts->size();

            // is a there a previous and next point, or is the allowed interpolation span empty?
            if ( i == 0 || i + 1 >= n || p.max_scan_timespan == utctimespan{0} ) {
                return shyft::nan;
            }

            // find previous _valid_ point
            std::size_t j = i;
            while ( j-- ) {
                utctime t0 = ts->time(j);

                // point within allowed timespan?
                if ( t - t0 > p.max_scan_timespan ) {
                    return shyft::nan;
                }
                // is value at previous point valid?
                double x0 = ts->value(j);
                if ( p.is_ok_quality(j, x0, *this) ) {
                    // found previous point, now scan for next point
                    for ( std::size_t k = i + 1; k < n; ++k ) {
                        utctime t1 = ts->time(k);

                        // distance between previous and next within allowed span?
                        if ( t1 - t0 > p.max_scan_timespan ) {
                            return shyft::nan;
                        }
                        // is value at next valid?
                        double x1 = ts->value(k);
                        if ( p.is_ok_quality(k, x1, *this) ) {
                            // do interpolation!
                            const double a = (x1 - x0) / to_seconds(t1 - t0);
                            const double b = x0 - a * to_seconds(t0);// x= a*t + b -> b= x- a*t
                            return a * to_seconds(t) + b;
                        }
                    }
                }
            }

            // failed to find a substitute value using interpolation
            return shyft::nan;
        }

        // ----------------------------------------


        double qac_ts::value(size_t i) const {

            // construct repeating values cache if neccessary
            if ( this->p.repeating_qa_enabled ) {
                if ( this->_repeating_cache.size() != this->size() ) {  // simple check to determine if repeating values are cached
                    auto r=ts->values();
                    this->_scan_repeating_values(r);
                }
            }

            // do we need to correct this value?
            double x = ts->value(i);
            if ( p.is_ok_quality(i, x, *this) )
                return x;

            // try to fill in replacement value
            return _fill_value(i);
        }
        double qac_ts::value_at(utctime t) const {
            size_t i = index_of(t);

            // time in ts span?
            if (i == string::npos) {
                return shyft::nan;
            }

            // if ts is average, then the value is constant over each period
            double x0 = value(i);
            if (ts->point_interpretation() == ts_point_fx::POINT_AVERAGE_VALUE) {
                return x0;
            }
            
            // else do linear interpolation over period
            // * this is not the same interpolation (_fill_linear_interpolation)
            //   as for correcting periods determined to be not-ok

            // - is value exactly at the start of the period?
            utctime t0 = time(i);
            if (t0 == t) {
                return x0;
            }
            // - is there a next point to interpolate to?
            if ( i + 1 >= size() ) {
                return shyft::nan;
            }
            // - is the next point a valid number?
            double x1 = value(i + 1);
            if ( ! isfinite(x1) ) {
                return shyft::nan;
            }
            // - interpolate!
            utctime t1 = ts->time(i + 1);
            double a = (x1 - x0) / to_seconds(t1 - t0);
            double b = x0 - a * to_seconds(t0);

            return a * to_seconds(t) + b; // otherwise linear interpolation
        }
        vector<double> qac_ts::values() const {
            auto r=ts->values();

            // construct repeating values cache if neccessary
            if ( this->p.repeating_qa_enabled ) {
                if ( this->_repeating_cache.size() != this->size() ) {  // simple check to determine if repeating values are cached
                    this->_scan_repeating_values(r);
                }
            }

            // compute values for this ts
            for (size_t i = 0; i < r.size(); ++i) {
                if ( ! p.is_ok_quality(i, r[i], *this) )
                    r[i]=_fill_value(i);
            }

            return r;
        }

        // ================================================================================

        double inside_ts::value(size_t i) const {
            return p.inside_value(ts->value(i));
        }

        double inside_ts::value_at(utctime t) const {
            size_t i = index_of(t);
            if (i == string::npos)
                return shyft::nan;
            return value(i);
        }

        vector<double> inside_ts::values() const {
            const size_t n{ size() };
            vector<double> r; r.reserve(n);
            for (size_t i = 0; i < n; ++i)
                r.emplace_back(value(i));
            return r;
        }

        double decode_ts::value(size_t i) const {
            return p.decode(ts->value(i));
        }

        double decode_ts::value_at(utctime t) const {
            size_t i = index_of(t);
            if (i == string::npos)
                return shyft::nan;
            return value(i);
        }

        vector<double> decode_ts::values() const {
            const size_t n{ size() };
            vector<double> r; r.reserve(n);
            for (size_t i = 0; i < n; ++i)
                r.emplace_back(value(i));
            return r;
        }
        }
    }
}

