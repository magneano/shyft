/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#include "boostpython_pch.h"
#include "core/utctime_utilities.h"
#include "core/time_axis.h"
#include "core/time_series.h"
#include "core/predictions.h"
#include "api/api.h"
#include "core/time_series_dd.h"
#include "py_convertible.h"


namespace expose {
    using namespace shyft;
    using namespace shyft::core;
    using namespace boost::python;
    using namespace std;
	using shyft::time_series::dd::ats_vector;
	using shyft::time_series::dd::quantile_map_forecast;
	using shyft::time_series::dd::rts_t;
	using namespace shyft::time_series::dd;
    namespace py = boost::python;
    
    ats_vector quantile_map_forecast_5(vector<ats_vector> const & forecast_set, vector<double> const& set_weights, ats_vector const& historical_data, shyft::time_series::dd::gta_t const&time_axis, utctime interpolation_start ) {
        return quantile_map_forecast(forecast_set, set_weights, historical_data, time_axis, interpolation_start);
    }
    ats_vector quantile_map_forecast_6(vector<ats_vector> const & forecast_set, vector<double> const& set_weights, ats_vector const& historical_data, shyft::time_series::dd::gta_t const&time_axis, utctime interpolation_start, utctime interpolation_end) {
        return quantile_map_forecast(forecast_set, set_weights, historical_data, time_axis, interpolation_start, interpolation_end);
    }
    ats_vector quantile_map_forecast_7(vector<ats_vector> const & forecast_set, vector<double> const& set_weights, ats_vector const& historical_data, shyft::time_series::dd::gta_t const&time_axis, utctime interpolation_start, utctime interpolation_end, bool interpolated_quantiles) {
        return quantile_map_forecast(forecast_set, set_weights, historical_data, time_axis, interpolation_start, interpolation_end, interpolated_quantiles);
    }
	static string nice_str(const gta_t & ta) {
		char s[100]; s[0] = 0;
		switch (ta.gt) {
		case gta_t::generic_type::FIXED:sprintf(s, "TaF[%s,%g,%zd]", calendar().to_string(ta.f.t).c_str(), to_seconds(ta.f.dt), ta.f.n); break;
		case gta_t::generic_dt::CALENDAR:sprintf(s, "TaC[%s,%s,%g,%zd]", ta.c.cal->tz_info->name().c_str(), ta.c.cal->to_string(ta.c.t).c_str(), to_seconds(ta.c.dt), ta.c.n); break;
		case gta_t::generic_dt::POINT:sprintf(s, "TaP[%s,%zd]", ta.p.total_period().to_string().c_str(), ta.p.size()); break;
		}
		return s;
	}

	static string nice_str(const string& lhs,iop_t op,const string &rhs) {
        switch (op) {
        case iop_t::OP_ADD:return "("+ lhs + " + " + rhs + ")";
        case iop_t::OP_SUB:return "(" + lhs + " - " + rhs + ")";
        case iop_t::OP_DIV:return "(" + lhs + "/" + rhs + ")";
        case iop_t::OP_MUL:return "(" + lhs + "*" + rhs + ")";
        case iop_t::OP_MAX:return "max(" + lhs + ", " + rhs + ")";
        case iop_t::OP_MIN:return "min(" + lhs + ", " + rhs + ")";
        case iop_t::OP_POW:return "pow(" + lhs + ", " + rhs + ")";
        case iop_t::OP_LOG:return "log(" + lhs + ")";
        case iop_t::OP_NONE:break;// just fall to exception
        }
        return "unsupported_op(" + lhs + "," + rhs + ")";
	}
	static string nice_str(derivative_method dm) {
        switch(dm) {
            case derivative_method::default_diff:return "default";
            case derivative_method::forward_diff:return "forward";
            case derivative_method::backward_diff:return "backward";
            case derivative_method::center_diff:return "center";            
        }
        return "unknown_method";
    }
	static string nice_str(const apoint_ts&ats);
	static string nice_str(double x) { return to_string(x); }
	static string nice_str(const shared_ptr<gpoint_ts>& g) {return g? "Ts{" + nice_str(g->rep.ta) + "v[..]}":string("null");}
	static string nice_str(const shared_ptr<aref_ts>&r ) {return r->id + (r->rep?"(" +nice_str(r->rep)+ ")":string(""));}
	static string nice_str(const shared_ptr<abin_op_ts> &b) { return nice_str(nice_str(b->lhs), b->op, nice_str(b->rhs)); }
	static string nice_str(const shared_ptr<abin_op_ts_scalar> &b) { return nice_str(nice_str(b->lhs), b->op, nice_str(b->rhs)); }
	static string nice_str(const shared_ptr<abin_op_scalar_ts> &b) { return nice_str(nice_str(b->lhs), b->op, nice_str(b->rhs)); }
	static string nice_str(const shared_ptr<abs_ts> &b) { return "abs(" +nice_str(apoint_ts(b->ts)) + ")"; }
	static string nice_str(const shared_ptr<time_series::dd::average_ts>&b) { return "average(" + nice_str(apoint_ts(b->ts)) + "," + nice_str(b->ta) + ")"; }
	static string nice_str(const shared_ptr<time_series::dd::integral_ts>&b) { return "integral(" + nice_str(apoint_ts(b->ts)) + "," + nice_str(b->ta) + ")"; }
	static string nice_str(const shared_ptr<time_series::dd::accumulate_ts>&b) { return "accumulate(" + nice_str(apoint_ts(b->ts)) + "," + nice_str(b->ta) + ")"; }
	static string nice_str(const shared_ptr<time_series::dd::time_shift_ts>&b) { return "time_shift(" + nice_str(apoint_ts(b->ts)) + "," + to_string(to_seconds(b->dt)) + ")"; }
	static string nice_str(const shared_ptr<time_series::dd::periodic_ts>&b) { return "periodic_ts("+nice_str(b->ts.ta) + ")"; }
	static string nice_str(const shared_ptr<time_series::dd::convolve_w_ts>&b) { return "convolve_w_ts(" + nice_str(b->ts_impl.ts) + ",..)"; }
	static string nice_str(const shared_ptr<time_series::dd::extend_ts>&b) { return "extend_ts(" + nice_str(b->lhs)+","+nice_str(b->rhs)+",..)"; }
	static string nice_str(const shared_ptr<time_series::dd::use_time_axis_from_ts>&b) { return "(" + nice_str(b->lhs)+".use_time_axis_from("+nice_str(b->rhs)+",..))"; }
    static string nice_str(const shared_ptr<time_series::dd::rating_curve_ts>&b) { return "rating_curve_ts(" + nice_str(b->ts.level_ts) + ",..)"; }
    static string nice_str(const shared_ptr<time_series::dd::ice_packing_ts>&b) { return "ice_packing_ts(" + nice_str(b->ts.temp_ts) + ",..)"; }
    static string nice_str(const shared_ptr<time_series::dd::ice_packing_recession_ts>&b) { return "ice_packing_recession_ts(" + nice_str(b->flow_ts) + "," + nice_str(b->ice_packing_ts) + ",..)"; }
	static string nice_str(const shared_ptr<time_series::dd::krls_interpolation_ts>&b) { return "krls(" + nice_str(b->ts) + ",..)"; }
	static string nice_str(const shared_ptr<time_series::dd::qac_ts>&b) { return "qac_ts(" + nice_str(apoint_ts(b->ts)) + ", "+nice_str(apoint_ts(b->cts))+"..)"; }
	static string nice_str(const shared_ptr<time_series::dd::inside_ts>&b) { return "inside_ts(" + nice_str(apoint_ts(b->ts)) + ", "+to_string(b->p.min_x)+", "+to_string(b->p.max_x)+", ..)"; }
	static string nice_str(const shared_ptr<time_series::dd::decode_ts>&b) { return "decode_ts(" + nice_str(apoint_ts(b->ts)) + ",start_bit="+to_string(b->p.start_bit)+",n_bits="+to_string(b->p.n_bits())+")"; }
	static string nice_str(const shared_ptr<time_series::dd::bucket_ts>&b) { return "bucket_ts(" + nice_str(apoint_ts(b->ts)) + ",start_hour_utc="+to_string(static_cast<int>(to_seconds64(b->p.hour_offset)/3600)) + ",bucket_emptying_limit="+to_string(b->p.bucket_empty_limit)+")"; }
	static string nice_str(const shared_ptr<time_series::dd::derivative_ts>&b) { return "derivative(" + nice_str(apoint_ts(b->ts)) + ",dm="+nice_str(b->dm)+")"; }


	static string nice_str(const apoint_ts&ats) {
		if (!ats.ts)
			return "null";
		if (const auto& g = dynamic_pointer_cast<gpoint_ts>(ats.ts))return nice_str(g);
		if (const auto& r = dynamic_pointer_cast<aref_ts>(ats.ts)) return nice_str(r);
		if (const auto& b = dynamic_pointer_cast<abin_op_ts>(ats.ts)) return nice_str(b);
		if (const auto& b = dynamic_pointer_cast<abin_op_ts_scalar>(ats.ts)) return nice_str(b);
		if (const auto& b = dynamic_pointer_cast<abin_op_scalar_ts>(ats.ts)) return nice_str(b);
		if (const auto& b = dynamic_pointer_cast<abs_ts>(ats.ts)) return nice_str(b);
		if (const auto& b = dynamic_pointer_cast<time_series::dd::average_ts>(ats.ts)) return nice_str(b);
		if (const auto& b = dynamic_pointer_cast<time_series::dd::integral_ts>(ats.ts)) return nice_str(b);
		if (const auto& b = dynamic_pointer_cast<time_series::dd::accumulate_ts>(ats.ts)) return nice_str(b);
		if (const auto& b = dynamic_pointer_cast<time_series::dd::time_shift_ts>(ats.ts)) return nice_str(b);
		if (const auto& b = dynamic_pointer_cast<time_series::dd::periodic_ts>(ats.ts)) return nice_str(b);
		if (const auto& b = dynamic_pointer_cast<time_series::dd::convolve_w_ts>(ats.ts)) return nice_str(b);
		if (const auto& b = dynamic_pointer_cast<time_series::dd::extend_ts>(ats.ts)) return nice_str(b);
        if (const auto& b = dynamic_pointer_cast<time_series::dd::rating_curve_ts>(ats.ts)) return nice_str(b);
        if (const auto& b = dynamic_pointer_cast<time_series::dd::ice_packing_ts>(ats.ts)) return nice_str(b);
        if (const auto& b = dynamic_pointer_cast<time_series::dd::ice_packing_recession_ts>(ats.ts)) return nice_str(b);
		if (const auto& b = dynamic_pointer_cast<time_series::dd::krls_interpolation_ts>(ats.ts)) return nice_str(b);
		if (const auto& b = dynamic_pointer_cast<time_series::dd::qac_ts>(ats.ts)) return nice_str(b);
        if (const auto& b = dynamic_pointer_cast<time_series::dd::inside_ts>(ats.ts)) return nice_str(b);
        if (const auto& b = dynamic_pointer_cast<time_series::dd::decode_ts>(ats.ts)) return nice_str(b);
        if (const auto& b = dynamic_pointer_cast<time_series::dd::use_time_axis_from_ts>(ats.ts)) return nice_str(b);
        if (const auto& b = dynamic_pointer_cast<time_series::dd::bucket_ts>(ats.ts)) return nice_str(b);
        if (const auto& b = dynamic_pointer_cast<time_series::dd::derivative_ts>(ats.ts)) return nice_str(b);
        
		return "not_yet_stringified_ts";
	}
    static string ts_stringify(const apoint_ts&ats) { return nice_str(ats); }
	static void expose_core_ts_vector() {
		typedef vector<pts_t> core_ts_vector;
		class_<core_ts_vector>("CoreTsVector",
			doc_intro("A raw vector of core time-series.")
			doc_intro("")
			)
			.def(vector_indexing_suite<core_ts_vector>());
	}

    /** ats_vector python expose helper for constructor variants so that
     * it's easy for to use. We can't use py_convertible since
     * python get confused by 3*tsv etc. (trying to convert 3 -> tsvector)
     */
	struct ats_vector_ext {
        static ats_vector *create_default() {return new ats_vector{};}
        static ats_vector *create_from_clone(const ats_vector& c) {return new ats_vector(c);}
        static ats_vector *create_from_ts_list(py::list tsl) {
            size_t n = py::len(tsl);
            if(n==0) return new ats_vector{};
            auto r= new ats_vector();
            r->reserve(n);
            for(size_t i=0;i<n;++i) {
                py::object ts= tsl[i];
                py::extract<apoint_ts> xts(ts);
                if(xts.check()) {
                    r->push_back(xts());
                } else {
                    throw runtime_error("Failed to convert "+ to_string(i)+" element to TimeSeries");
                }
            }
            return r;
        }
    };

    static void expose_ats_vector() {
        using namespace shyft::api;
        using shyft::time_series::dd::pow;
        typedef ats_vector(ats_vector::*m_double)(double)const;
        typedef ats_vector(ats_vector::*m_ts)(apoint_ts const&)const;
        typedef ats_vector(ats_vector::*m_tsv) (ats_vector const&)const;
        typedef ats_vector(ats_vector::*m_na) () const;

        class_<ats_vector>("TsVector",
                doc_intro("A vector of time-series that supports ts-math operations.")
                doc_intro("")
                doc_intro("Like:")
                doc_intro("  number bin_op ts_vector -> ts_vector")
                doc_intro("  ts_vector bin_op ts_vector -> ts_vector ")
                doc_intro("  ts bin_op ts_vector -> ts_vector")
                doc_intro("  where bin_op is any of (*,/,+,-)")
                doc_intro("")
                doc_intro("In addition, .average(..),.integral(..),.accumulate(..),.time_shift(..), .percentiles(..)")
                doc_intro("  is also supported")
                doc_intro("")
                doc_intro("All operation returns a *new* ts-vector, containing the resulting expressions"),
                no_init
            )
            .def(vector_indexing_suite<ats_vector>())
            .def("__init__",make_constructor(&ats_vector_ext::create_default,default_call_policies()),
                doc_intro("Create an empty TsVector")
            )
            .def("__init__",make_constructor(&ats_vector_ext::create_from_clone,default_call_policies(),(py::arg("cloneme"))),
                doc_intro("Create a shallow clone of  the TsVector")
                doc_parameters()
                doc_parameter("cloneme","TsVector","The TsVector to be cloned")
            )
            .def("__init__",make_constructor(&ats_vector_ext::create_from_ts_list,default_call_policies(),(py::arg("ts_list"))),
                doc_intro("Create a TsVector from a python list of TimeSeries")
                doc_parameters()
                doc_parameter("ts_list","List[TimeSeries]","A list of TimeSeries")
            )
            .def("values_at",&ats_vector::values_at_time,args("t"),
                 doc_intro("Computes the value at specified time t for all time-series")
                 doc_parameters()
                 doc_parameter("t","utctime","seconds since epoch 1970 UTC")
            )
			.def("values_at", &ats_vector::values_at_time_i, args("t"),
				doc_intro("Computes the value at specified time t for all time-series")
				doc_parameters()
				doc_parameter("t", "int", "seconds since epoch 1970 UTC")
			)
            .def("percentiles",&ats_vector::percentiles,args("time_axis","percentiles"),
                doc_intro("Calculate the percentiles, NIST R7, excel,R definition, of the timeseries")
                doc_intro("over the specified time-axis.")
                doc_intro("The time-series point_fx interpretation is used when performing")
                doc_intro("the true-average over the time_axis periods.")
                doc_parameters()
                doc_parameter("percentiles","IntVector","A list of numbers,[ 0, 25,50,-1,75,100] will return 6 time-series,\n -1 -> arithmetic average\n -1000 -> min extreme value\n +1000 max extreme value")
                doc_parameter("time_axis","TimeAxis","The time-axis used when applying true-average to the time-series")
                doc_returns("calculated_percentiles","TsVector","Time-series list with evaluated percentile results, same length as input")
            )
            .def("percentiles",&ats_vector::percentiles_f,args("time_axis","percentiles"),
                doc_intro("Calculate the percentiles, NIST R7, excel,R definition, of the timeseries")
                doc_intro("over the specified time-axis.")
                doc_intro("The time-series point_fx interpretation is used when performing")
                doc_intro("the true-average over the time_axis periods.")
                doc_parameters()
                doc_parameter("percentiles","IntVector","A list of numbers,[ 0, 25,50,-1,75,100] will return 6 time-series,\n -1 -> arithmetic average\n -1000 -> min extreme value\n +1000 max extreme value")
                doc_parameter("time_axis","TimeAxisFixedDeltaT","The time-axis used when applying true-average to the time-series")
                doc_returns("calculated_percentiles","TsVector","Time-series list with evaluated percentile results, same length as input")
            )
            .def("slice",&ats_vector::slice,args("indexes"),
                 doc_intro("returns a slice of self, specified by indexes")
                 doc_parameters()
                 doc_parameter("indexes","IntVector","the indicies to pick out from self, if indexes is empty, then all is returned")
                 doc_returns("slice","TsVector","a new TsVector, with content according to indexes specified")
            )
            .def("abs", &ats_vector::abs,
                doc_intro("create a new ts-vector, with all members equal to abs(self")
                doc_returns("tsv", "TsVector", "a new TsVector expression, that will provide the abs-values of self.values")
            )

            .def("average", &ats_vector::average, args("ta"),
                doc_intro("create a new vector of ts that is the true average of self")
                doc_intro("over the specified time-axis ta.")
                doc_parameters()
                doc_parameter("ta","TimeAxis","time-axis that specifies the periods where true-average is applied")
                doc_returns("tsv","TsVector","a new time-series expression, that will provide the true-average when requested")
                doc_notes()
                doc_note("the self point interpretation policy is used when calculating the true average")
			)
            .def("integral", &ats_vector::integral, args("ta"),
                doc_intro("create a new vector of ts that is the true integral of self")
                doc_intro("over the specified time-axis ta.")
                doc_intro(" defined as integral of the non-nan part of each time-axis interval")
                doc_parameters()
                doc_parameter("ta", "TimeAxis", "time-axis that specifies the periods where true-integral is applied")
                doc_returns("tsv", "TsVector", "a new time-series expression, that will provide the true-integral when requested")
                doc_notes()
                doc_note("the self point interpretation policy is used when calculating the true average")
            )
            .def("accumulate", &ats_vector::accumulate, args("ta"),
                doc_intro("create a new  vector of ts where each i'th element is the ")
                doc_intro("    integral f(t) *dt, from t0..ti,")
                doc_intro("given the specified time-axis ta")
                doc_parameters()
                doc_parameter("ta","TimeAxis","time-axis that specifies the periods where accumulated integral is applied")
                doc_returns("tsv","TsVector","a new time-series expression, that will provide the accumulated values when requested")
                doc_notes()
                doc_note("the self point interpretation policy is used when calculating the accumulated values")
            )
            .def("derivative", &ats_vector::derivative, (py::arg("self"),py::arg("method")=derivative_method::default_diff),
                doc_intro("create a new  vector of ts where each i'th element is the ")
                doc_intro("    derivative of f(t)")
                doc_parameters()
                doc_parameter("method","derivative_method","what derivative_method variant to use")
                doc_returns("tsv","TsVector","where each member is the derivative of the source")
            )

			.def("time_shift", &ats_vector::time_shift,args("delta_t"),
				doc_intro("create a new vector of ts that is a the time-shift'ed  version of self")
				doc_parameters()
                doc_parameter("delta_t","int","number of seconds to time-shift, positive values moves forward")
				doc_returns("tsv","TsVector",	"a new time-series, that appears as time-shifted version of self")
			)
            .def("use_time_axis_from",&ats_vector::use_time_axis_from,(py::arg("self"),py::arg("other")),
                doc_intro("Create a new ts-vector applying .use_time_axis_from on each member")
                doc_see_also("TimeSeries.use_time_axis_from")
                doc_intro("")
				doc_parameters()
                doc_parameter("other","TimeSeries","time-series that provides the wanted time-axis")
				doc_returns("tsv","TsVector",	"time-series vector, where each element have time-axis from other")
			) 
            .def("extend_ts", &ats_vector::extend_ts, (py::arg("ts"), py::arg("split_policy") = extend_ts_split_policy::EPS_LHS_LAST, py::arg("fill_policy") = extend_ts_fill_policy::EPF_NAN, py::arg("split_at") = utctime(seconds(0)), py::arg("fill_value") = shyft::nan),
                doc_intro("create a new ats_vector where all time-series are extended by ts")
                doc_parameters()
                doc_parameter("ts", "TimeSeries", "time-series to extend each time-series in self with")
                doc_parameter("split_policy", "extend_ts_split_policy", "policy determining where to split between self and ts")
                doc_parameter("fill_policy", "extend_ts_fill_policy", "policy determining how to fill any gap between self and ts")
                doc_parameter("split_at", "utctime", "time at which to split if split_policy == EPS_VALUE")
                doc_parameter("fill_value", "float", "value to fill any gap with if fill_policy == EPF_FILL")
                doc_returns("new_ts_vec" ,"TsVector", "a new time-series vector where all time-series in self have been extended by ts")
            )
            .def("extend_ts", &ats_vector::extend_vec, (py::arg("ts"), py::arg("split_policy") = extend_ts_split_policy::EPS_LHS_LAST, py::arg("fill_policy") = extend_ts_fill_policy::EPF_NAN, py::arg("split_at") = utctime(seconds(0)), py::arg("fill_value") = shyft::nan),
                doc_intro("create a new ats_vector where all ts' are extended by the matching ts from ts_vec")
                doc_parameters()
                doc_parameter("ts_vec", "TsVector", "time-series vector to extend time-series in self with")
                doc_parameter("split_policy", "extend_ts_split_policy", "policy determining where to split between self and ts")
                doc_parameter("fill_policy", "extend_ts_fill_policy", "policy determining how to fill any gap between self and ts")
                doc_parameter("split_at", "utctime", "time at which to split if split_policy == EPS_VALUE")
                doc_parameter("fill_value", "float", "value to fill any gap with if fill_policy == EPF_FILL")
                doc_returns("new_ts_vec" ,"TsVector", "a new time-series vector where all time-series in self have been extended by the corresponding time-series in ts_vec")
            )
            .def("min",(m_double)&ats_vector::min,(py::arg("self"),py::arg("number")),"returns min of vector and a number")
            .def("min", (m_ts)&ats_vector::min, (py::arg("self"),py::arg("ts")), "returns min of ts-vector and a ts")
            .def("min", (m_tsv)&ats_vector::min, args("tsv"), "returns min of ts-vector and another ts-vector")
            .def("max", (m_double)&ats_vector::max, (py::arg("self"),py::arg("number")), "returns max of vector and a number")
            .def("max", (m_ts)&ats_vector::max, (py::arg("self"),py::arg("ts")), "returns max of ts-vector and a ts")
            .def("max", (m_tsv)&ats_vector::max, (py::arg("self"),py::arg("tsv")), "returns max of ts-vector and another ts-vector")
            .def("pow", (m_double)&ats_vector::pow, (py::arg("self"),py::arg("number")), "returns TsVector pow(self,number)")
            .def("pow", (m_ts)&ats_vector::pow, (py::arg("self"),py::arg("ts")), "returns TsVector pow(self,ts)")
            .def("pow", (m_tsv)&ats_vector::pow, (py::arg("self"), py::arg("tsv")), "returns TsVector pow(self,tsv)")
            .def("log", (m_na)&ats_vector::log, (py::arg("self")), "returns TsVector log(self)")
            .def("forecast_merge",&ats_vector::forecast_merge,args("lead_time","fc_interval"),
                 doc_intro("merge the forecasts in this vector into a time-series that is constructed")
                 doc_intro("taking a slice of length fc_interval starting lead_time into each of the forecasts")
                 doc_intro("of this time-series vector.")
                 doc_intro("The content of the vector should be ordered in forecast-time, each entry at least")
                 doc_intro("fc_interval separated from the previous.")
                 doc_intro("If there is missing forecasts (larger than fc_interval between two forecasts) this is")
                 doc_intro("automagically repaired using extended slices from the existing forecasts")
                 doc_parameters()
                 doc_parameter("lead_time","int","start slice number of seconds from t0 of each forecast")
                 doc_parameter("fc_interval","int","length of each slice in seconds, and thus also gives the forecast-interval separation")
                 doc_returns("merged time-series","TimeSeries","A merged forecast time-series")
                 )
             .def("nash_sutcliffe",&ats_vector::nash_sutcliffe,args("observation_ts","lead_time","delta_t","n"),
                doc_intro("Computes the nash-sutcliffe (wiki nash-sutcliffe) criteria between the")
                doc_intro("observation_ts over the slice of each time-series in the vector.")
                doc_intro("The slice for each ts is specified by the lead_time, delta_t and n")
                doc_intro("parameters. The function is provided to ease evaluation of forecast performance")
                doc_intro("for different lead-time periods into each forecast.")
                doc_intro("The returned value range is 1.0 for perfect match -oo for no match, or nan if observations is constant or data missing")
                doc_parameters()
                doc_parameter("observation_ts","TimeSeries","the observation time-series")
                doc_parameter("lead_time","int","number of seconds lead-time offset from each ts .time(0)")
                doc_parameter("delta_t","int","delta-time seconds to average as basis for n.s. simulation and observation values")
                doc_parameter("n","int","number of time-steps of length delta_t to slice out of each forecast/simulation ts")
                doc_returns("nash-sutcliffe value","double","the nash-sutcliffe criteria evaluated over all time-series in the TsVector for the specified lead-time, delta_t and number of elements")
                  doc_notes()
                  doc_see_also("nash_sutcliffe_goal_function")
             )
             .def("average_slice",&ats_vector::average_slice,args("lead_time","delta_t","n"),
                doc_intro("Returns a ts-vector with the average time-series of the specified slice")
                doc_intro("The slice for each ts is specified by the lead_time, delta_t and n")
                doc_intro("parameters. ")
                doc_parameters()
                doc_parameter("lead_time","int","number of seconds lead-time offset from each ts .time(0)")
                doc_parameter("delta_t","int","delta-time seconds to average as basis for n.s. simulation and observation values")
                doc_parameter("n","int","number of time-steps of length delta_t to slice out of each forecast/simulation ts")
                doc_returns("ts_vector_sliced","TsVector","a ts-vector with average ts of each slice specified.")
                  doc_notes()
                  doc_see_also("nash_sutcliffe,forecast_merge")
             )
             .def("inside", &ats_vector::inside,
                (py::arg("self"),py::arg("min_v"),py::arg("max_v"),py::arg("nan_v")=shyft::nan,py::arg("inside_v")=1.0,py::arg("outside_v")=0.0),
                doc_intro(
                    "Create an inside min-max range ts-vector, that transforms the point-values\n"
                    "that falls into the half open range [min_v .. max_v > to \n"
                    "the value of inside_v(default=1.0), or outside_v(default=0.0),\n"
                    "and if the value considered is nan, then that value is represented as nan_v(default=nan)\n"
                    "You would typically use this function to form a true/false series (inside=true, outside=false)\n"
                )
                doc_parameters()
                doc_parameter("min_v","float","minimum range, values <  min_v are not inside min_v==NaN means no lower limit")
                doc_parameter("max_v","float","maximum range, values >= max_v are not inside. max_v==NaN means no upper limit")
                doc_parameter("nan_v","float","value to return if the value is nan")
                doc_parameter("inside_v","float","value to return if the ts value is inside the specified range")
                doc_parameter("outside_v","float","value to return if the ts value is outside the specified range")
                doc_returns("inside_tsv","TsVector","New TsVector where each element is an evaluated-on-demand inside time-series")
            )
            // defining vector math-operations goes here
            .def(-self)
            .def(self*double())
            .def(double()*self)
            .def(self*self)
            .def(apoint_ts()*self)
            .def(self*apoint_ts())

            .def(self/double())
            .def(double()/self)
            .def(self/self)
            .def(apoint_ts()/self)
            .def(self/apoint_ts())

            .def(self+double())
            .def(double()+self)
            .def(self+self)
            .def(apoint_ts()+self)
            .def(self+apoint_ts())

            .def(self-double())
            .def(double()-self)
            .def(self-self)
            .def(apoint_ts()-self)
            .def(self-apoint_ts())
            .def(operator!(self))
            ;
            // expose min-max functions:
            typedef ats_vector(*f_atsv_double)(ats_vector const &, double b);
            typedef ats_vector(*f_double_atsv)(double b, ats_vector const &a);
            typedef ats_vector(*f_atsv_ats)(ats_vector const &, apoint_ts const& );
            typedef ats_vector(*f_ats_atsv)(apoint_ts const &b, ats_vector const& a);
            typedef ats_vector(*f_atsv_atsv)(ats_vector const &b, ats_vector const &a);

            def("min", (f_ats_atsv)min, args("ts", "ts_vector"), "return minimum of ts and ts_vector");
            def("min", (f_atsv_ats)min, args("ts_vector", "ts"), "return minimum of ts_vector and ts");
            def("min", (f_atsv_double)min, args("ts_vector", "number"), "return minimum of ts_vector and number");
            def("min", (f_double_atsv)min, args("number","ts_vector"), "return minimum of number and ts_vector");
            def("min", (f_atsv_atsv)min, args("a", "b"), "return minimum of ts_vectors a and b (requires equal size!)");
            def("max", (f_ats_atsv)max, args("ts", "ts_vector"), "return max of ts and ts_vector");
            def("max", (f_atsv_ats)max, args("ts_vector", "ts"), "return max of ts_vector and ts");
            def("max", (f_atsv_double)max, args("ts_vector", "number"), "return max of ts_vector and number");
            def("max", (f_double_atsv)max, args("number", "ts_vector"), "return max of number and ts_vector");
            def("max", (f_atsv_atsv)max, args("a", "b"), "return max of ts_vectors a and b (requires equal size!)");

            def("pow", (f_ats_atsv)pow, args("ts", "ts_vector"), "return pow(ts,ts_vector)->TsVector");
            def("pow", (f_atsv_ats)pow, args("ts_vector", "ts"), "return pow(ts_vector,ts)->TsVector");
            def("pow", (f_atsv_double)pow, args("ts_vector", "number"), "return pow(ts_vector,number)->TsVector");
            def("pow", (f_double_atsv)pow, args("number", "ts_vector"), "return pow(number,ts_vector)->TsVector");
            def("pow", (f_atsv_atsv)pow, args("a", "b"), "return pow(a,b)->TsVector, (requires equal size!)");

            // we also need a vector of ats_vector for quantile_map_forecast function
            typedef std::vector<ats_vector> TsVectorSet;
            class_<TsVectorSet>("TsVectorSet",
                doc_intro("A set of TsVector")
                doc_intro("")
                doc_see_also("quantile_map_forecast,TsVector")
            )
            .def(vector_indexing_suite<TsVectorSet>())
            .def(init<TsVectorSet const&>(args("clone_me")))
            ;
            const char* qm_doc =
                doc_intro("Computes the quantile-mapped forecast from the supplied input.")
                doc_intro(" TBD:detailed description with references ")
                doc_parameters()
                doc_parameter("forecast_sets", "TsVectorSet", "forecast sets, each of them a TsVector with n forecasts (might differ in size and length)")
                doc_parameter("set_weights", "DoubleVector", "a weight for each of the forecast set in forecast-sets,correlated by same index)")
                doc_parameter("historical_data", "TsVector", "historical time-series that should cover the requested time-axis")
                doc_parameter("time_axis", "TimeAxis", "the time-axis that the resulting time-series are mapped into")
                doc_parameter("interpolation_start", "int", "time where the historical to forecast interpolation should start, 1970 utc seconds since epoch")
                doc_parameter("interpolation_end", "int", "time where the interpolation should end, if no_utctime, use end of forecast-set")
                doc_parameter("interpolated_quantiles", "bool", "whether the quantile values should be interpolated or assigned the values lower than or equal to the current quantile")
                doc_returns("qm_forecast", "TsVector", "quantile mapped forecast with the requested time-axis")
                ;

            def("quantile_map_forecast",quantile_map_forecast_5,
                args("forecast_sets","set_weights","historical_data","time_axis","interpolation_start"),
                qm_doc
                );
            def("quantile_map_forecast", quantile_map_forecast_6,
                args("forecast_sets", "set_weights", "historical_data", "time_axis", "interpolation_start", "interpolation_end"),
                qm_doc
            );
            def("quantile_map_forecast", quantile_map_forecast_7,
                args("forecast_sets", "set_weights", "historical_data", "time_axis", "interpolation_start", "interpolation_end", "interpolated_quantiles"),
                qm_doc
            );

	}

#define DEF_STD_TS_STUFF() \
            .def("point_interpretation",&pts_t::point_interpretation,(py::arg("self")),"returns the point interpretation policy")\
            .def("set_point_interpretation",&pts_t::set_point_interpretation,(py::arg("self"),py::arg("policy")),"set new policy")\
            .def("value",&pts_t::value,(py::arg("self"),py::arg("i")),"returns the value at the i'th time point")\
            .def("time",&pts_t::time,(py::arg("self"),py::arg("i")),"returns the time at the i'th point")\
            .def("get",&pts_t::get,(py::arg("self"),py::arg("t")),"returns the point(t,v) at time t ")\
            .def("set",&pts_t::set,(py::arg("self"),py::arg("i"),py::arg("v")),"set the i'th value")\
            .def("fill",&pts_t::fill,(py::arg("self"),py::arg("v")),"fill all values with v")\
            .def("scale_by",&pts_t::scale_by,(py::arg("self"),py::arg("v")),"scale all values by the specified factor v")\
            .def("size",&pts_t::size,(py::arg("self")),"returns number of points")\
            .def("index_of", (size_t (pts_t::*)(utctime t) const) &pts_t::index_of,(py::arg("self"),py::arg("t")),"return the index of the intervall that contains t, or npos if not found")\
            .def("total_period",&pts_t::total_period,(py::arg("self")),"returns the total period covered by the time-axis of this time-series")\
            .def("__call__",&pts_t::operator(),(py::arg("self"),py::arg("t")),"return the f(t) value for the time-series")


    template <class TA>
    static void point_ts(const char *ts_type_name,const char *doc) {
        typedef time_series::point_ts<TA> pts_t;
        class_<pts_t,bases<>,shared_ptr<pts_t>,boost::noncopyable>(ts_type_name, doc)
            .def(init<const TA&,const vector<double>&,time_series::ts_point_fx>(
				(py::arg("self"),py::arg("ta"),py::arg("v"),py::arg("policy")),
				doc_intro("constructs a new timeseries from timeaxis, points and policy (how the points are to be interpreted, instant, or average of the interval)")
				)
			)
            .def(init<const TA&,double,time_series::ts_point_fx>(
				(py::arg("self"),py::arg("ta"),py::arg("fill_value"),py::arg("policy")),
				doc_intro("constructs a new timeseries from timeaxis, fill-value and policy")
				)
			)
            DEF_STD_TS_STUFF()
            .def_readonly("v",&pts_t::v,
				doc_intro("the point vector<double>, same as .values, kept around for backward compatibility")
			)
			.def("get_time_axis", &pts_t::time_axis,(py::arg("self")),
				"returns the time-axis", return_internal_reference<>()
			) // have to use func plus init.py fixup due to boost py policy
            ;
    }

    using shyft::core::utctime;
    using shyft::core::seconds;
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(point_ts_overloads     ,shyft::api::TsFactory::create_point_ts,4,5);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(time_point_ts_overloads,shyft::api::TsFactory::create_time_point_ts,3,4);
    extern utctime  x_kwarg_utctime(const py::tuple& args, const py::dict& kwargs,size_t i,const char *kw);
    template<class T>
	static T x_arg(const py::tuple& args, size_t i) {
		if (py::len(args) + 1 < (int)i)
			throw std::runtime_error("missing arg #" + std::to_string(i));
		py::object o = args[i];
		py::extract<T> xtract_arg(o);
		return xtract_arg();
	}

    template<class T>
    static T  x_kwarg(const py::tuple& args, const py::dict& kwargs,size_t i,const char *kw) {
        if (py::len(args)  > (int)i) {
            return x_arg<T>(args,i);
        }
        if(kwargs.has_key(kw)) {
            py::object o = kwargs[kw];
            py::extract<T> xtract_arg(o);
            return xtract_arg();
        }
        throw std::runtime_error("missing kw arg #" + std::string(kw) );
	}
    static vector<utctime>  x_kwarg_utctime_vector(const py::tuple& args, const py::dict& kwargs,size_t i,const char *kw) {
        py::object o;
        if (py::len(args)  > (int)i) {
            o = args[i];
        } else if(kwargs.has_key(kw)) {
            o = kwargs[kw];
        } else {
            throw std::runtime_error("missing kw arg #" + std::string(kw) );
        }

        py::extract<vector<utctime>> xtract_utctime_vector(o);
        if(xtract_utctime_vector.check()) {
            return xtract_utctime_vector();
        }
        // TODO: this does not work, since py do a convertible lookup above, that fails
        py::extract<vector<int64_t>> xtract_int64_vector(o);
        if(xtract_int64_vector.check()) {
            auto v= xtract_int64_vector();
            vector<utctime> r;r.reserve(v.size());
            for(const auto &vi:v)r.emplace_back(seconds(vi));
            return r;
        }
        throw std::runtime_error("Expected UtcTimeVector, or Int64Vector  for kw arg #" + std::string(kw) );
	}

    template<class T>
    static T  x_kwarg_default(const py::tuple& args, const py::dict& kwargs,size_t i,const char *kw,const T& default_value) {
        if (py::len(args)  > (int)i) {
            return x_arg<T>(args,i);
        }
        if(kwargs.has_key(kw)) {
            py::object o = kwargs[kw];
            py::extract<T> xtract_arg(o);
            return xtract_arg();
        }
        return default_value;
	}
    using shyft::core::utcperiod;

    /** helper to ensure we can take any time-rep as input args */
    struct ts_factory_ext {

        static shyft::api::TsFactory x_self(const py::tuple& args) {
			if (py::len(args) == 0)
				throw std::runtime_error("self is null in UtcTime");
			py::object self = args[0];
			py::extract<shyft::api::TsFactory> xtract_self(self);
			return xtract_self();
		}

        static py::object create_point_ts(const py::tuple& args, const py::dict& kwargs) {
            auto self=x_self(args);
            size_t n=x_kwarg<int>(args,kwargs,1,"n");
            auto  t= x_kwarg_utctime(args,kwargs,2,"t");
            auto dt=x_kwarg_utctime(args,kwargs,3,"dt");
            auto v = x_kwarg<vector<double>>(args,kwargs,4,"values");
            auto ct= x_kwarg_default<shyft::time_series::ts_point_fx>(args,kwargs,5,"interpretation",shyft::time_series::POINT_INSTANT_VALUE);
            return py::object(self.create_point_ts(n,t,dt,v,ct));
        }
        static py::object create_time_point_ts(const py::tuple& args, const py::dict& kwargs) {
            auto self=x_self(args);
            auto p=x_kwarg<utcperiod>(args,kwargs,1,"period");
            auto  t= x_kwarg_utctime_vector(args,kwargs,2,"times");
            auto v = x_kwarg<vector<double>>(args,kwargs,3,"values");
            auto ct= x_kwarg_default<shyft::time_series::ts_point_fx>(args,kwargs,4,"interpretation",shyft::time_series::POINT_INSTANT_VALUE);
            return py::object(self.create_time_point_ts(p,t,v,ct));
        }
    };

    static void TsFactory() {
        class_<shyft::api::TsFactory>("TsFactory",
			doc_intro("TsFactory is used in specific contexts, to create point time-series that exposes the ITimeSeriesOfPoint interface, using the internal ts-implementations")
			doc_intro("This class is intended for internal shyft-use, related to calibration/running etc.")
			doc_intro("For geneneral time-series, please use TimeSeries() that have plenty well defined constructors")
			,
			init<>(py::arg("self"))
			)
            .def("create_point_ts",raw_function(ts_factory_ext::create_point_ts,4),
                doc_intro("creates a fixed interval time-series based on input parameters")
                doc_parameters()
                doc_parameter("n","int","number of points in time-series")
                doc_parameter("t","utctime","start of first point, as seconds since epoch")
                doc_parameter("dt","utctime","interval, in seconds")
                doc_parameter("values","DoubleVector","values as DoubleVector")
                doc_parameter("interpretation","ts_point_fx","point interpretation, default POINT_INSTANT_VALUE, other value is POINT_AVERAGE_VALUE ")
                doc_returns("ts","TimeSeries","constructed time-series")
            )
            .def("create_time_point_ts",raw_function(ts_factory_ext::create_time_point_ts,3),
                 //args("period","times","values","interpretation"),"return a point ts from specified arguments")
                doc_intro("creates a variable interval time-series based on input parameters")
                doc_parameters()
                doc_parameter("period","utcperiod"," where .start should be equal to the first point in the supplied time-vector,t, and .end should be the end-time of the last interval")
                doc_parameter("times","UtcTimeVector","start time of each interval")
                doc_parameter("values","DoubleVector","values as DoubleVector, one for each interval")
                doc_parameter("interpretation","ts_point_fx","point interpretation, default POINT_INSTANT_VALUE, other value is POINT_AVERAGE_VALUE ")
                doc_returns("ts","TimeSeries","constructed time-series")
             )
            ;
    }


    static void expose_apoint_ts() {
        using namespace shyft::api;
        using shyft::time_series::dd::qac_parameter;
        typedef apoint_ts pts_t;
        typedef pts_t (pts_t::*self_na_t)() const;  // no arg function type
        typedef pts_t (pts_t::*self_dbl_t)(double) const;  // double arg function type
        typedef pts_t (pts_t::*self_ts_t)(const pts_t &)const;  // ts arg function type
        
        self_dbl_t min_double_f=&pts_t::min;
        self_ts_t  min_ts_f =&pts_t::min;
        self_dbl_t max_double_f=&pts_t::max;
        self_ts_t  max_ts_f =&pts_t::max;
        self_dbl_t pow_double_f=&pts_t::pow;
        self_ts_t  pow_ts_f = &pts_t::pow;
        self_na_t  log_na_f =&pts_t::log;

        typedef ts_bind_info TsBindInfo;
        class_<TsBindInfo>("TsBindInfo",
            doc_intro("TsBindInfo gives information about the time-series and it's binding")
			doc_intro("represented by encoded string reference")
			doc_intro("Given that you have a concrete ts,")
			doc_intro("you can bind that the bind_info.ts")
			doc_intro("using bind_info.ts.bind()")
			doc_intro("see also Timeseries.find_ts_bind_info() and Timeseries.bind()"),
			init<>(py::arg("self"))
            )
            .def_readwrite("id", &ts_bind_info::reference, "a unique id/url that identifies a time-series in a ts-database/file-store/service")
            .def_readwrite("ts", &ts_bind_info::ts,"the ts, provides .bind(another_ts) to set the concrete values")
            ;

        typedef vector<TsBindInfo> TsBindInfoVector;
        class_<TsBindInfoVector>("TsBindInfoVector",
			doc_intro("A vector of TsBindInfo")
			doc_intro("see also TsBindInfo"),
			init<>(py::arg("self"))
			)
            .def(vector_indexing_suite<TsBindInfoVector>())
            ;
        py_api::iterable_converter().from_python<TsBindInfoVector>();

        apoint_ts (apoint_ts::*min_max_check_linear_fill_t)(double ,double ,utctimespan ) const = &apoint_ts::min_max_check_linear_fill;
        apoint_ts  (apoint_ts::*min_max_check_ts_fill_t)(double,double,utctimespan,const apoint_ts& ) const=&apoint_ts::min_max_check_ts_fill;
        apoint_ts  (apoint_ts::*min_max_check_linear_fill_i)(double,double ,int64_t ) const= &apoint_ts::min_max_check_linear_fill;
        apoint_ts  (apoint_ts::*min_max_check_ts_fill_i)(double,double,int64_t ,const apoint_ts&) const =&apoint_ts::min_max_check_ts_fill;


		class_<apoint_ts>("TimeSeries",
                doc_intro("A time-series providing mathematical and statistical operations and functionality.")
                doc_intro("")
                doc_intro("A time-series can be an expression, or a concrete point time-series.")
                doc_intro("All time-series do have a time-axis, values, and a point fx policy.")
                doc_intro("")
                doc_intro("The time-series can provide a value for all the intervals, and the point_fx policy")
                doc_intro("defines how the values should be interpreted:")
                doc_intro("POINT_INSTANT_VALUE:")
                doc_intro("    the point value is valid at the start of the period, linear between points")
                doc_intro("    -extend flat from last point to +oo, nan before first value")
                doc_intro("    typical for state-variables, like water-level, temperature measured at 12:00 etc.")
                doc_intro("POINT_AVERAGE_VALUE:")
                doc_intro("    the point represents an average or constant value over the period")
                doc_intro("    typical for model-input and results, precipitation mm/h, discharge m^3/s")
                doc_intro("")
                doc_intro("Example:")
                doc_intro("import numpy as np")
                doc_intro("from shyft.api import Calendar,deltahours,TimeAxis,TimeSeries,POINT_AVERAGE_VALUE as fx_avg,DoubleVector as dv")
                doc_intro("")
                doc_intro("utc = Calendar()  # ensure easy consistent explicit handling of calendar and time")
                doc_intro("ta = TimeAxis(utc.time(2016, 9, 1, 8, 0, 0), deltahours(1), 10)  # create a time-axis to use")
                doc_intro("a = TimeSeries(ta, dv.from_numpy(np.linspace(0, 10, num=len(ta))), fx_avg)")
                doc_intro("b = TimeSeries(ta, dv.from_numpy(np.linspace(0,  1, num=len(ta))), fx_avg)")
                doc_intro("c = a + b*3.0  # c is now an expression, time-axis is the overlap of a and b, lazy evaluation")
                doc_intro("c_values = c.values.to_numpy()  # compute and extract the values, as numpy array")
                doc_intro("c_evaluated=c.evaluate() # computes the expression, return a new concrete point-ts equal to the expression")
                doc_intro("")
                doc_intro("The TimeSeries functionality includes ")
                doc_intro(" resampling:average,accumulate,time_shift")
                doc_intro(" statistics: min/max,correlation by nash-sutcliffe, kling-gupta")
                doc_intro(" filtering: convolution,average,derivative")
                doc_intro(" quality and correction: min-max limits, replace by linear interpolation or replacement ts")
                doc_intro(" partitioning and percentiles ")
                doc_intro("Please check notebooks, examples and api-tests for usage.")
                doc_see_also("TimeAxis,DoubleVector,Calendar,point_interpretation_policy"),
			    init<>( (py::arg("self")), doc_intro("constructs and empty time-series"))
            )

			.def(init<const time_axis::generic_dt&, const std::vector<double>&, time_series::ts_point_fx >(
				(py::arg("self"),py::arg("ta"), py::arg("values"), py::arg("point_fx")),
				doc_intro("construct a timeseries time-axis ta, corresponding values and point interpretation policy point_fx")
				)
			)
			.def(init<const time_axis::generic_dt&, double, time_series::ts_point_fx >(
				(py::arg("self"),py::arg("ta"), py::arg("fill_value"), py::arg("point_fx")),
				doc_intro("construct a time-series with time-axis ta, specified fill-value, and point interpretation policy point_fx")
				)
			)

			.def(init<const time_axis::fixed_dt&, const std::vector<double>&, time_series::ts_point_fx >(
				(py::arg("self"),py::arg("ta"),py::arg("values"), py::arg("point_fx")),
				doc_intro("construct a timeseries timeaxis ta with corresponding values, and point interpretation policy point_fx")
				)
			)
			.def(init<const time_axis::fixed_dt&, double, time_series::ts_point_fx >(
				(py::arg("self"),py::arg("ta"),py::arg("fill_value"),py::arg("point_fx")),
				doc_intro("construct a timeseries with fixed-delta-t time-axis ta, specified fill-value, and point interpretation policy point_fx")
				)
			)

			.def(init<const time_axis::point_dt&, const std::vector<double>&, time_series::ts_point_fx >(
				(py::arg("self"),py::arg("ta"), py::arg("values"), py::arg("point_fx")),
				doc_intro("construct a time-series with a point-type time-axis ta, corresponding values, and point-interpretation point_fx")
				)
			)
			.def(init<const time_axis::point_dt&, double, time_series::ts_point_fx>(
				(py::arg("self"),py::arg("ta"),py::arg("fill_value"), py::arg("point_fx")),
				doc_intro("construct a time-series with a point-type time-axis ta, specified fill-value, and point-interpretation point_fx")
				)
			)

            .def(init<const rts_t &>(
				(py::arg("self"),py::arg("core_result_ts")),
				doc_intro("construct a time-series from a shyft core time-series, to ease working with core-time-series in user-interface/scripting")
				)
			)
			.def(init<const apoint_ts&>(
				(py::arg("self"),py::arg("clone")),
				doc_intro("creates a shallow copy of the clone time-series")
				)
			)
			.def(init<const vector<double>&, utctimespan, const time_axis::generic_dt&>(
				(py::arg("self"),py::arg("pattern"), py::arg("dt"), py::arg("ta")),
				doc_intro("construct a repeated pattern time-series given a equally spaced dt pattern and a time-axis ta")
				doc_parameters()
				doc_parameter("pattern","DoubleVector","a list of numbers giving the pattern")
				doc_parameter("dt","int","number of seconds between each of the pattern-values")
				doc_parameter("ta","TimeAxis","time-axis that forms the resulting time-series time-axis")
				)
			)
			.def(init<const vector<double>&, utctimespan,utctime, const time_axis::generic_dt&>(
				(py::arg("self"),py::arg("pattern"), py::arg("dt"),py::arg("t0"),py::arg("ta")),
				doc_intro("construct a time-series given a equally spaced dt pattern, starting at t0, and a time-axis ta")
				)
			)
            .def(init<std::string>(
				(py::arg("self"),py::arg("ts_id")),
                doc_intro("constructs a bind-able ts,")
                doc_intro("providing a symbolic possibly unique id that at a later time")
                doc_intro("can be bound, using the .bind(ts) method to concrete values")
                doc_intro("if the ts is used as ts, like size(),.value(),time() before it")
                doc_intro("is bound, then a runtime-exception is raised")
                doc_parameters()
                doc_parameter("ts_id","str","url-like identifier for the time-series,notice that shyft://<container>/<path> is for shyft-internal store")
                )
            )
            .def(init<std::string,const apoint_ts&>(
				(py::arg("self"),py::arg("ts_id"),py::arg("bts")),
                doc_intro("constructs a ready bound ts,")
                doc_intro("providing a symbolic possibly unique id that at a later time")
                doc_intro("can be used to correlate with back-end store\n")
                doc_parameters()
                doc_parameter("ts_id","str","url-type of id, notice that shyft://<container>/<path> is for shyft-internal store")
                doc_parameter("bts","TimeSeries","A time-series, that is either a concrete ts, or an expression that can be evaluated to form a concrete ts")
                )
            )

            .def("ts_id",&apoint_ts::id,(py::arg("self")),
                doc_intro("returns ts_id of symbolic ts, or empty string if not symbolic ts")
                doc_returns("ts_id","str","url-like ts_id as passed to constructor or empty if the ts is not a ts with ts_id")
                doc_see_also("TimeSeries('url://like/id'),TimeSeries('url://like/id',ts_with_values)")
            )
            .def("bucket_to_hourly",&apoint_ts::bucket_to_hourly,(py::arg("self"),py::arg("start_hour_utc"),py::arg("bucket_emptying_limit")),
                doc_intro(
                "Precipitation bucket measurements have a lot of tweaks that needs to be resolved,\n"
                "including negative variations over the day due to faulty temperature-dependent\n"
                "volume/weight sensors attached.\n"
                "\n"
                "A precipitation bucket accumulates precipitation, so the readings should be strictly\n"
                "increasing by time, until the bucket is emptied (full, or as part of maintenance).\n"
                "\n" 
                "The goal for the bucket_to_hourly algorithm is to provide *hourly* precipitation, based on some input signal\n"
                "that usually is hourly(averaging is used if not hourly).\n"
                "\n"
                "The main strategy is to use 24 hour differences (typically at hours in a day where the\n"
                "temperature is low, like early in the morning.), to adjust the hourly volume.\n"
                "\n"
                "Differences in periods of 24hour are distributed on all positive hourly evenets, the \n"
                "negative derivatives are zeroed out, so that the hourly result for each 24 hour\n"
                "is steady increasing, and equal to the difference of the 24hour area.\n"
                "\n"
                "The derivative is then used to compute the hourly precipitation rate in mm/h\n"
                )
                doc_parameters()
                doc_parameter("start_hour_utc","int","valid range [0..24], usually set to early morning(low-stable temperature)")
				doc_parameter("bucket_emptying_limit", "float", "a negative number, range[-oo..0>, limit of when to detect an emptying of a bucket in the unit of the measurements series")
                doc_returns("ts","TimeSeries","a new hourly rate ts, that transforms the accumulated series, compensated for the described defects")
            )
            .def("evaluate",&apoint_ts::evaluate,(py::arg("self")),
                doc_intro("Forces evaluation of the expression, returns a new concrete time-series")
                doc_intro("that is detached from the expression")
                doc_returns("ts","TimeSeries","the evaluated copy of the expression that self represents")
                doc_see_also("TimeSeries.values")
            )
			DEF_STD_TS_STUFF()
			//--
			// expose time_axis sih: would like to use property, but no return value policy, so we use get_ + fixup in init.py

			.def("get_time_axis", &apoint_ts::time_axis,(py::arg("self")), "returns the time-axis", return_internal_reference<>())
			.add_property("values", &apoint_ts::values,"return the values (possibly calculated on the fly)")
			// operators
			.def(self * self)
			.def(double() * self)
			.def(self * double())

			.def(self + self)
			.def(double() + self)
			.def(self + double())

			.def(self / self)
			.def(double() / self)
			.def(self / double())

			.def(self - self)
			.def(double() - self)
			.def(self - double())
			.def(-self)
            .def(operator!(self))
            .def(self==self)
            .def("abs", &apoint_ts::abs,(py::arg("self")),
                doc_intro("create a new ts, abs(self")
                doc_returns("ts", "TimeSeries", "a new time-series expression, that will provide the abs-values of self.values")
            )
			.def("average", &apoint_ts::average, (py::arg("self"),py::arg("ta")),
                doc_intro("create a new ts that is the true average of self")
                doc_intro("over the specified time-axis ta.")
                doc_intro("Notice that same definition as for integral applies; non-nan parts goes into the average")
                doc_parameters()
                doc_parameter("ta","TimeAxis","time-axis that specifies the periods where true-average is applied")
                doc_returns("ts","TimeSeries","a new time-series expression, that will provide the true-average when requested")
                doc_notes()
                doc_note("the self point interpretation policy is used when calculating the true average")
			)
            .def("integral", &apoint_ts::integral,( py::arg("self"),py::arg("ta")),
                doc_intro("create a new ts that is the true integral of self")
                doc_intro("over the specified time-axis ta.")
                doc_intro(" defined as integral of the non-nan part of each time-axis interval")
                doc_parameters()
                doc_parameter("ta", "TimeAxis", "time-axis that specifies the periods where true-integral is applied")
                doc_returns("ts", "TimeSeries", "a new time-series expression, that will provide the true-integral when requested")
                doc_notes()
                doc_note("the self point interpretation policy is used when calculating the true average")
            )
            .def("accumulate", &apoint_ts::accumulate, (py::arg("self"),py::arg("ta")),
                doc_intro("create a new ts where each i'th value is the ")
                doc_intro("    integral f(t) *dt, from t0..ti,")
                doc_intro("given the specified time-axis ta")
                doc_parameters()
                doc_parameter("ta","TimeAxis","time-axis that specifies the periods where accumulated integral is applied")
                doc_returns("ts","TimeSeries","a new time-series expression, that will provide the accumulated values when requested")
                doc_notes()
                doc_note("the self point interpretation policy is used when calculating the accumulated values")
            )
            .def("derivative",&apoint_ts::derivative,(py::arg("self"),py::arg("method")=derivative_method::default_diff),
                 doc_intro("Compute the derivative of the ts, according to the method specified.")
                 doc_intro("For linear(POINT_INSTANT_VALUE), it is always the derivative of the straight line between points,")
                 doc_intro(" - using nan for the interval starting at the last point until end of time-axis.")
                 doc_intro("Default for stair-case(POINT_AVERAGE_VALUE) is the average derivative over each time-step,")
                 doc_intro(" - using 0 as rise for the first/last half of the intervals at the boundaries.")
                 doc_intro("  here you can influence the method used, selecting .forward_diff, .backward_diff")
                 doc_parameters()
                 doc_parameter("method","derivative_method","default value gives center/average derivative .(DEFAULT|FORWARD|BACKWARD|CENTER)")
                 doc_returns("derivative","TimeSeries","The derivative ts")
             )
			.def("time_shift", &apoint_ts::time_shift,(py::arg("self"),py::arg("delta_t")),
				doc_intro("create a new ts that is a the time-shift'ed  version of self")
				doc_parameters()
                doc_parameter("delta_t","int","number of seconds to time-shift, positive values moves forward")
				doc_returns("ts","TimeSeries",	"a new time-series, that appears as time-shifted version of self")
			)
            .def("use_time_axis_from",&apoint_ts::use_time_axis_from,(py::arg("self"),py::arg("other")),
                doc_intro("Create a new ts that have the same values as self, but filtered to the time-axis points from")
                doc_intro("from the other supplied time-series.")
                doc_intro("This function migth be useful for making new time-series, that exactly matches")
                doc_intro("the time-axis of another series.")
                doc_intro("Values of the resulting time-series is like like: ")
                doc_intro(" [self(t) for t in other.time_axis.time_points[:-1]")
                doc_intro("")
				doc_parameters()
                doc_parameter("other","TimeSeries","time-series that provides the wanted time-axis")
				doc_returns("ts","TimeSeries",	"a new time-series, that appears as time-shifted version of self")
			) 
            .def("convolve_w", &apoint_ts::convolve_w, (py::arg("self"),py::arg("weights"), py::arg("policy")),
                doc_intro("create a new ts that is the convolved ts with the given weights list")
                doc_parameters()
                doc_parameter("weights","DoubleVector","the weights profile, use DoubleVector.from_numpy(...) to create these.\n"
                                "\t it's the callers responsibility to ensure the sum of weights are 1.0\n")
                doc_parameter("policy","convolve_policy","(USE_NEAREST|USE_ZERO|USE_NAN + BACKWARD|FORWARD|CENTER)\n"
                "\t Specifies how to handle boundary values\n")
                doc_returns("ts","TimeSeries","a new time-series that is evaluated on request to the convolution of self")
                doc_see_also("ConvolvePolicy")
            )
            .def("krls_interpolation", &apoint_ts::krls_interpolation,
                ( py::arg("self"), py::arg("dt"), py::arg("gamma") = 1.E-3, py::arg("tolerance") = 0.01, py::arg("size") = 1000000u ),
                doc_intro("Compute a new TS that is a krls interpolation of self.")
                doc_intro("")
                doc_intro("The KRLS algorithm is a kernel regression algorithm for aproximating data, the implementation")
                doc_intro("used here is from DLib: http://dlib.net/ml.html#krls")
                doc_intro("The new time-series has the same time-axis as self, and the values vector contain no `nan` entries.")
                doc_intro("")
                doc_intro("If you also want the mean-squared error of the interpolation use get_krls_predictor instead, and")
                doc_intro("use the predictor api to generate a interpolation and a mse time-series.")
                doc_parameters()
                doc_parameter("dt", "float", "The time-step in seconds the underlying predictor is specified for.\n"
                    "    Note that this does not put a limit on time-axes used, but for best results it should be\n"
                    "    approximatly equal to the time-step of time-axes used with the predictor. In addition it\n"
                    "    should not be to long, else you will get poor results. Try to keep the dt less than a day,\n"
                    "    3-8 hours is usually fine." )
                doc_parameter("gamma", "float (optional)", "Determines the width of the radial basis functions for the KRLS algorithm.\n"
                    "    Lower values mean wider basis functions, wider basis functions means faster computation but lower\n"
                    "    accuracy. Note that the tolerance parameter also affects speed and accurcy.\n"
                    "    A large value is around `1E-2`, and a small value depends on the time step. By using values larger\n"
                    "    than `1E-2` the computation will probably take to long. Testing have reveled that `1E-3` works great\n"
                    "    for a time-step of 3 hours, while a gamma of `1E-2` takes a few minutes to compute. Use `1E-4` for a\n"
                    "    fast and tolerably accurate prediction.\n"
                    "    Defaults to `1E-3`" )
                doc_parameter("tolerance", "float (optional)", "The krls training tolerance. Lower values makes the prediction more accurate,\n"
                    "    but slower. This typically have less effect than gamma, but is usefull for tuning. Usually it should be\n"
                    "    either `0.01` or `0.001`.\n"
                    "    Defaults to `0.01`" )
                doc_parameter("size", "int (optional)", "The size of the \"memory\" of the underlying predictor. The default value is\n"
                    "    usually enough. Defaults to `1000000`." )
                doc_intro("")
                doc_intro("Examples\n--------\n")
                doc_intro(">>> import numpy as np")
                doc_intro(">>> import scipy.stats as stat")
                doc_intro(">>> from shyft.api import (")
                doc_intro("...     Calendar, utctime_now, deltahours,")
                doc_intro("...     TimeAxis, TimeSeries")
                doc_intro("... )")
                doc_intro(">>>")
                doc_intro(">>> cal = Calendar()")
                doc_intro(">>> t0 = utctime_now()")
                doc_intro(">>> dt = deltahours(1)")
                doc_intro(">>> n = 365*24  # one year")
                doc_intro(">>>")
                doc_intro(">>> # generate random bell-shaped data")
                doc_intro(">>> norm = stat.norm()")
                doc_intro(">>> data = np.linspace(0, 20, n)")
                doc_intro(">>> data = stat.norm(10).pdf(data) + norm.pdf(np.random.rand(*data.shape))")
                doc_intro(">>> # -----")
                doc_intro(">>> ta = TimeAxis(cal, t0, dt, n)")
                doc_intro(">>> ts = TimeSeries(ta, data)")
                doc_intro(">>>")
                doc_intro(">>> # compute the interpolation")
                doc_intro(">>> ts_ipol = ts.krls_interpolation(deltahours(3))")
                doc_returns("krls_ts", "TimeSeries", "A new time series being the KRLS interpolation of self.")
                doc_see_also("TimeSeries.get_krls_predictor, KrlsRbfPredictor")
            )
            .def("get_krls_predictor", &apoint_ts::get_krls_predictor,
                ( py::arg("self"), py::arg("dt"), py::arg("gamma") = 1.E-3, py::arg("tolerance") = 0.01, py::arg("size") = 1000000u ),
                doc_intro("Get a KRLS predictor trained on this time-series.")
                doc_intro("")
                doc_intro("If you only want a interpolation of self use krls_interpolation instead, this method")
                doc_intro("return the underlying predictor instance that can be used to generate mean-squared error")
                doc_intro("estimates, or can be further trained on more data.")
                doc_notes()
                doc_note("A predictor can only be generated for a bound time-series.")
                doc_parameters()
                doc_parameter("dt", "float", "The time-step in seconds the underlying predictor is specified for.\n"
                    "    Note that this does not put a limit on time-axes used, but for best results it should be\n"
                    "    approximatly equal to the time-step of time-axes used with the predictor. In addition it\n"
                    "    should not be to long, else you will get poor results. Try to keep the dt less than a day,\n"
                    "    3-8 hours is usually fine." )
                doc_parameter("gamma", "float (optional)", "Determines the width of the radial basis functions for the KRLS algorithm.\n"
                    "    Lower values mean wider basis functions, wider basis functions means faster computation but lower\n"
                    "    accuracy. Note that the tolerance parameter also affects speed and accurcy.\n"
                    "    A large value is around `1E-2`, and a small value depends on the time step. By using values larger\n"
                    "    than `1E-2` the computation will probably take to long. Testing have reveled that `1E-3` works great\n"
                    "    for a time-step of 3 hours, while a gamma of `1E-2` takes a few minutes to compute. Use `1E-4` for a\n"
                    "    fast and tolerably accurate prediction.\n"
                    "    Defaults to `1E-3`" )
                doc_parameter("tolerance", "float (optional)", "The krls training tolerance. Lower values makes the prediction more accurate,\n"
                    "    but slower. This typically have less effect than gamma, but is usefull for tuning. Usually it should be\n"
                    "    either `0.01` or `0.001`.\n"
                    "    Defaults to `0.01`" )
                doc_parameter("size", "int (optional)", "The size of the \"memory\" of the underlying predictor. The default value is\n"
                    "    usually enough. Defaults to `1000000`." )
                doc_intro("")
                doc_intro("Examples\n--------\n")
                doc_intro(">>> import numpy as np")
                doc_intro(">>> import scipy.stats as stat")
                doc_intro(">>> from shyft.api import (")
                doc_intro("...     Calendar, utctime_now, deltahours,")
                doc_intro("...     TimeAxis, TimeSeries")
                doc_intro("... )")
                doc_intro(">>>")
                doc_intro(">>> cal = Calendar()")
                doc_intro(">>> t0 = utctime_now()")
                doc_intro(">>> dt = deltahours(1)")
                doc_intro(">>> n = 365*24  # one year")
                doc_intro(">>>")
                doc_intro(">>> # generate random bell-shaped data")
                doc_intro(">>> norm = stat.norm()")
                doc_intro(">>> data = np.linspace(0, 20, n)")
                doc_intro(">>> data = stat.norm(10).pdf(data) + norm.pdf(np.random.rand(*data.shape))")
                doc_intro(">>> # -----")
                doc_intro(">>> ta = TimeAxis(cal, t0, dt, n)")
                doc_intro(">>> ts = TimeSeries(ta, data)")
                doc_intro(">>>")
                doc_intro(">>> # create a predictor")
                doc_intro(">>> pred = ts.get_krls_predictor()")
                doc_intro(">>> total_mse = pred.predictor_mse(ts)  # compute mse relative to ts")
                doc_intro(">>> krls_ts = pred.predict(ta)  # generate a prediction, this is the result from ts.krls_interpolation")
                doc_intro(">>> krls_mse_ts = pred.mse_ts(ts, points=6)  # compute a mse time-series using 6 points around each sample")
                doc_returns("krls_predictor", "KrlsRbfPredictor", "A KRLS predictor pre-trained once on self.")
                doc_see_also("TimeSeries.krls_interpolation, KrlsRbfPredictor")
            )
			.def("rating_curve", &apoint_ts::rating_curve,(py::arg("self"), py::arg("rc_param")),
				doc_intro("Create a new TimeSeries that is computed using a RatingCurveParameter instance.")
				doc_intro("")
                doc_intro("Examples\n--------\n")
                doc_intro("")
				doc_intro(">>> import numpy as np")
				doc_intro(">>> from shyft.api import (")
				doc_intro("...     utctime_now, deltaminutes,")
				doc_intro("...     TimeAxis, TimeSeries,")
				doc_intro("...     RatingCurveFunction, RatingCurveParameters")
				doc_intro("... )")
				doc_intro(">>>")
				doc_intro(">>> # parameters")
				doc_intro(">>> t0 = utctime_now()")
				doc_intro(">>> dt = deltaminutes(30)")
				doc_intro(">>> n = 48*2")
				doc_intro(">>>")
				doc_intro(">>> # make rating function, each with two segments")
				doc_intro(">>> rcf_1 = RatingCurveFunction()")
				doc_intro(">>> rcf_1.add_segment(0, 2, 0, 1)    # add segment from level 0, computing f(h) = 2*(h - 0)**1")
				doc_intro(">>> rcf_1.add_segment(5.3, 1, 1, 1.4)  # add segment from level 5.3, computing f(h) = 1.3*(h - 1)**1.4")
				doc_intro(">>> rcf_2 = RatingCurveFunction()")
				doc_intro(">>> rcf_2.add_segment(0, 1, 1, 1)    # add segment from level 0, computing f(h) = 1*(h - 1)**1")
				doc_intro(">>> rcf_2.add_segment(8.0, 0.5, 0, 2)  # add segment from level 8.0, computing f(h) = 0.5*(h - 0)**2")
				doc_intro(">>>")
				doc_intro(">>> # add rating curves to a parameter pack")
				doc_intro(">>> rcp = RatingCurveParameters()")
				doc_intro(">>> rcp.add_curve(t0, rcf_1)  # rcf_1 is active from t0")
				doc_intro(">>> rcp.add_curve(t0+dt*n//2, rcf_2)  # rcf_2 takes over from t0 + dt*n/2")
				doc_intro(">>>")
				doc_intro(">>> # create a time-axis/-series")
				doc_intro(">>> ta = TimeAxis(t0, dt, n)")
				doc_intro(">>> ts = TimeSeries(ta, np.linspace(0, 12, n))")
				doc_intro(">>> rc_ts = ts.rating_curve(rcp)  # create a new time series computed using the rating curve functions")
				doc_intro(">>>")
				doc_parameters()
				doc_parameter("rc_param", "RatingCurveParameter", "RatingCurveParameter instance.")
				doc_returns("rcts", "TimeSeries", "A new TimeSeries computed using self and rc_param.")
			)
            .def("ice_packing", &apoint_ts::ice_packing, (py::arg("self"), py::arg("ip_params"), py::arg("ipt_policy")),
                doc_intro("Create a binary time-series indicating whether ice-packing is occuring or not.")
                doc_intro("")
                doc_intro("Notes\n-----")
                doc_intro("    `self` is interpreted and assumed to be a temperature time-series.")
                doc_intro("")
                doc_intro("The ice packing detection is based on the mean temperature in a predetermined time")
                doc_intro("window before the time-point of interrest (see `IcePackingParameters.window`.")
                doc_intro("The algorithm determines there to be ice packing when the mean temperature is below")
                doc_intro("a given threshold temperature (see `IcePackingParameters.threshold_temp`).")
                doc_intro("")
                doc_parameters()
                doc_parameter("ip_param", "IcePackingParameters",
                              "Parameter container controlling the ice packing detection.")
                doc_parameter("ipt_policy", "ice_packing_temperature_policy",
                              "Policy flags for determining how to deal with missing temperature values.")
                doc_intro("")
                doc_returns("ice_packing_ts", "TimeSeries", "A time-series indicating wheter ice packing occurs or not")
                doc_intro("")
                doc_intro("Example\n-------")
                doc_intro(">>> import numpy as np")
                doc_intro(">>> from shyft.api import (")
                doc_intro("...     IcePackingParameters, ice_packing_temperature_policy,")
                doc_intro("...     TimeAxis, TimeSeries, point_interpretation_policy, DoubleVector,")
                doc_intro("...     utctime_now, deltahours, deltaminutes,")
                doc_intro("... )")
                doc_intro(">>> ")
                doc_intro(">>> t0 = utctime_now()")
                doc_intro(">>> dt = deltaminutes(15)")
                doc_intro(">>> n = 100")
                doc_intro(">>> ")
                doc_intro(">>> # generate jittery data")
                doc_intro(">>> # - first descending from +5 to -5 then ascending back to +5")
                doc_intro(">>> # - include a NaN hole at the bottom of the V")
                doc_intro(">>> n_ = n if (n//2)*2 == n else n+1  # assure even")
                doc_intro(">>> data = np.concatenate((")
                doc_intro("...     np.linspace(5, -5, n_//2), np.linspace(-5, 5, n_//2)")
                doc_intro("... )) + np.random.uniform(-0.75, 0.75, n_)  # add uniform noise")
                doc_intro(">>> data[n_//2 - 1:n_//2 + 2] = float('nan')  # add some missing data")
                doc_intro(">>> ")
                doc_intro(">>> # create Shyft data structures")
                doc_intro(">>> ta = TimeAxis(t0, dt, n_)")
                doc_intro(">>> temperature_ts = TimeSeries(ta, DoubleVector.from_numpy(data),")
                doc_intro("...                             point_interpretation_policy.POINT_AVERAGE_VALUE)")
                doc_intro(">>> ")
                doc_intro(">>> # do the ice packing detection")
                doc_intro(">>> ip_param = IcePackingParameters(")
                doc_intro("...     threshold_window=deltahours(5),")
                doc_intro("...     threshold_temperature=-1.0)")
                doc_intro(">>> # try all the different temperature policies")
                doc_intro(">>> ice_packing_ts_disallow = temperature_ts.ice_packing(ip_param, ice_packing_temperature_policy.DISALLOW_MISSING)")
                doc_intro(">>> ice_packing_ts_initial = temperature_ts.ice_packing(ip_param, ice_packing_temperature_policy.ALLOW_INITIAL_MISSING)")
                doc_intro(">>> ice_packing_ts_any = temperature_ts.ice_packing(ip_param, ice_packing_temperature_policy.ALLOW_ANY_MISSING)")
                doc_intro(">>> ")
                doc_intro(">>> # plotting")
                doc_intro(">>> from matplotlib import pyplot as plt")
                doc_intro(">>> from shyft.api import time_axis_extract_time_points")
                doc_intro(">>> ")
                doc_intro(">>> # NOTE: The offsets below are added solely to be able to distinguish between the different time-axes")
                doc_intro(">>> ")
                doc_intro(">>> plt.plot(time_axis_extract_time_points(ta)[:-1], temperature_ts.values, label='Temperature')")
                doc_intro(">>> plt.plot(time_axis_extract_time_points(ta)[:-1], ice_packing_ts_disallow.values.to_numpy() + 1,")
                doc_intro("...          label='Ice packing? [DISALLOW_MISSING]')")
                doc_intro(">>> plt.plot(time_axis_extract_time_points(ta)[:-1], ice_packing_ts_initial.values.to_numpy() - 1,")
                doc_intro("...          label='Ice packing? [ALLOW_INITIAL_MISSING]')")
                doc_intro(">>> plt.plot(time_axis_extract_time_points(ta)[:-1], ice_packing_ts_any.values.to_numpy() - 3,")
                doc_intro("...          label='Ice packing? [ALLOW_ANY_MISSING]')")
                doc_intro(">>> plt.legend()")
                doc_intro(">>> plt.show()")
            )
            .def("ice_packing_recession", &apoint_ts::ice_packing_recession, (py::arg("self"), py::arg("ip_ts"), py::arg("ipr_params")),
                doc_intro("Create a new time series where segments are replaced by recession curves.")
                doc_intro("")
                doc_intro("Notes\n-----")
                doc_intro("    The total period (`TimeSeries.total_period`) of `self` needs to be equal to,")
                doc_intro("    or contained in the total period of `ip_ts`.")
                doc_intro("")
                doc_parameters()
                doc_parameter("ip_ts", "TimeSeries",
                    "A binary time-series indicating if ice packing occurring. See `TimeSeries.ice_packing`.")
                doc_parameter("ip_param", "IcePackingParameters",
                    "Parameter container controlling the ice packing recession curve.")
                doc_intro("")
                doc_returns("ice_packing_recession_ts", "TimeSeries",
                    "A time-series where sections in `self` is replaced by recession curves as indicated by `ip_ts`.")
                doc_intro("")
                doc_intro("Example\n-------")
                doc_intro(">>> import numpy as np")
                doc_intro(">>> from shyft.api import (")
                doc_intro("...     IcePackingParameters, IcePackingRecessionParameters, ice_packing_temperature_policy,")
                doc_intro("...     TimeAxis, TimeSeries, point_interpretation_policy, DoubleVector,")
                doc_intro("...     utctime_now, deltahours, deltaminutes,")
                doc_intro("... )")
                doc_intro(">>> ")
                doc_intro(">>> t0 = utctime_now()")
                doc_intro(">>> dt = deltaminutes(15)")
                doc_intro(">>> n = 100")
                doc_intro(">>> ")
                doc_intro(">>> # generate jittery temperature data")
                doc_intro(">>> # - first descending from +5 to -5 then ascending back to +5")
                doc_intro(">>> # - include a NaN hole at the bottom of the V")
                doc_intro(">>> n_ = n if (n//2)*2 == n else n+1  # assure even")
                doc_intro(">>> temperature_data = np.concatenate((")
                doc_intro("...     np.linspace(5, -5, n_//2), np.linspace(-5, 5, n_//2)")
                doc_intro("... )) + np.random.uniform(-0.75, 0.75, n_)  # add uniform noise")
                doc_intro(">>> temperature_data[n_ // 2 - 1:n_ // 2 + 2] = float('nan')  # add some missing data")
                doc_intro(">>> ")
                doc_intro(">>> # create Shyft data structures for temperature")
                doc_intro(">>> ta = TimeAxis(t0, dt, n_)")
                doc_intro(">>> temperature_ts = TimeSeries(ta, DoubleVector.from_numpy(temperature_data),")
                doc_intro("...                             point_interpretation_policy.POINT_AVERAGE_VALUE)")
                doc_intro(">>> ")
                doc_intro(">>> # generate jittery waterflow data")
                doc_intro(">>> # - an upwards curving parabola")
                doc_intro(">>> x0 = ta.total_period().start")
                doc_intro(">>> x1 = ta.total_period().end")
                doc_intro(">>> x = np.linspace(x0, x1, n_)")
                doc_intro(">>> flow_data = -0.0000000015*(x - x0)*(x - x1) + 1 + np.random.uniform(-0.5, 0.5, n_)")
                doc_intro(">>> del x0, x1, x")
                doc_intro(">>> ")
                doc_intro(">>> # create Shyft data structures for temperature")
                doc_intro(">>> flow_ts = TimeSeries(ta, DoubleVector.from_numpy(flow_data),")
                doc_intro("...                      point_interpretation_policy.POINT_AVERAGE_VALUE)")
                doc_intro(">>> ")
                doc_intro(">>> # do the ice packing detection")
                doc_intro(">>> ip_param = IcePackingParameters(")
                doc_intro("...     threshold_window=deltahours(5),")
                doc_intro("...     threshold_temperature=-1.0)")
                doc_intro(">>> # compute the detection time-series")
                doc_intro(">>> # ice_packing_ts = temperature_ts.ice_packing(ip_param, ice_packing_temperature_policy.DISALLOW_MISSING)")
                doc_intro(">>> # ice_packing_ts = temperature_ts.ice_packing(ip_param, ice_packing_temperature_policy.ALLOW_INITIAL_MISSING)")
                doc_intro(">>> ice_packing_ts = temperature_ts.ice_packing(ip_param, ice_packing_temperature_policy.ALLOW_ANY_MISSING)")
                doc_intro(">>> ")
                doc_intro(">>> # setup for the recession curve")
                doc_intro(">>> ipr_param = IcePackingRecessionParameters(")
                doc_intro("...     alpha=0.00009,")
                doc_intro("...     recession_minimum=2.)")
                doc_intro(">>> # compute a recession curve based on the ice packing ts")
                doc_intro(">>> ice_packing_recession_ts_initial = flow_ts.ice_packing_recession(ice_packing_ts, ipr_param)")
                doc_intro(">>> ")
                doc_intro(">>> # plotting")
                doc_intro(">>> from matplotlib import pyplot as plt")
                doc_intro(">>> from shyft.api import time_axis_extract_time_points")
                doc_intro(">>> ")
                doc_intro(">>> plt.plot(time_axis_extract_time_points(ta)[:-1], temperature_ts.values, label='Temperature')")
                doc_intro(">>> plt.plot(time_axis_extract_time_points(ta)[:-1], flow_ts.values, label='Flow')")
                doc_intro(">>> plt.plot(time_axis_extract_time_points(ta)[:-1], ice_packing_ts.values.to_numpy(),")
                doc_intro("...          label='Ice packing?')")
                doc_intro(">>> plt.plot(time_axis_extract_time_points(ta)[:-1], ice_packing_recession_ts_initial.values.to_numpy(),")
                doc_intro("...          label='Recession curve')")
                doc_intro(">>> plt.legend()")
                doc_intro(">>> plt.show()")
            )
            .def("extend", &apoint_ts::extend, (py::arg("self"), py::arg("ts"), py::arg("split_policy") = extend_ts_split_policy::EPS_LHS_LAST, py::arg("fill_policy") = extend_ts_fill_policy::EPF_NAN, py::arg("split_at") = utctime(seconds(0)), py::arg("fill_value") = shyft::nan),
                doc_intro("create a new time-series that is self extended with ts")
                doc_parameters()
                doc_parameter("ts", "TimeSeries", "time-series to extend self with, only values after both the start of self, and split_at is used")
                doc_parameter("split_policy", "extend_split_policy", "policy determining where to split between self and ts")
                doc_parameter("fill_policy", "extend_fill_policy", "policy determining how to fill any gap between self and ts")
                doc_parameter("split_at", "utctime", "time at which to split if split_policy == EPS_VALUE")
                doc_parameter("fill_value", "float", "value to fill any gap with if fill_policy == EPF_FILL")
                doc_returns("extended_ts" ,"TimeSeries", "a new time-series that is the extension of self with ts")
            )
            .def("merge_points",&apoint_ts::merge_points,(py::arg("self"),py::arg("ts")),
                 doc_intro("Given that self is a concrete point-ts(not an expression), or empty ts,")
                 doc_intro("this function modifies the point-set of self, with points, (time,value) from other ts")
                 doc_intro("The result of the merge operation is the distinct set of time-points from self and other ts")
                 doc_intro("where values from other ts overwrites values of self if they happen")
                 doc_intro("to be at the same time-point")
                 doc_parameters()
                 doc_parameter("ts","TimeSeries","time-series to merge the time,value points from")
                 doc_returns("self","TimeSeries","self modified with the merged points from other ts")
            )
            .def("slice",&apoint_ts::slice,(py::arg("self"),py::arg("i0"),py::arg("n")),
                 doc_intro("Given that self is a concrete point-ts(not an expression), or empty ts,")
                 doc_intro("return a new TimeSeries containing the n values starting from index i0.")
                 doc_parameters()
                 doc_parameter("i0", "int", "Index of first element to include in the slice")
                 doc_parameter("n", "int", "Number of elements to include in the slice")
            )
            .def("pow",pow_double_f,(py::arg("self"),py::arg("number")),"create a new ts that contains pow(self,number)")
            .def("pow",pow_ts_f,(py::arg("self"),py::arg("ts_other")),"create a new ts that contains pow(self,ts_other)")
            .def("min",min_double_f,(py::arg("self"),py::arg("number")),"create a new ts that contains the min of self and number for each time-step")
            .def("min",min_ts_f,(py::arg("self"),py::arg("ts_other")),"create a new ts that contains the min of self and ts_other")
            .def("max",max_double_f,(py::arg("self"),py::arg("number")),"create a new ts that contains the max of self and number for each time-step")
            .def("max",max_ts_f,(py::arg("self"),py::arg("ts_other")),"create a new ts that contains the max of self and ts_other")
            .def("log", log_na_f, (py::arg("self")), "create a new ts that contains log(self)")
            .def("min_max_check_linear_fill",min_max_check_linear_fill_i,
                 (py::arg("self"),py::arg("v_min"),py::arg("v_max"),py::arg("dt_max")=shyft::core::max_utctime),
                 doc_intro("Create a min-max range checked ts with linear-fill-values if value is NaN or outside range")
                 doc_parameters()
                 doc_parameter("v_min","float","minimum range, values < v_min are considered NaN. v_min==NaN means no lower limit")
                 doc_parameter("v_max","float","maximum range, values > v_max are considered NaN. v_max==NaN means no upper limit")
                 doc_parameter("dt_max","int","maximum time-range in seconds allowed for interpolating values, default= max_utctime")
                 doc_returns("min_max_check_linear_fill","TimeSeries","Evaluated on demand time-series with NaN, out of range values filled in")
            )
            .def("min_max_check_ts_fill",min_max_check_ts_fill_i,
                 (py::arg("self"),py::arg("v_min"),py::arg("v_max"),py::arg("dt_max"),py::arg("cts")),
                 doc_intro("Create a min-max range checked ts with cts-filled-in-values if value is NaN or outside range")
                 doc_parameters()
                 doc_parameter("v_min","float","minimum range, values < v_min are considered NaN. v_min==NaN means no lower limit")
                 doc_parameter("v_max","float","maximum range, values > v_max are considered NaN. v_max==NaN means no upper limit")
                 doc_parameter("dt_max","int","maximum time-range in seconds allowed for interpolating values")
                 doc_parameter("cts","TimeSeries","time-series that keeps the values to be filled in at points that are NaN or outside min-max-limits")
                 doc_returns("min_max_check_ts_fill","TimeSeries","Evaluated on demand time-series with NaN, out of range values filled in")
            )
            .def("min_max_check_linear_fill",min_max_check_linear_fill_t,
                 (py::arg("self"),py::arg("v_min"),py::arg("v_max"),py::arg("dt_max")=shyft::core::max_utctime),
                 doc_intro("Create a min-max range checked ts with linear-fill-values if value is NaN or outside range")
                 doc_parameters()
                 doc_parameter("v_min","float","minimum range, values < v_min are considered NaN. v_min==NaN means no lower limit")
                 doc_parameter("v_max","float","maximum range, values > v_max are considered NaN. v_max==NaN means no upper limit")
                 doc_parameter("dt_max","int","maximum time-range in seconds allowed for interpolating values, default= max_utctime")
                 doc_returns("min_max_check_linear_fill","TimeSeries","Evaluated on demand time-series with NaN, out of range values filled in")
            )
            .def("min_max_check_ts_fill",min_max_check_ts_fill_t,
                 (py::arg("self"),py::arg("v_min"),py::arg("v_max"),py::arg("dt_max"),py::arg("cts")),
                 doc_intro("Create a min-max range checked ts with cts-filled-in-values if value is NaN or outside range")
                 doc_parameters()
                 doc_parameter("v_min","float","minimum range, values < v_min are considered NaN. v_min==NaN means no lower limit")
                 doc_parameter("v_max","float","maximum range, values > v_max are considered NaN. v_max==NaN means no upper limit")
                 doc_parameter("dt_max","int","maximum time-range in seconds allowed for interpolating values")
                 doc_parameter("cts","TimeSeries","time-series that keeps the values to be filled in at points that are NaN or outside min-max-limits")
                 doc_returns("min_max_check_ts_fill","TimeSeries","Evaluated on demand time-series with NaN, out of range values filled in")
            )
            .def("quality_and_self_correction",&apoint_ts::quality_and_self_correction,
                 (py::arg("self"),py::arg("parameters")),
                 doc_intro("returns a new time-series that applies quality checks accoring to parameters")
                 doc_intro("and fills in values according to rules specified in parameters.")
                 doc_parameters()
                 doc_parameter("parameter","QacParameter","Parameter with rules for quality and corrections")
                 doc_returns("ts","TimeSeries","a new time-series where the values are subject to quality and correction as specified")
             )
            .def("quality_and_ts_correction",&apoint_ts::quality_and_ts_correction,
                 (py::arg("self"),py::arg("parameters"),py::arg("cts")),
                 doc_intro("returns a new time-series that applies quality checks accoring to parameters")
                 doc_intro("and fills in values from the cts, according to rules specified in parameters.")
                 doc_parameters()
                 doc_parameter("parameter","QacParameter","Parameter with rules for quality and corrections")
                 doc_parameter("cts","TimeSeries","is used to fill in correct values, as f(t) for values that fails quality-checks")
                 doc_returns("ts","TimeSeries","a new time-series where the values are subject to quality and correction as specified")
             )
            .def("inside", &apoint_ts::inside,
                (py::arg("self"),py::arg("min_v"),py::arg("max_v"),py::arg("nan_v")=shyft::nan,py::arg("inside_v")=1.0,py::arg("outside_v")=0.0),
                doc_intro(
                    "Create an inside min-max range ts, that transforms the point-values\n"
                    "that falls into the half open range [min_v .. max_v > to \n"
                    "the value of inside_v(default=1.0), or outside_v(default=0.0),\n"
                    "and if the value considered is nan, then that value is represented as nan_v(default=nan)\n"
                    "You would typically use this function to form a true/false series (inside=true, outside=false)\n"
                )
                doc_parameters()
                doc_parameter("min_v","float","minimum range, values <  min_v are not inside min_v==NaN means no lower limit")
                doc_parameter("max_v","float","maximum range, values >= max_v are not inside. max_v==NaN means no upper limit")
                doc_parameter("nan_v","float","value to return if the value is nan")
                doc_parameter("inside_v","float","value to return if the ts value is inside the specified range")
                doc_parameter("outside_v","float","value to return if the ts value is outside the specified range")
                doc_returns("inside_ts","TimeSeries","Evaluated on demand inside time-series")
            )
            .def("decode", &apoint_ts::decode,
                (py::arg("self"),py::arg("start_bit"),py::arg("n_bits")),
                doc_intro(
                    "Create an time-series that decodes the source using provided\n"
                    "specification start_bit and n_bits.\n"
                    "This function can typically be used to decode status-signals from sensors stored as \n"
                    "binary encoded bits, using integer representation\n"
                    "The floating point format allows up to 52 bits to be precisely stored as integer\n"
                    "- thus there are restrictions to start_bit and n_bits accordingly.\n"
                    "Practical sensors quality signals have like 32 bits of status information encoded\n"
                    "If the value in source time-series is \n"
                    " a) negative\n"
                    " b) nan\n"
                    " c) larger than 52 bits\n"
                    "Then nan is returned for those values\n"
                    "\n"
                    "ts.decode(start_bit=1,n_bits=1) will return values [0,1,nan]\n"
                    "similar:\n"
                    "ts.decode(start_bit=1,n_bits=2) will return values [0,1,2,3,nan]\n"
                    "etc..\n"
                )
                doc_parameters()
                doc_parameter("start_bit","int","where in the n-bits integer the value is stored, range[0..51]")
                doc_parameter("n_bits","int","how many bits are encoded, range[0..51], but start_bit +n_bits < 51")
                doc_returns("decode_ts","TimeSeries","Evaluated on demand decoded time-series")
            )

			.def("partition_by",&apoint_ts::partition_by,
                (py::arg("self"),py::arg("calendar"), py::arg("t"), py::arg("partition_interval"), py::arg("n_partitions"), py::arg("common_t0")),
				doc_intro("from a time-series, construct a TsVector of n time-series partitions.")
				doc_intro("The partitions are simply specified by calendar, delta_t(could be symbolic, like YEAR : MONTH:DAY) and n.")
				doc_intro("To make yearly partitions, just pass Calendar.YEAR as partition_interval.")
				doc_intro("The t - parameter set the start - time point in the source-time-series, e.g. like 1930.09.01")
				doc_intro("The common_t0 - parameter set the common start - time of the new partitions, e.g. 2017.09.01")
				doc_intro("")
				doc_intro("The typical usage will be to use this function to partition years into a vector with")
				doc_intro("80 years, where we can do statistics, percentiles to compare and see the different effects of")
				doc_intro("yearly season variations.")
				doc_intro("Note that the function is more general, allowing any periodic partition, like daily, weekly, monthly etc.")
				doc_intro("that allows you to study any pattern or statistics that might be periodic by the partition pattern.")
				doc_parameters()
				doc_parameter("cal","Calendar","The calendar to use, typically utc")
				doc_parameter("t","utctime","specifies where to pick the first partition")
				doc_parameter("partition_interval","utctimespan","the length of each partition, Calendar.YEAR,Calendar.DAY etc.")
				doc_parameter("n_partitions","int","number of partitions")
				doc_parameter("common_t0","utctime","specifies the time to correlate all the partitions")
				doc_returns("ts-partitions","TsVector","with length n_partitions, each ts is time-shifted to common_t0 expressions")
                doc_see_also("time_shift,average,TsVector")
				)
            .def("bind",&apoint_ts::bind,(py::arg("self"),py::arg("bts")),
                doc_intro("given that this ts,self, is a bind-able ts (aref_ts)")
                doc_intro("and that bts is a concrete point TimeSeries, or something that can be evaluated to one,")
                doc_intro("use it as representation")
                doc_intro("for the values of this ts")
                doc_parameters()
                doc_parameter("bts","TimeSeries","a concrete point ts, or ready-to-evaluate expression, with time-axis, values and fx_policy")
                doc_notes()
                doc_note("raises runtime_error if any of preconditions is not true")
                doc_see_also("find_ts_bind_info,TimeSeries('a-ref-string')")
            )
            .def("bind_done",&apoint_ts::do_bind,(py::arg("self")),
                 doc_intro("after bind operations on unbound time-series of an expression is done, call bind_done()")
                 doc_intro("to prepare the expression for use")
                 doc_notes()
                 doc_note("Usually this is done automatically by the dtss framework, but if not using dtss")
                 doc_note("this function is needed *after* the symbolic ts's are bound")
                 doc_see_also(".bind(), .find_ts_bind_info(), needs_bind()")
            )
            .def("needs_bind",&apoint_ts::needs_bind,(py::arg("self")),
                 doc_intro("returns true if there are any unbound time-series in the expression")
                 doc_intro("this time-series represent")
                 doc_see_also(".find_ts_bind_info(),bind() and bind_done()")

            )
            .def("find_ts_bind_info",&apoint_ts::find_ts_bind_info,(py::arg("self")),
                doc_intro("recursive search through the expression that this ts represents,")
                doc_intro("and return a list of TsBindInfo that can be used to")
                doc_intro("inspect and possibly 'bind' to ts-values \ref bind.")
                doc_returns("bind_info","TsBindInfoVector","A list of BindInfo where each entry contains a symbolic-ref and a ts that needs binding")
                doc_see_also("bind() method")

            )
            .def("serialize",&apoint_ts::serialize_to_bytes,(py::arg("self")),
                "convert ts (expression) into a binary blob\n"
            )
            .def("deserialize",&apoint_ts::deserialize_from_bytes,py::args("blob"),
               "convert a blob, as returned by .serialize() into a Timeseries"
            ).staticmethod("deserialize")
        ;
        typedef apoint_ts (*avg_func_t)(const apoint_ts&,const shyft::time_axis::generic_dt&);
        typedef apoint_ts(*int_func_t)(const apoint_ts&, const shyft::time_axis::generic_dt&);
        avg_func_t avg=average;
        int_func_t intfnc = integral;
		avg_func_t acc = accumulate;
        def("average",avg,args("ts","time_axis"),"creates a true average time-series of ts for intervals as specified by time_axis");
        def("integral", intfnc, args("ts", "time_axis"), "creates a true integral time-series of ts for intervals as specified by time_axis");
		def("accumulate", acc, args("ts", "time_axis"), "create a new ts that is the integral f(t) *dt, t0..ti, the specified time-axis");
        //def("max",time_series::dd::max,(boost::python::arg("ts_a"),boost::python::arg("ts_b")),"creates a new time-series that is the max of the supplied ts_a and ts_b");
		def("ts_stringify",ts_stringify , (py::arg("ts")),
			doc_intro("Given a TimeSeries, return a string showing the details/expression")
		);

        typedef apoint_ts(*ts_op_t)(const apoint_ts&a);
        typedef apoint_ts (*ts_op_ts_t)(const apoint_ts&a, const apoint_ts&b);
        typedef apoint_ts (*double_op_ts_t)(double, const apoint_ts&b);
        typedef apoint_ts (*ts_op_double_t)(const apoint_ts&a, double);

        ts_op_ts_t max_ts_ts         = time_series::dd::max;
        double_op_ts_t max_double_ts = time_series::dd::max;
        ts_op_double_t max_ts_double = time_series::dd::max;
        def("max",max_ts_ts    ,args("ts_a","ts_b"),"returns a new ts as max(ts_a,ts_b)");
        def("max",max_double_ts,args("a"   ,"ts_b"),"returns a new ts as max(a,ts_b)");
        def("max",max_ts_double,args("ts_a","b"   ),"returns a new ts as max(ts_a,b)");

        ts_op_ts_t min_ts_ts         = time_series::dd::min;
        double_op_ts_t min_double_ts = time_series::dd::min;
        ts_op_double_t min_ts_double = time_series::dd::min;
        def("min",min_ts_ts    ,args("ts_a","ts_b"),"returns a new ts as min(ts_a,ts_b)");
        def("min",min_double_ts,args("a"   ,"ts_b"),"returns a new ts as min(a,ts_b)");
        def("min",min_ts_double,args("ts_a","b"   ),"returns a new ts as min(ts_a,b)");

        ts_op_ts_t pow_ts_ts         = time_series::dd::pow;
        double_op_ts_t pow_double_ts = time_series::dd::pow;
        ts_op_double_t pow_ts_double = time_series::dd::pow;
        def("pow",pow_ts_ts    ,args("ts_a","ts_b"),"returns a new ts as pow(ts_a,ts_b)");
        def("pow",pow_double_ts,args("a"   ,"ts_b"),"returns a new ts as pow(a,ts_b)");
        def("pow",pow_ts_double,args("ts_a","b"   ),"returns a new ts as pow(ts_a,b)");

        ts_op_t log_ts_ts = time_series::dd::log;
        def("log", log_ts_ts, args("ts"), "returns a new ts as log(ts)");

        def("time_shift", time_series::dd::time_shift,args("timeseries","delta_t"),
            "returns a delta_t time-shifted time-series\n"
            " the values are the same as the original,\n"
            " but the time_axis equals the original + delta_t\n");

        def("create_glacier_melt_ts_m3s", create_glacier_melt_ts_m3s, args("temperature", "sca_m2", "glacier_area_m2", "dtf"),
            doc_intro("create a ts that provide the glacier-melt algorithm based on the inputs")
            doc_parameters()
            doc_parameter("temperature", "TimeSeries", "a temperature time-series, unit [deg.Celcius]")
            doc_parameter("sca_m2", "TimeSeries", "a snow covered area (sca) time-series, unit [m2]")
            doc_parameter("glacier_area_m2", "float", "the glacier area, unit[m2]")
            doc_parameter("dtf","float","degree timestep factor [mm/day/deg.C]; lit. values for Norway: 5.5 - 6.4 in Hock, R. (2003), J. Hydrol., 282, 104-115")
            doc_returns("glacier_melt","TimeSeries","an expression computing the glacier melt based on the inputs")
        );
		/* local scope */ {

			typedef shyft::time_axis::fixed_dt ta_t;
			typedef shyft::time_series::average_accessor<pts_t, ta_t> AverageAccessorTs;
			class_<AverageAccessorTs>("AverageAccessorTs", "Accessor to get out true average for the time-axis intervals for a point time-series", no_init)
				.def(init<const pts_t&, const ta_t&>((py::arg("self"),py::arg("ts"),py::arg("ta")), "construct accessor from ts and time-axis ta"))
				.def(init<shared_ptr<pts_t>, const ta_t&>(( py::arg("self"),py::arg("ts"), py::arg("ta") ), "constructor from ref ts and time-axis ta"))
				.def("value", &AverageAccessorTs::value, (py::arg("self"),py::arg("i")), "returns the i'th true average value")
				.def("size", &AverageAccessorTs::size, (py::arg("self")),"returns number of intervals in the time-axis for this accessor")
				;
		}
		class_<qac_parameter>(
            "QacParameter",
            doc_intro("The qac parameter controls how quailty checks are done, providing min-max range, plus repeated values checks")
            doc_intro("It also provides parameters that controls how the replacement/correction values are filled in,")
            doc_intro("like maximum time-span between two valid neighbour points that allows for linear/extension filling")
            )
            .def(init<utctimespan,double,double,utctimespan,double,double>((py::arg("self"),py::arg("max_timespan"),py::arg("min_x"),py::arg("max_x"),
                                                                            py::arg("repeat_timespan"),py::arg("repeat_tolerance"),py::arg("constant_filler")=shyft::nan),
                                                                           doc_intro("a quite complete qac, only lacks repeat_allowed value(s)"))
            )
            .def(init<utctimespan,double,double,utctimespan,double,double,double>((py::arg("self"),py::arg("max_timespan"),py::arg("min_x"),py::arg("max_x"),
                                                                            py::arg("repeat_timespan"),py::arg("repeat_tolerance"),py::arg("repeat_allowed"),py::arg("constant_filler")=shyft::nan),
                                                                           doc_intro("a quite complete qac, including one repeat_allowed value"))
            )
            .def_readwrite("min_v",&qac_parameter::min_x,"minimum value or nan for no minimum value limit")
            .def_readwrite("max_v",&qac_parameter::max_x,"maximum value or nan for no maximum value limit")
            .def_readwrite("max_timespan",&qac_parameter::max_timespan,"maximum timespan between two ok values that allow interpolation, or extension of values.If zero, no linear/extend correction")
            .def_readwrite("repeat_timespan",&qac_parameter::repeat_timespan,"maximum timespan the same value can be repeated (within repeat_tolerance).If zero, no repeat validation done")
            .def_readwrite("repeat_tolerance",&qac_parameter::repeat_tolerance,"values are considered repeated if they differ by less than repeat_tolerance")
            .def_readwrite("repeat_allowed",&qac_parameter::repeat_allowed,"values that are allowed to repeat, within repeat-tolerance")
            .def_readwrite("constant_filler",&qac_parameter::constant_filler,"this is applied to values that fails quality checks, if no correction ts, and no interpolation/extension is active")
            ;
        
    }

    /** python api need some helper classes to make this more elegant */
    struct rating_curve_t_f {
        utctime t{ no_utctime };
        rating_curve_function f{};
        rating_curve_t_f() = default;
        rating_curve_t_f(utctime t,const rating_curve_function&f) :t{ t },f{ f } {}
        bool operator==(const rating_curve_t_f& o) const { return t==o.t; } // satisfy boost py indexing find
    };

    /** helper to fix custom constructors taking the above f_t as list */
    struct rcp_ext {
        static rating_curve_parameters* create_default() {return new rating_curve_parameters();}
        static rating_curve_parameters* create_from_t_f_list(const vector<rating_curve_t_f>& f_t_list) {
            auto rcp = new rating_curve_parameters();
            for (const auto&e:f_t_list)
                rcp->add_curve(e.t, e.f);
            return rcp;
        }
    };

	static void expose_rating_curve_classes() {

		// overloads for rating_curve_segment::flow
		double (shyft::core::rating_curve_segment::*rcs_flow_1)(double) const = &shyft::core::rating_curve_segment::flow;
		std::vector<double> (shyft::core::rating_curve_segment::*rcs_flow_2)(const std::vector<double> &, std::size_t, std::size_t) const = &shyft::core::rating_curve_segment::flow;

		class_<shyft::core::rating_curve_segment>("RatingCurveSegment",
				doc_intro("Represent a single rating-curve equation.")
				doc_intro("")
				doc_intro("The rating curve function is `a*(h - b)^c` where `a`, `b`, and `c` are parameters")
				doc_intro("for the segment and `h` is the water level to compute flow for. Additionally there")
				doc_intro("is a `lower` parameter for the least water level the segment is valid for. Seen")
				doc_intro("separatly a segment is considered valid for any level greater than `lower`.")
				doc_intro("")
				doc_intro("The function segments are gathered into `RatingCurveFunction`s to represent a")
				doc_intro("set of different rating functions for different levels.")
				doc_see_also("RatingCurveFunction, RatingCurveParameters"),
				init<>(py::arg("self"))
			)
			.def_readonly("lower", &shyft::core::rating_curve_segment::lower,
						  "Least valid water level. Not mutable after constructing a segment.")
			.def_readwrite("a", &shyft::core::rating_curve_segment::a, "Parameter a")
			.def_readwrite("b", &shyft::core::rating_curve_segment::b, "Parameter b")
			.def_readwrite("c", &shyft::core::rating_curve_segment::c, "Parameter c")
			.def(init<double, double, double, double>( (py::arg("self"),py::arg("lower"),py::arg("a"),py::arg("b"),py::arg("c") ), "Defines a new RatingCurveSegment with the specified parameters"))
			.def("valid", &shyft::core::rating_curve_segment::valid, (py::arg("self"),py::arg("level")),
					doc_intro("Check if a water level is valid for the curve segment")
					doc_parameter("level", "float", "water level")
					doc_returns("valid", "bool", "True if level is greater or equal to lower")
				)
            //NOTE: For some reason boost 1.65 needs this def *before* the other simpler def, otherwise it fails finding the simple one
            .def("flow", rcs_flow_2, (py::arg("self"),py::arg("levels"), py::arg("i0") = 0u, py::arg("iN") = std::numeric_limits<std::size_t>::max()),
                doc_intro("Compute the flow for a range of water levels")
                doc_parameters()
                doc_parameter("levels", "DoubleVector", "Vector of water levels")
                doc_parameter("i0", "int", "first index to use from levels, defaults to 0")
                doc_parameter("iN", "int", "first index _not_ to use from levels, defaults to std::size_t maximum.")
                doc_returns("flow", "DoubleVector", "Vector of flow values.")
            )
            .def("flow", rcs_flow_1, (py::arg("self"),py::arg("level")),
					doc_intro("Compute the flow for the given water level.")
					doc_notes()
					doc_note("There is _no_ check to see if level is valid. It's up to the user to call")
					doc_note("with a correct level.")
					doc_parameters()
					doc_parameter("level", "float", "water level")
					doc_returns("flow", "double", "the flow for the given water level")
				)
			.def("__str__", &shyft::core::rating_curve_segment::operator std::string, "Stringify the segment.")
            .def(self==self)
            .def(self!=self)
			;

        typedef std::vector<shyft::core::rating_curve_segment> RatingCurveSegmentVector;
        class_<RatingCurveSegmentVector>("RatingCurveSegments",
              doc_intro("A typed list of RatingCurveSegment, used to construct RatingCurveParameters."),
              init<>(py::arg("self"))
            )
            .def(vector_indexing_suite<RatingCurveSegmentVector>())
            .def(init<const RatingCurveSegmentVector&>(args("clone_me")))
        ;
        py_api::iterable_converter().from_python<RatingCurveSegmentVector>();
		// overloads for rating_curve_function::flow
		double (shyft::core::rating_curve_function::*rcf_flow_val)(double) const = &shyft::core::rating_curve_function::flow;
		std::vector<double> (shyft::core::rating_curve_function::*rcf_flow_vec)(const std::vector<double> & ) const = &shyft::core::rating_curve_function::flow;
		// overloads for rating_curve_function::add_segment
		void (shyft::core::rating_curve_function::*rcf_add_args)(double, double, double, double) = &shyft::core::rating_curve_function::add_segment;
		void (shyft::core::rating_curve_function::*rcf_add_obj)(const rating_curve_segment & ) = &shyft::core::rating_curve_function::add_segment;

		class_<shyft::core::rating_curve_function>("RatingCurveFunction",
				doc_intro("Combine multiple RatingCurveSegments into a rating function.")
				doc_intro("")
				doc_intro("RatingCurveFunction aggregates multiple RatingCurveSegments and routes.")
				doc_intro("computation calls to the correct segment based on the water level to compute for.")
				doc_see_also("RatingCurveSegment, RatingCurveParameters"),
				init<>(py::arg("self"),doc_intro("Defines a new empty rating curve function."))
			)
            .def(init<const RatingCurveSegmentVector&,bool>(
                (py::arg("self"),py::arg("segments"),py::arg("is_sorted")=true),
                doc_intro("constructs a function from a segment-list")
                )
            )
			.def("size", &shyft::core::rating_curve_function::size,(py::arg("self")), "Get the number of RatingCurveSegments composing the function.")
			.def("add_segment", rcf_add_args, (py::arg("self"),py::arg("lower"), py::arg("a"), py::arg("b"), py::arg("c")),
					doc_intro("Add a new curve segment with the given parameters.")
					doc_see_also("RatingCurveSegment")
				)
			.def("add_segment", rcf_add_obj, (py::arg("self"),py::arg("segment")),
					doc_intro("Add a new curve segment as a copy of an exting.")
					doc_see_also("RatingCurveSegment")
				)
            // ref. note above regarding the order of overloaded member functions
            .def("flow", rcf_flow_vec, (py::arg("self"), py::arg("levels")),
                doc_intro("Compute flow for a range of water levels.")
                doc_parameters()
                doc_parameter("levels", "DoubleVector", "Range of water levels to compute flow for.")
            )
            .def("flow", rcf_flow_val,(py::arg("self"), py::arg("level")),
					doc_intro("Compute flow for the given level.")
					doc_parameters()
					doc_parameter("level", "float", "Water level to compute flow for.")
				)
			.def("__iter__", py::range(&shyft::core::rating_curve_function::cbegin,
									   &shyft::core::rating_curve_function::cend),
				 "Constant iterator. Invalidated on calls to .add_segment")
			.def("__str__", &shyft::core::rating_curve_function::operator std::string, "Stringify the function.")
			;

        class_<rating_curve_t_f>("RatingCurveTimeFunction",
                doc_intro("Composed of time t and RatingCurveFunction"),
                init<>(py::arg("self"), doc_intro("Defines empty pair t,f"))
            )
            .def(init<utctime,const rating_curve_function&>((py::arg("self"), py::arg("t"), py::arg("f")),
                doc_intro("Construct an object with function f valid from time t")
                doc_parameters()
                doc_parameter("t", "int", "epoch time in 1970 utc [s]")
                doc_parameter("f", "RatingCurveFunction", "the function")
                )
            )
            .def_readwrite("t",&rating_curve_t_f::t,doc_intro(".f is valid from t, the epoch 1970[s] time "))
            .def_readwrite("f", &rating_curve_t_f::f, doc_intro("the rating curve function"))
            ;

        typedef vector<rating_curve_t_f> RatingCurveTimeFunctionVector;
        class_<RatingCurveTimeFunctionVector>("RatingCurveTimeFunctions",
            doc_intro("A typed list of RatingCurveTimeFunction elements"),
            init<>(py::arg("self"), doc_intro("Defines empty list pair t,f"))
            )
            .def(vector_indexing_suite<RatingCurveTimeFunctionVector>())
            .def(init<const RatingCurveTimeFunctionVector&>(args("clone_me")))
            ;
        py_api::iterable_converter().from_python<RatingCurveTimeFunctionVector>();

		// overloads for rating_curve_function::flow
		double (shyft::core::rating_curve_parameters::*rcp_flow_val)(utctime, double) const = &shyft::core::rating_curve_parameters::flow;
		std::vector<double> (shyft::core::rating_curve_parameters::*rcp_flow_ts)(const apoint_ts & ) const = &shyft::core::rating_curve_parameters::flow<apoint_ts>;
		// overloads for rating_curve_function::add_segment
        rating_curve_parameters& (shyft::core::rating_curve_parameters::*rcp_add_obj)(utctime, const rating_curve_function & ) = &shyft::core::rating_curve_parameters::add_curve;

		class_<shyft::core::rating_curve_parameters>("RatingCurveParameters",
				doc_intro("Parameter pack controlling rating level computations.")
				doc_intro("")
				doc_intro("A parameter pack encapsulates multiple RatingCurveFunction's with time-points.")
				doc_intro("When used with a TimeSeries representing level values it maps computations for")
				doc_intro("each level value onto the correct RatingCurveFunction, which again maps onto the")
				doc_intro("correct RatingCurveSegment for the level value.")
				doc_see_also("RatingCurveSegment, RatingCurveFunction, TimeSeries.rating_curve"),
				no_init //init<>((py::arg("self")),"Defines a empty RatingCurveParameter instance")
			)
            .def("__init__",make_constructor(&rcp_ext::create_default),"Defines a empty RatingCurveParameter instance")
            .def("__init__",make_constructor(
                    &rcp_ext::create_from_t_f_list,
                    default_call_policies(),
                    (py::arg("t_f_list"))
                ),
                doc_intro("create parameters in one go from list of RatingCurveTimeFunction elements")
                doc_parameters()
                doc_parameter("t_f_list","RatingCurveTimeFunctions","a list of RatingCurveTimeFunction elements")
            )
			.def("add_curve", rcp_add_obj, (py::arg("self"),py::arg("t"), py::arg("curve")),
					doc_intro("Add a curve to the parameter pack.")
					doc_parameters()
					doc_parameter("t", "RatingCurveFunction", "First time-point the curve is valid for.")
					doc_parameter("curve", "RatingCurveFunction", "RatingCurveFunction to add at t.")
                    doc_returns("self","RatingCurveParameters"," to allow chaining building functions"),
                 return_value_policy<reference_existing_object>()

				)
			.def("flow", rcp_flow_val, (py::arg("self"), py::arg("t"), py::arg("level")),
					doc_intro("Compute the flow at a specific time point.")
					doc_parameters()
					doc_parameter("t", "utctime", "Time-point of the level value.")
					doc_parameter("level", "float", "Level value at t.")
					doc_returns("flow", "float", "Flow correcponding to input level at t, `nan` if level is less than the least water level of the first segment or before the time of the first rating curve function.")
				)
			.def("flow", rcp_flow_ts,(py::arg("self"), py::arg("ts")),
					doc_intro("Compute the flow at a specific time point.")
					doc_parameters()
					doc_parameter("ts", "TimeSeries", "Time series of level values.")
					doc_returns("flow", "DoubleVector", "Flow correcponding to the input levels of the time-series, `nan` where the level is less than the least water level of the first segment and for time-points before the first rating curve function.")
				)
			.def("__iter__", py::range(&shyft::core::rating_curve_parameters::cbegin,
									   &shyft::core::rating_curve_parameters::cend),
				 "Constant iterator. Invalidated on calls to .add_curve")
			.def("__str__", &shyft::core::rating_curve_parameters::operator std::string, "Stringify the parameters.")
			;
	}

    static void expose_ice_packing_parameters() {
        enum_<shyft::core::ice_packing_temperature_policy>("ice_packing_temperature_policy",
            doc_intro("Policy enum to specify how `TimeSeries.ice_packing` handles missing temperature values.")
            doc_intro("")
            doc_intro("The enum defines three values:")
            doc_intro(" * `DISALLOW_MISSING` disallows any missing values. With this policy whenever a NaN value is encountered,")
            doc_intro("   or the window of values to consider extends outside the range of the time series, a NaN value will be")
            doc_intro("   written to the result time-series.")
            doc_intro(" * `ALLOW_INITIAL_MISSING` disallows explicit NaN values, but allows the window of values to consider")
            doc_intro("   to expend past the range of the time-series for the initial values.")
            doc_intro(" * `ALLOW_ANY_MISSING` allow the window of values to contain NaN values, averaging what it can.")
            doc_intro("   Only if all the values in the window is NaN, the result wil be NaN.")
            )
            .value("DISALLOW_MISSING", shyft::core::ice_packing_temperature_policy::DISALLOW_MISSING)
            .value("ALLOW_INITIAL_MISSING", shyft::core::ice_packing_temperature_policy::ALLOW_INITIAL_MISSING)
            .value("ALLOW_ANY_MISSING", shyft::core::ice_packing_temperature_policy::ALLOW_ANY_MISSING)
            .export_values()
            ;

        class_<shyft::core::ice_packing_parameters>("IcePackingParameters",
            doc_intro("Parameter pack controlling ice packing computations.")
            doc_intro("See `TimeSeries.ice_packing` for usage."),
            init<core::utctimespan, double>((py::arg("self"), py::arg("threshold_window"), py::arg("threshold_temperature")),
                doc_intro("Defines a paramter pack for ice packing detection.")
                doc_intro("")
                doc_parameters()
                doc_parameter("threshold_window", "utctime", "Positive,  seconds for the lookback window.")
                doc_parameter("threshold_temperature", "float", "Floating point threshold temperature."))
            )
            .def(init<int64_t, double>((py::arg("self"), py::arg("threshold_window"), py::arg("threshold_temperature")),
                doc_intro("Defines a paramter pack for ice packing detection.")
                doc_intro("")
                doc_parameters()
                doc_parameter("threshold_window", "int", "Positive integer seconds for the lookback window.")
                doc_parameter("threshold_temperature", "float", "Floating point threshold temperature."))
            )
            .def_readwrite("threshold_window",
                &shyft::core::ice_packing_parameters::window,
                doc_intro("The period back in seconds for which the average temperature is computed when")
                doc_intro("looking for ice packing.")
            )
            .def_readwrite("threshold_temperature",
                &shyft::core::ice_packing_parameters::threshold_temp,
                doc_intro("The threshold temperature for ice packing to occur. Ice packing will occur")
                doc_intro("when the average temperature in the `window` period is less than the threshold.")
            )
            .def(self==self)
            ;

        class_<shyft::time_series::dd::ice_packing_recession_parameters>("IcePackingRecessionParameters",
            doc_intro("Parameter pack controlling ice packing recession computations.")
            doc_intro("See `TimeSeries.ice_packing_recession` for usage."),
            init<double, double>((py::arg("self"), py::arg("alpha"), py::arg("recession_minimum")),
                doc_intro("Defines a parameter pack for ice packing reduction using a simple recession for the water-flow.")
                doc_intro("")
                doc_parameters()
                doc_parameter("alpha", "float", "Recession curve curving parameter.")
                doc_parameter("recession_minimum", "float", "Minimum value for the recession."))
            )
            .def_readwrite("alpha",
                &shyft::time_series::dd::ice_packing_recession_parameters::alpha,
                doc_intro("Parameter controlling the curving of the recession curve.")
            )
            .def_readwrite("recession_minimum",
                &shyft::time_series::dd::ice_packing_recession_parameters::recession_minimum,
                doc_intro("The minimum value of the recession curve.")
            )
            .def(self==self)
            ;
    }

	static void expose_correlation_functions() {
		const char * kg_doc =
			doc_intro("Computes the kling-gupta KGEs correlation for the two time-series over the specified time_axis")
			doc_parameters()
			doc_parameter("observed_ts","TimeSeries","the observed time-series")
			doc_parameter("model_ts","TimeSeries","the time-series that is the model simulated / calculated ts")
			doc_parameter("time_axis","TimeAxis","the time-axis that is used for the computation")
			doc_parameter("s_r","float","the kling gupta scale r factor(weight the correlation of goal function)")
			doc_parameter("s_a","float","the kling gupta scale a factor(weight the relative average of the goal function)")
			doc_parameter("s_b","float","the kling gupta scale b factor(weight the relative standard deviation of the goal function)")
            doc_returns("KGEs","float","The  KGEs= 1-EDs that have a maximum at 1.0");
		def("kling_gupta", kling_gupta, args("observation_ts", "model_ts", "time_axis", "s_r", "s_a", "s_b"),
			kg_doc
		);

		const char *ns_doc =
            doc_intro("Computes the Nash-Sutcliffe model effiency coefficient (n.s) ")
			doc_intro("for the two time-series over the specified time_axis\n")
			doc_intro("Ref:  http://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient \n")
            doc_parameters()
			doc_parameter("observed_ts","TimeSeries","the observed time-series")
			doc_parameter("model_ts","TimeSeries","the time-series that is the model simulated / calculated ts")
			doc_parameter("time_axis","TimeAxis","the time-axis that is used for the computation")
			doc_returns("ns","float","The  n.s performance, that have a maximum at 1.0");

		def("nash_sutcliffe", nash_sutcliffe, args("observation_ts", "model_ts", "time_axis"),
			ns_doc
		);


	}
	static void expose_periodic_ts() {
		const char *docstr =
			doc_intro("Create a Timeseries by repeating the pattern-specification")
			doc_parameters()
			doc_parameter("pattern","DoubleVector","the value-pattern as a sequence of values")
            doc_parameter("dt","int","number of seconds between the pattern values, e.g. deltahours(3)")
            doc_parameter("t0","utctime","specifies the start-time of the pattern")
            doc_parameter("ta","TimeAxis","the time-axis for which the pattern is repeated\n\te.g. your pattern might be 8 3h values,and you could supply\n\ta time-axis 'ta' at hourly resolution")
			;
		def("create_periodic_pattern_ts", create_periodic_pattern_ts, args("pattern","dt","t0","ta"), docstr);

	}

	static void expose_krls() {

        using krls_rbf_predictor = shyft::prediction::krls_rbf_predictor;

		class_<krls_rbf_predictor>("KrlsRbfPredictor",
			    doc_intro("Time-series predictor using the KRLS algorithm with radial basis functions.")
                doc_intro("")
                doc_intro("The KRLS (Kernel Recursive Least-Squares) algorithm is a kernel regression")
                doc_intro("algorithm for aproximating data, the implementation used here is from")
                doc_intro("DLib: http://dlib.net/ml.html#krls.")
                doc_intro("This predictor uses KRLS with radial basis functions (RBF).")
                doc_intro("")
                doc_intro("Examples\n--------\n")
                doc_intro(">>>")
                doc_intro(">>> import numpy as np")
                doc_intro(">>> import matplotlib.pyplot as plt")
                doc_intro(">>> from shyft.api import (")
                doc_intro("...     Calendar, utctime_now, deltahours,")
                doc_intro("...     TimeAxis, TimeSeries,")
                doc_intro("...     KrlsRbfPredictor")
                doc_intro("... )")
                doc_intro(">>>")
                doc_intro(">>> # setup")
                doc_intro(">>> cal = Calendar()")
                doc_intro(">>> t0 = utctime_now()")
                doc_intro(">>> dt = deltahours(3)")
                doc_intro(">>> n = 365*8  # one year")
                doc_intro(">>>")
                doc_intro(">>> # ready plot")
                doc_intro(">>> fig, ax = plt.subplots()")
                doc_intro(">>> ")
                doc_intro(">>> # shyft objects")
                doc_intro(">>> ta = TimeAxis(t0, dt, n)")
                doc_intro(">>> pred = KrlsRbfPredictor(")
                doc_intro("...     dt=deltahours(8),")
                doc_intro("...     gamma=1e-5,  # NOTE: this should be 1e-3 for real data")
                doc_intro("...     tolerance=0.001")
                doc_intro("... )")
                doc_intro(">>>")
                doc_intro(">>> # generate data")
                doc_intro(">>> total_series = 4")
                doc_intro(">>> data_range = np.linspace(0, 2*np.pi, n)")
                doc_intro(">>> ts = None  # to store the final data-ts")
                doc_intro(">>> # -----")
                doc_intro(">>> for i in range(total_series):")
                doc_intro(">>>     data = np.sin(data_range) + (np.random.random(data_range.shape) - 0.5)/5")
                doc_intro(">>>     ts = TimeSeries(ta, data)")
                doc_intro(">>>     # -----")
                doc_intro(">>>     training_mse = pred.train(ts)  # train the predictor")
                doc_intro(">>>     # -----")
                doc_intro(">>>     print(f'training step {i+1:2d}: mse={training_mse}')")
                doc_intro(">>>     ax.plot(ta.time_points[:-1], ts.values, 'bx')  # plot data")
                doc_intro(">>>")
                doc_intro(">>> # prediction")
                doc_intro(">>> ts_pred = pred.predict(ta)")
                doc_intro(">>> ts_mse = pred.mse_ts(ts, points=3)  # mse using 7 point wide filter")
                doc_intro(">>>                                     # (3 points before/after)")
                doc_intro(">>>")
                doc_intro(">>> # plot interpolation/predicton on top of results")
                doc_intro(">>> ax.plot(ta.time_points[:-1], ts_mse.values, '0.6', label='mse')")
                doc_intro(">>> ax.plot(ta.time_points[:-1], ts_pred.values, 'r-', label='prediction')")
                doc_intro(">>> ax.legend()")
                doc_intro(">>> plt.show()")
                doc_see_also("TimeSeries.krls_interpolation, TimeSeries.get_krls_predictor")
			)
			.def(init<core::utctimespan, double, double, std::size_t>(
                    ( py::arg("self"), py::arg("dt"), py::arg("gamma") = 1.E-3, py::arg("tolerance") = 0.01, py::arg("size") = 1000000u ),
                    doc_intro("Construct a new predictor.")
                    doc_parameters()
                    doc_parameter("dt", "float", "The time-step in seconds the predictor is specified for.\n"
                        "    Note that this does not put a limit on time-axes used, but for best results it should be\n"
                        "    approximatly equal to the time-step of time-axes used with the predictor. In addition it\n"
                        "    should not be to long, else you will get poor results. Try to keep the dt less than a day,\n"
                        "    3-8 hours is usually fine." )
                    doc_parameter("gamma", "float (optional)", "Determines the width of the radial basis functions for\n"
                        "    the KRLS algorithm. Lower values mean wider basis functions, wider basis functions means faster\n"
                        "    computation but lower accuracy. Note that the tolerance parameter also affects speed and accurcy.\n"
                        "    A large value is around `1E-2`, and a small value depends on the time step. By using values larger\n"
                        "    than `1E-2` the computation will probably take to long. Testing have reveled that `1E-3` works great\n"
                        "    for a time-step of 3 hours, while a gamma of `1E-2` takes a few minutes to compute. Use `1E-4` for a\n"
                        "    fast and tolerably accurate prediction.\n"
                        "    Defaults to `1E-3`" )
                    doc_parameter("tolerance", "float (optional)", "The krls training tolerance. Lower values makes the prediction more accurate,\n"
                        "    but slower. This typically have less effect than gamma, but is usefull for tuning. Usually it should be\n"
                        "    either `0.01` or `0.001`.\n"
                        "    Defaults to `0.01`" )
                    doc_parameter("size", "int (optional)", "The size of the \"memory\" of the predictor. The default value is\n"
                        "    usually enough. Defaults to `1000000`." )
                ) )
			.def("train", &krls_rbf_predictor::train<apoint_ts>,
				    ( py::arg("self"),
					  py::arg("ts"),
                      py::arg("offset") = 0u, py::arg("count") = std::numeric_limits<std::size_t>::max(), py::arg("stride") = 1u,
                      py::arg("iterations") = 1u, py::arg("mse_tol") = 0.001 ),
				    doc_intro("Train the predictor using samples from ts.")
                    doc_parameters()
                    doc_parameter("ts", "TimeSeries", "Time-series to train on.")
                    doc_parameter("offset", "int (optional)", "Positive offset from the start of the time-series. Default to 0.")
                    doc_parameter("count", "int (optional)", "Positive number of samples to to use. Default to the maximum value.")
                    doc_parameter("stride", "int (optional)", "Positive stride between samples from the time-series. Defaults to 1.")
                    doc_parameter("iterations", "int (optional)", "Positive maximum number of times to train on the samples. Defaults to 1.")
                    doc_parameter("mse_tol", "float (optional)", "Positive tolerance for the mean-squared error over the training data.\n"
                        "    If the mse after a training session is less than this skip training further. Defaults to `1E-9`." )
                    doc_returns("mse", "float (optional)", "Mean squared error of the predictor relative to the time-series trained on.")
                )
			.def("predict", &krls_rbf_predictor::predict<apoint_ts,shyft::time_series::dd::gta_t>,
                    (py::arg("self"), py::arg("ta")),
				    doc_intro("Predict a time-series for for time-axis.")
                    doc_notes()
                    doc_note("The predictor will predict values outside the range of the values it is trained on, but these")
                    doc_note("values will often be zero. This may also happen if there are long gaps in the training data")
                    doc_note("and you try to predict values for the gap. Using wider basis functions partly remedies this,")
                    doc_note("but makes the prediction overall less accurate.")
                    doc_parameters()
                    doc_parameter("ta", "TimeAxis", "Time-axis to predict values for.")
                    doc_returns("ts", "TimeSeries", "Predicted time-series.")
                    doc_see_also("KrlsRbfPredictor.mse_ts, KrlsRbfPredictor.predictor_mse")
			    )
            .def("mse_ts", &krls_rbf_predictor::mse_ts<apoint_ts, apoint_ts>,
                    ( py::arg("self"), py::arg("ts"), py::arg("points") = 0u ),
                    doc_intro("Compute a mean-squared error time-series of the predictor relative to the supplied ts.")
                    doc_parameters()
                    doc_parameter("ts", "TimeSeries", "Time-series to compute mse against.")
                    doc_parameter("points", "int (optional)", "Positive number of extra points around each point to use for mse.\n"
                        "    Defaults to 0." )
                    doc_returns("mse_ts", "TimeSeries", "Time-series with mean-squared error values.")
                    doc_see_also("KrlsRbfPredictor.predictor_mse, KrlsRbfPredictor.predict")
                )
			.def("predictor_mse", &krls_rbf_predictor::predictor_mse<apoint_ts>,
                    ( py::arg("self"),
					  py::arg("ts"),
                      py::arg("offset") = 0u, py::arg("count") = std::numeric_limits<std::size_t>::max(), py::arg("stride") = 1u ),
                    doc_intro("Compute the predictor mean-squared prediction error for count first from ts.")
                    doc_parameters()
                    doc_parameter("ts", "TimeSeries", "Time-series to compute mse against.")
                    doc_parameter("offset", "int (optional)", "Positive offset from the start of the time-series. Default to 0.")
                    doc_parameter("count", "int (optional)", "Positive number of samples from the time-series to to use.\n"
                        "    Default to the maximum value." )
                    doc_parameter("stride", "int (optional)", "Positive stride between samples from the time-series. Defaults to 1.")
                    doc_see_also("KrlsRbfPredictor.predict, KrlsRbfPredictor.mse_ts")
            )
            .def("clear", &krls_rbf_predictor::clear,(py::arg("self")),
                    doc_intro("Clear all training data from the predictor.")
                )
			;
	}

    void timeseries() {
        enum_<time_series::ts_point_fx>("point_interpretation_policy")
            .value("POINT_INSTANT_VALUE",time_series::POINT_INSTANT_VALUE)
            .value("POINT_AVERAGE_VALUE",time_series::POINT_AVERAGE_VALUE)
            .export_values()
            ;
        enum_<time_series::statistics_property>("statistics_property")
            .value("AVERAGE",time_series::statistics_property::AVERAGE)
            .value("MIN_EXTREME",time_series::statistics_property::MIN_EXTREME)
            .value("MAX_EXTREME",time_series::statistics_property::MAX_EXTREME)
            ;

        enum_<time_series::dd::extend_ts_fill_policy>(
            "extend_fill_policy",
            "Ref TimeSeries.extend function, this policy determines how to represent values in a gap\n"
            "EPF_NAN : use nan values in the gap\n"
            "EPF_LAST: use the last value before the gap\n"
            "EPF_FILL: use a supplied value in the gap\n"
            )
            .value("FILL_NAN",   time_series::dd::extend_ts_fill_policy::EPF_NAN)
            .value("USE_LAST",   time_series::dd::extend_ts_fill_policy::EPF_LAST)
            .value("FILL_VALUE", time_series::dd::extend_ts_fill_policy::EPF_FILL)
            .export_values()
            ;
        enum_<time_series::dd::extend_ts_split_policy>(
            "extend_split_policy",
            "Ref TimeSeries.extend function, this policy determines where to split/shift from one ts to the other\n"
            "EPS_LHS_LAST : use nan values in the gap\n"
            "EPS_RHS_FIRST: use the last value before the gap\n"
            "EPS_VALUE    : use a supplied value in the gap\n"
            )
            .value("LHS_LAST",  time_series::dd::extend_ts_split_policy::EPS_LHS_LAST)
            .value("RHS_FIRST", time_series::dd::extend_ts_split_policy::EPS_RHS_FIRST)
            .value("AT_VALUE",  time_series::dd::extend_ts_split_policy::EPS_VALUE)
            .export_values()
            ;

        enum_<time_series::convolve_policy>(
            "convolve_policy",
            "Ref Timeseries.convolve_w function, this policy determinte how to handle initial conditions\n"
            "USE_NEAREST: value(0) is used for all values before value(0), and\n"
            "             value(n-1) is used for all values after value(n-1) == 'mass preserving'\n"
            "USE_ZERO : use zero for all values before value(0) or after value(n-1) == 'shape preserving'\n"
            "USE_NAN  : nan is used for all values outside the ts\n"
            "BACKWARD : filter is 'backward looking' == boundary handling in the beginning of ts\n"
            "FORWARD  : filter is 'forward looking' == boundary handling in the end of ts\n"
            "CENTER   : filter is centered == boundary handling in both ends\n"
            )
            .value("USE_NEAREST", time_series::convolve_policy::USE_NEAREST)
            .value("USE_ZERO", time_series::convolve_policy::USE_ZERO)
            .value("USE_NAN", time_series::convolve_policy::USE_NAN)
            .value("BACKWARD", time_series::convolve_policy::BACKWARD)
            .value("FORWARD", time_series::convolve_policy::FORWARD)
            .value("CENTER", time_series::convolve_policy::CENTER)
            .export_values()
            ;
        enum_<time_series::dd::derivative_method>(
            "derivative_method",
            doc_intro("Ref. the .derivative time-series function, this defines how to compute the")
            doc_intro("derivative of a given time-series")
            )
            .value("DEFAULT",time_series::dd::derivative_method::default_diff)
            .value("FORWARD",time_series::dd::derivative_method::forward_diff)
            .value("BACKWARD",time_series::dd::derivative_method::backward_diff)
            .value("CENTER",time_series::dd::derivative_method::center_diff)
            .export_values()
        ;
        class_<time_series::point> ("Point", "A timeseries point specifying utctime t and value v")
            .def(init<utctime,double>(args("t","v")))
            .def_readwrite("t",&time_series::point::t)
            .def_readwrite("v",&time_series::point::v)
            ;
        point_ts<time_axis::fixed_dt>("TsFixed","A time-series with a fixed delta t time-axis, used by the Shyft core,see also TimeSeries for end-user ts");
        point_ts<time_axis::point_dt>("TsPoint","A time-series with a variable delta time-axis, used by the Shyft core,see also TimeSeries for end-user ts");
        TsFactory();
		expose_rating_curve_classes();
        expose_ice_packing_parameters();
        expose_apoint_ts();
		expose_periodic_ts();
		expose_correlation_functions();
		expose_core_ts_vector();
		expose_ats_vector();
		expose_krls();
    }
}
