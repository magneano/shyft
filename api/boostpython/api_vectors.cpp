/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#include "boostpython_pch.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include "numpy_boost_python.hpp"

#include "py_convertible.h"
#include "core/utctime_utilities.h"
#include "core/geo_point.h"
#include "core/geo_cell_data.h"
#include "core/time_series_dd.h"
#include "api/api.h"

//#include <boost/python/raw_function.hpp>

namespace expose {
    using namespace shyft::core;
    using namespace boost::python;
    using namespace std;
    namespace sa=shyft::api;
    namespace sc=shyft::core;
    namespace ts=shyft::time_series;
    namespace py=boost::python;
	using shyft::time_series::dd::ats_vector;
	using shyft::time_series::dd::gta_t;

    static void* np_import() {
        import_array();
        return nullptr;
    }

    ats_vector create_tsv_from_np(const gta_t& ta, const numpy_boost<double,2>& a ,ts::ts_point_fx point_fx) {
        ats_vector r;
        size_t n_ts = a.shape()[0];
        size_t n_pts = a.shape()[1];
        if(ta.size() != n_pts)
            throw std::runtime_error("time-axis should have same length as second dim in numpy array");
        r.reserve(n_ts);
        for(size_t i=0;i<n_ts;++i) {
            std::vector<double> v;v.reserve(n_pts);
            for(size_t j=0;j<n_pts;++j) v.emplace_back(a[i][j]);
            r.emplace_back(ta,v ,point_fx);
        }
        return r;
    }

    template <class S>
    vector<S> create_from_geo_tsv_from_np(const gta_t& ta,const vector<sc::geo_point>&gpv ,const numpy_boost<double,2>& a ,ts::ts_point_fx point_fx) {
        vector<S> r;
        size_t n_ts = a.shape()[0];
        size_t n_pts = a.shape()[1];
        if(ta.size() != n_pts)
            throw std::runtime_error("time-axis should have same length as second dim in numpy array");
        if(n_ts != gpv.size())
            throw std::runtime_error("geo-point vector should have same size as first dim (n_ts) in numpy array");
        r.reserve(n_ts);
        for(size_t i=0;i<n_ts;++i) {
            std::vector<double> v;v.reserve(n_pts);
            for(size_t j=0;j<n_pts;++j) v.emplace_back(a[i][j]);
            r.emplace_back(gpv[i], sa::apoint_ts(ta,v ,point_fx));
        }
        return r;
    }

    template<class T>
    static vector<T> FromNdArray(const numpy_boost<T,1>& npv) {
        vector<T> r;r.reserve(npv.shape()[0]);
        for(size_t i=0;i<npv.shape()[0];++i) {
            r.push_back(npv[i]);
        }
        return r;
    }

    template<class T>
    static numpy_boost<T,1> ToNpArray(const vector<T>&v) {
        int dims[]={int(v.size())};
        numpy_boost<T,1> r(dims);
        for(size_t i=0;i<r.size();++i) {
            r[i]=v[i];
        }
        return r;
    }

    template <class T>
    static void expose_vector(const char *name) {
        typedef std::vector<T> XVector;

        class_<XVector>(name)
        .def(vector_indexing_suite<XVector>()) // meaning it get all it needs to appear as python list
        .def(init<const XVector&>(args("const_ref_v"))) // so we can copy construct
        .def("FromNdArray",FromNdArray<T>).staticmethod("FromNdArray") // BW compatible
        .def("from_numpy",FromNdArray<T>).staticmethod("from_numpy")// static construct from numpy TODO: fix __init__
        .def("to_numpy",ToNpArray<T>,"convert to numpy") // Ok, to make numpy 1-d arrays
        ;
        numpy_boost_python_register_type<T, 1>(); // register the numpy object so we can access it in C++
        py_api::iterable_converter().from_python<XVector>();
    }
    static void expose_str_vector(const char *name) {
        typedef std::vector<std::string> XVector;
        class_<XVector>(name)
            .def(vector_indexing_suite<XVector>()) // meaning it get all it needs to appear as python list
            .def(init<const XVector&>(args("const_ref_v"))) // so we can copy construct
            ;
        py_api::iterable_converter().from_python<XVector>();
    }

    typedef std::vector<shyft::core::geo_point> GeoPointVector;
    static GeoPointVector create_from_x_y_z_vectors(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z) {
        if(!(x.size()==y.size() && y.size()==z.size()))
            throw std::runtime_error("x,y,z vectors need to have same number of elements");
        GeoPointVector r;r.reserve(x.size());
        for(size_t i=0;i<x.size();++i)
            r.emplace_back(x[i],y[i],z[i]);
        return r;
    }

    static void expose_geo_point_vector() {

        class_<GeoPointVector>("GeoPointVector", "A vector, list, of GeoPoints")
            .def(vector_indexing_suite<GeoPointVector>())
            .def(init<const GeoPointVector&>(args("const_ref_v")))
            .def("create_from_x_y_z",create_from_x_y_z_vectors,args("x","y","z"),"Create a GeoPointVector from x,y and z DoubleVectors of equal length")
            .staticmethod("create_from_x_y_z")
            ;
            py_api::iterable_converter().from_python<GeoPointVector>();
    }
    static void expose_geo_cell_data_vector() {
        typedef std::vector<shyft::core::geo_cell_data> GeoCellDataVector;
        class_<GeoCellDataVector>("GeoCellDataVector", "A vector, list, of GeoCellData")
            .def(vector_indexing_suite<GeoCellDataVector>())
            .def(init<const GeoCellDataVector&>(args("const_ref_v")))
            ;
        py_api::iterable_converter().from_python<GeoCellDataVector>();
    }
    static void expose_ts_vector_create() {
        def("create_ts_vector_from_np_array",
            &create_tsv_from_np,(py::arg("time_axis"),py::arg("np_array"),py::arg("point_fx")),
            doc_intro("Create a TsVector from specified time_axis, 2-d np_array and point_fx.")
            doc_parameters()
            doc_parameter("time_axis","TimeAxis","time-axis that matches in length to 2nd dim of np_array")
            doc_parameter("np_array","np.ndarray","numpy array of dtype=np.float64, and shape(n_ts,n_points)")
            doc_parameter("point_fx","point interpretation", "one of POINT_AVERAGE_VALUE|POINT_INSTANT_VALUE")
            doc_returns("tsv","TsVector","a TsVector of length first np_array dim, n_ts, each with time-axis, values and point_fx")
        );

        def("create_temperature_source_vector_from_np_array",
            &create_from_geo_tsv_from_np<sa::TemperatureSource>,(py::arg("time_axis"),py::arg("geo_points"),py::arg("np_array"),py::arg("point_fx")),
            doc_intro("Create a TemperatureSourceVector from specified time_axis,geo_points, 2-d np_array and point_fx.")
            doc_parameters()
            doc_parameter("time_axis","TimeAxis","time-axis that matches in length to 2nd dim of np_array")
            doc_parameter("geo_points","GeoPointVector","the geo-positions for the time-series, should be of length n_ts")
            doc_parameter("np_array","np.ndarray","numpy array of dtype=np.float64, and shape(n_ts,n_points)")
            doc_parameter("point_fx","point interpretation", "one of POINT_AVERAGE_VALUE|POINT_INSTANT_VALUE")
            doc_returns("tsv","TemperatureSourceVector","a TemperatureSourceVector of length first np_array dim, n_ts, each with geo-point and time-series with time-axis, values and point_fx")
        );

        def("create_precipitation_source_vector_from_np_array",
            &create_from_geo_tsv_from_np<sa::PrecipitationSource>,(py::arg("time_axis"),py::arg("geo_points"),py::arg("np_array"),py::arg("point_fx")),
            doc_intro("Create a PrecipitationSourceVector from specified time_axis,geo_points, 2-d np_array and point_fx.")
            doc_parameters()
            doc_parameter("time_axis","TimeAxis","time-axis that matches in length to 2nd dim of np_array")
            doc_parameter("geo_points","GeoPointVector","the geo-positions for the time-series, should be of length n_ts")
            doc_parameter("np_array","np.ndarray","numpy array of dtype=np.float64, and shape(n_ts,n_points)")
            doc_parameter("point_fx","point interpretation", "one of POINT_AVERAGE_VALUE|POINT_INSTANT_VALUE")
            doc_returns("tsv","PrecipitationSourceVector","a PrecipitationSourceVector of length first np_array dim, n_ts, each with geo-point and time-series with time-axis, values and point_fx")
        );

        def("create_wind_speed_source_vector_from_np_array",
            &create_from_geo_tsv_from_np<sa::WindSpeedSource>,(py::arg("time_axis"),py::arg("geo_points"),py::arg("np_array"),py::arg("point_fx")),
            doc_intro("Create a WindSpeedSourceVector from specified time_axis,geo_points, 2-d np_array and point_fx.")
            doc_parameters()
            doc_parameter("time_axis","TimeAxis","time-axis that matches in length to 2nd dim of np_array")
            doc_parameter("geo_points","GeoPointVector","the geo-positions for the time-series, should be of length n_ts")
            doc_parameter("np_array","np.ndarray","numpy array of dtype=np.float64, and shape(n_ts,n_points)")
            doc_parameter("point_fx","point interpretation", "one of POINT_AVERAGE_VALUE|POINT_INSTANT_VALUE")
            doc_returns("tsv","WindSpeedSourceVector","a WindSpeedSourceVector of length first np_array dim, n_ts, each with geo-point and time-series with time-axis, values and point_fx")
        );

        def("create_rel_hum_source_vector_from_np_array",
            &create_from_geo_tsv_from_np<sa::RelHumSource>,(py::arg("time_axis"),py::arg("geo_points"),py::arg("np_array"),py::arg("point_fx")),
            doc_intro("Create a RelHumSourceVector from specified time_axis,geo_points, 2-d np_array and point_fx.")
            doc_parameters()
            doc_parameter("time_axis","TimeAxis","time-axis that matches in length to 2nd dim of np_array")
            doc_parameter("geo_points","GeoPointVector","the geo-positions for the time-series, should be of length n_ts")
            doc_parameter("np_array","np.ndarray","numpy array of dtype=np.float64, and shape(n_ts,n_points)")
            doc_parameter("point_fx","point interpretation", "one of POINT_AVERAGE_VALUE|POINT_INSTANT_VALUE")
            doc_returns("tsv","RelHumSourceVector","a RelHumSourceVector of length first np_array dim, n_ts, each with geo-point and time-series with time-axis, values and point_fx")
        );

        def("create_radiation_source_vector_from_np_array",
            &create_from_geo_tsv_from_np<sa::RadiationSource>,(py::arg("time_axis"),py::arg("geo_points"),py::arg("np_array"),py::arg("point_fx")),
            doc_intro("Create a RadiationSourceVector from specified time_axis,geo_points, 2-d np_array and point_fx.")
            doc_parameters()
            doc_parameter("time_axis","TimeAxis","time-axis that matches in length to 2nd dim of np_array")
            doc_parameter("geo_points","GeoPointVector","the geo-positions for the time-series, should be of length n_ts")
            doc_parameter("np_array","np.ndarray","numpy array of dtype=np.float64, and shape(n_ts,n_points)")
            doc_parameter("point_fx","point interpretation", "one of POINT_AVERAGE_VALUE|POINT_INSTANT_VALUE")
            doc_returns("tsv","RadiationSourceVector","a RadiationSourceVector of length first np_array dim, n_ts, each with geo-point and time-series with time-axis, values and point_fx")
        );
        numpy_boost_python_register_type<double, 2>();
    }
	typedef std::vector<utctime> UtcTimeVector;
    extern utctime  x_kwarg_utctime(const py::tuple& args, const py::dict& kwargs,size_t i,const char *kw) ;
	struct utc_ext {
		static UtcTimeVector *create_default() {
			return new UtcTimeVector();
		}
		static UtcTimeVector *create_from_intv(const std::vector<int>&v) {
			auto r = new UtcTimeVector(); r->reserve(v.size());
			for (auto s : v)
				r->emplace_back(utctime{ std::chrono::seconds(s) });
			return r;
		}

		static UtcTimeVector *create_from_doublev(const std::vector<double>&v) {
			auto r = new UtcTimeVector(); r->reserve(v.size());
			for (auto s : v)
				r->emplace_back(from_seconds(s));
			return r;
		}
		template<class tp>
		static UtcTimeVector *create_from_np_tp(const numpy_boost<tp, 1>& npv) {
			auto r = new UtcTimeVector();
			size_t n = npv.shape()[0];
			r->reserve(n);
			for (size_t i = 0; i<n; ++i) {
				r->emplace_back(from_seconds(npv[i]));
			}
			return r;
		}
		static UtcTimeVector *create_from_clone(const std::vector<utctime>&v) {
			return new UtcTimeVector(v);
		}
		static UtcTimeVector *create_from_list(py::list times) {
			if(py::len(times)==0)
				return new UtcTimeVector();
			auto r = new UtcTimeVector();
			size_t n = py::len(times);
			r->reserve(py::len(times));
			for(size_t i=0;i<n;++i) {
				py::object oi = times[i];
				py::extract<utctime> xtract_utctime(oi);
				if (xtract_utctime.check()) {
					r->push_back(xtract_utctime());
				} else {
					py::extract<int64_t> xtract_int(oi);
					if (xtract_int.check()) {
						r->push_back(utctime{ std::chrono::seconds(xtract_int()) });
					} else {
						py::extract<double> xtract_double(oi);
						if (xtract_double.check()) {
							r->push_back(utctime{ from_seconds(xtract_double()) });
						} else {
							py::extract<string> xtract_string(oi);
							if (xtract_string.check()){
								r->push_back(create_from_iso8601_string(xtract_string()));
							} else {
								throw std::runtime_error(std::string("failed to convert ") + std::to_string(i) + " element to utctime");
							}
						}
					}
				}
			}
			return r;

		}

		static py::object push_back(const py::tuple& args,const py::dict& kwargs) {
            py::extract<vector<utctime>&> s(args[0]);
            if(!s.check()) throw std::runtime_error("UtcTimeVector: invalid self to push_back");
            auto& v=s();
            auto t=x_kwarg_utctime(args,kwargs,1,"t");
            v.push_back(t);
            return py::object();
        }

	};
	
    static numpy_boost<int64_t, 1> utctime_to_numpy(const vector<utctime>&v) {
		int dims[] = { int(v.size()) };
		numpy_boost<int64_t, 1> r(dims);
		for (size_t i = 0; i<r.size(); ++i) {
			r[i] = to_seconds64( v[i]);
		}
		return r;
	}
    static numpy_boost<double, 1> utctime_to_numpy_double(const vector<utctime>&v) {
		int dims[] = { int(v.size()) };
		numpy_boost<double, 1> r(dims);
		for (size_t i = 0; i<r.size(); ++i) {
			r[i] = to_seconds( v[i]);
		}
		return r;
	}
	
	static vector<utctime> utctime_from_numpy(const numpy_boost<int64_t, 1>& npv) {
		size_t n = npv.shape()[0];
		vector<utctime> r; r.reserve(n);
		for (size_t i = 0; i<npv.shape()[0]; ++i) {
			r.emplace_back(from_seconds(npv[i]));
		}
		return r;
	}

	static void expose_utctime_vector() {

		class_<UtcTimeVector>("UtcTimeVector",no_init)
			.def(vector_indexing_suite<UtcTimeVector>()) // meaning it get all it needs to appear as python list
			.def("__init__", make_constructor(&utc_ext::create_default,
				default_call_policies()
				),
				doc_intro("construct empty UtcTimeVecor")
			)
			.def("__init__", make_constructor(&utc_ext::create_from_clone,
				default_call_policies(),(py::arg("clone_me"))
				),
				doc_intro("construct a copy of supplied  UtcTimeVecor")
				doc_parameters()
				doc_parameter("clone_me", "UtcTimeVector", "to be cloned")
			)
			.def("__init__", make_constructor(&utc_ext::create_from_intv,
				default_call_policies(),
				(py::arg("seconds_vector"))
			    ),
				doc_intro("construct a from seconds epoch utc as integer")
				doc_parameters()
				doc_parameter("seconds", "IntVector", "seconds")
			)
			.def("__init__", make_constructor(&utc_ext::create_from_doublev,
				default_call_policies(),
				(py::arg("seconds_vector"))
			    ),
				doc_intro("construct a from seconds epoch utc as float")
				doc_parameters()
				doc_parameter("seconds", "DoubleVector", "seconds, up to us resolution epoch utc")
			)
			.def("__init__", make_constructor(&utc_ext::create_from_list,
				default_call_policies(),
				(py::arg("times"))
				),
				doc_intro("construct a from a list of something that is convertible to UtcTime")
				doc_parameters()
				doc_parameter("times", "list", "a list with convertible times")
			)
			.def("__init__", make_constructor(&utc_ext::create_from_np_tp<int64_t>,
				default_call_policies(),
				(py::arg("np_times"))
				),
				doc_intro("construct a from a numpy array of int64 s epoch")
				doc_parameters()
				doc_parameter("np_times", "list", "a list with convertible times")
			)
			.def("__init__", make_constructor(&utc_ext::create_from_np_tp<double>,
				default_call_policies(),
				(py::arg("np_times"))
			    ),
				doc_intro("construct a from a numpy array of float s epoch")
				doc_parameters()
				doc_parameter("np_times", "list", "a list with float convertible times")
			)
            .def("push_back",raw_function(utc_ext::push_back,2),
                doc_intro("appends a utctime like value to the vector")
                doc_parameters()
                doc_parameter("t","utctime","an int (seconds), or utctime")
            )
			.def("from_numpy", utctime_from_numpy).staticmethod("from_numpy") //bw.compatible
			.def("to_numpy", utctime_to_numpy,(py::arg("self")),
				doc_intro("convert to numpy array of type np.int64, seconds since epoch")
			)
			.def("to_numpy_double", utctime_to_numpy_double,(py::arg("self")),
				doc_intro("convert to numpy array of type np.float64, seconds since epoch")
			)
			;

		numpy_boost_python_register_type<utctime, 1>(); // register the numpy object so we can access it in C++
		numpy_boost_python_register_type<int64_t, 1>();
		py_api::iterable_converter().from_python<UtcTimeVector>();

	}
    void vectors() {
        np_import();
        expose_str_vector("StringVector");
        expose_vector<double>("DoubleVector");
		expose_vector<vector<double>>("DoubleVectorVector");
        expose_vector<int>("IntVector");
        expose_vector<int64_t>("Int64Vector");
        expose_vector<char>("ByteVector");
        //expose_vector<utctime>("UtcTimeVector");
        expose_utctime_vector();
        expose_geo_point_vector();
        expose_geo_cell_data_vector();
        expose_ts_vector_create();
    }
}

