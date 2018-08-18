/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#include "boostpython_pch.h"
#include "core/utctime_utilities.h"

/*
 * This section provides a  functional-convertible mechanism that was lacking in boost python.
 * The idea is that there is a function Fx that takes a type Source and construct inplace a type Target.
 *
 */
namespace boost { namespace python {

namespace converter {

template <class Source, class Target,class Fx>
struct fx_implicit
{
    static void* convertible(PyObject* obj)
    {
        // Find a converter which can produce a Source instance from
        // obj. The user has told us that Source can be converted to
        // Target, and instantiating construct() below, ensures that
        // at compile-time.
        return implicit_rvalue_convertible_from_python(obj, registered<Source>::converters)
            ? obj : 0;
    }

    static void construct(PyObject* obj, rvalue_from_python_stage1_data* data)
    {
        void* storage = ((rvalue_from_python_storage<Target>*)data)->storage.bytes;

        arg_from_python<Source> get_source(obj);
        bool convertible = get_source.convertible();
        BOOST_VERIFY(convertible);
        Fx::construct(get_source(),(Target*)storage);
        // record successful construction
        data->convertible = storage;
    }
};

} // namespace converter

template <class Source, class Target,class Fx>
void fx_implicitly_convertible(boost::type<Source>* = 0, boost::type<Target>* = 0)
{
    typedef converter::fx_implicit<Source,Target,Fx> functions;

    converter::registry::push_back(
          &functions::convertible
        , &functions::construct
        , type_id<Target>()
#ifndef BOOST_PYTHON_NO_PY_SIGNATURES
        , &converter::expected_from_python_type_direct<Source>::get_pytype
#endif
        );
}
}} // namespace boost::python

namespace expose {
    using namespace shyft::core;
    using namespace boost::python;
    namespace py=boost::python;
    using namespace std;

    struct ct_from_int64 {
        static void construct(int64_t src_seconds,utctime*dst) {
            *dst=utctime_from_seconds64(src_seconds);
        }
    };
    struct ct_from_utctime {
        static void construct(utctime src,int64_t*dst) {
            *dst=to_seconds64(src);
        }
    };
    struct ct_utctime_from_double {
        static void construct(double src_seconds,utctime*dst) {
            *dst=from_seconds(src_seconds);
        }
    };
    struct ct_double_from_utctime {
        static void construct(utctime src,double*dst) {
            *dst=to_seconds(src);
        }
    };

    static bool is_npos(size_t n) {
        return n==string::npos;
    }

	template<class T>
	static T x_arg(const py::tuple& args, size_t i) {
		if (py::len(args) + 1 < (int)i)
			throw std::runtime_error("missing arg #" + std::to_string(i) + std::string(" in time"));
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

	static void args_check(const py::tuple& args) {
		if (py::len(args) < 2) 	throw std::runtime_error("compare needs two args");
	}

	struct utctime_ext {
		static utctime* create_default() {
			return new utctime{ 0 };
		}
		static utctime* create_from_int(int sec) {
			return new utctime(seconds(sec));
		}
		static utctime* create_from_double(double sec) {
			return new utctime{ from_seconds(sec) };
		}
		// consider spirit for faster and accurate parsing !
		static utctime* create_from_string(const std::string& s) {
			return new utctime{ create_from_iso8601_string(s) };
		}

		static utctime x_self(const py::tuple& args) {
			if (py::len(args) == 0)
				throw std::runtime_error("self is null in time");
			py::object self = args[0];
			py::extract<utctime> xtract_self(self);
			return xtract_self();
		}
		static py::object get_seconds(py::tuple args, py::dict kwargs) {
			utctime dt = x_self(args);
			auto dt_s = std::chrono::duration_cast<std::chrono::seconds>(dt);
			if (dt_s == dt)
				return py::object(int64_t(dt_s.count()));
			return py::object(to_seconds(dt));
		}
		static py::object get_int_seconds(py::tuple args, py::dict kwargs) {
			utctime dt = x_self(args);
			auto dt_s = std::chrono::duration_cast<std::chrono::seconds>(dt);
			if (dt_s == dt)
				return py::object(int64_t(dt_s.count()));
			return py::object(int64_t(to_seconds(dt)));
		}

		static py::object str(py::tuple args, py::dict kwargs) {
		    // note: this is somewhat fuzzy logic to make nice- strings in the debugger...
		    // presenting it as a iso 8601 if > 1 year from 1970, otherwise just report seconds
			utctime dt = x_self(args);
			auto dt_s = std::chrono::duration_cast<std::chrono::seconds>(dt);
			if(dt > calendar::YEAR || dt < -calendar::YEAR) { // assume user want to see year..
                calendar utc;
                return py::str(utc.to_string(dt));
			} else {
                char s[100];
                if (dt_s == dt)
                    #ifndef _WIN32
                    sprintf(s, "%lds",int64_t(dt_s.count()));
                    #else
                    sprintf(s, "%llds", int64_t(dt_s.count()));
                    #endif
                else
                    sprintf(s, "%0.6lfs", to_seconds(dt));
                return py::str(std::string(s));
			}
		}
		static py::object repr(py::tuple args, py::dict kwargs) {
			utctime dt = x_self(args);
			auto dt_s = std::chrono::duration_cast<std::chrono::seconds>(dt);
			char s[100];
			if (dt_s == dt)
#ifndef _WIN32
				sprintf(s, "time(%ld)", dt_s.count());
#else
				sprintf(s, "time(%lld)", dt_s.count());
#endif
			else
				sprintf(s, "time(%0.6lf)", to_seconds(dt));
			return py::str(std::string(s));
		}

		static utctime abs_timespan(utctime x) {
			return x>=x.zero() ? x:-x;
		}
		static double _float_(utctime x) {
			return to_seconds(x);
		}
		static utctime _round_(utctime x) {
            return from_seconds(std::round(to_seconds(x)));
        }
		/// best effort convert from any number or string to utctime
		/// as long as it is well defined
		static utctime as_utctime(const py::object& po) {
			py::extract<utctime> o(po);
			if (o.check())
				return o();
			py::extract<int64_t> i(po);
			if (i.check())
				return  from_seconds(i());
			py::extract<double> d(po);
			if (d.check())
				return   from_seconds(d());
            py::extract<string> s(po);
			if (s.check())
				return create_from_iso8601_string(s);
			throw std::runtime_error("supplied argument not convertible to time");
		}
		// rel ops
		static py::object _lt_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctime>(args, 0) < as_utctime(args[1]));}
		static py::object _le_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctime>(args, 0) <= as_utctime(args[1]));}
		static py::object _gt_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctime>(args, 0) > as_utctime(args[1]));}
		static py::object _ge_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctime>(args, 0) >= as_utctime(args[1]));}
		static py::object _eq_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctime>(args, 0) == as_utctime(args[1]));}
		static py::object _nq_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctime>(args, 0) != as_utctime(args[1]));}
		// math
		static py::object _add_(py::tuple args, py::dict kwargs) { args_check(args); return py::object(x_arg<utctime>(args, 0) + as_utctime(args[1])); }
		static py::object _sub_(py::tuple args, py::dict kwargs) { args_check(args); return py::object(x_arg<utctime>(args, 0) - as_utctime(args[1])); }
		static py::object _rsub_(py::tuple args, py::dict kwargs) { args_check(args); return py::object(-(x_arg<utctime>(args, 0) - as_utctime(args[1]))); }
		static py::object _mult_(py::tuple args, py::dict kwargs) {	args_check(args); return py::object(from_seconds( to_seconds(x_arg<utctime>(args, 0)) * to_seconds(as_utctime(args[1])) ));}
		static py::object _div_(py::tuple args, py::dict kwargs) { args_check(args); return py::object(from_seconds( to_seconds(x_arg<utctime>(args, 0)) / to_seconds(as_utctime(args[1])) ));	}
		static py::object _mod_(py::tuple args, py::dict kwargs) {	args_check(args);return py::object( x_arg<utctime>(args, 0) % as_utctime(args[1]) );}
		static py::object _floordiv_(py::tuple args, py::dict kwargs) {	args_check(args);return py::object(from_seconds(truncf(to_seconds(x_arg<utctime>(args, 0)) / to_seconds(as_utctime(args[1])))));}
		static py::object _rfloordiv_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(from_seconds(truncf(to_seconds(as_utctime(args[1]))/to_seconds(x_arg<utctime>(args, 0)) )));}
		static py::object _rdiv_(py::tuple args, py::dict kwargs) {	args_check(args);return py::object(from_seconds( to_seconds(as_utctime(args[1]))/ to_seconds(x_arg<utctime>(args, 0))));}
		// just for backward(consider change code)
		static py::object _sqrt_(py::tuple args, py::dict kwargs) {return py::object(from_seconds(sqrt(to_seconds(as_utctime(args[0])) ))); }
	};
	utctime  x_kwarg_utctime(const py::tuple& args, const py::dict& kwargs,size_t i,const char *kw) {
        if (py::len(args)  > (int)i) {
            return utctime_ext::as_utctime(args[i]);
        }
        if(kwargs.has_key(kw)) {
            return utctime_ext::as_utctime(kwargs[kw]);
        }
        throw std::runtime_error("missing kw arg #" + std::string(kw) );
	}

    struct utctime_picklers : py::pickle_suite {
		static py::tuple getinitargs(utctime const& t) { return py::make_tuple(to_seconds(t)); }
	};

    struct calendar_ext {
        static py::object trim(const py::tuple& args,const py::dict & kwargs) {
            const auto& c=x_arg<const calendar&>(args,0);
            auto t= x_kwarg_utctime(args,kwargs,1,"t");
            auto dt= x_kwarg_utctime(args,kwargs,2,"delta_t");
            return py::object(c.trim(t,dt));
        }
        static py::object add(const py::tuple& args,const py::dict & kwargs) {
            const auto& c=x_arg<const calendar&>(args,0);
            auto t= x_kwarg_utctime(args,kwargs,1,"t");
            auto dt= x_kwarg_utctime(args,kwargs,2,"delta_t");
            auto  n = x_kwarg<int>(args,kwargs,3,"n");
            return py::object(c.add(t,dt,n));
        }
        static py::object diff_units(const py::tuple& args,const py::dict & kwargs) {
            const auto& c=x_arg<const calendar&>(args,0);
            auto t1= x_kwarg_utctime(args,kwargs,1,"t1");
            auto t2= x_kwarg_utctime(args,kwargs,2,"t2");
            auto dt= x_kwarg_utctime(args,kwargs,3,"delta_t");
            return py::object(c.diff_units(t1,t2,dt));
        }
        static py::object to_string(const py::tuple& args,const py::dict & kwargs) {
            const auto& c=x_arg<const calendar&>(args,0);
            return py::object(c.to_string(x_kwarg_utctime(args,kwargs,1,"t")));
        }
        static py::object calendar_units(const py::tuple& args,const py::dict & kwargs) {
            const auto& c=x_arg<const calendar&>(args,0);
            return py::object(c.calendar_units(x_kwarg_utctime(args,kwargs,1,"t")));
        }
        static py::object calendar_week_units(const py::tuple& args,const py::dict & kwargs) {
            const auto& c=x_arg<const calendar&>(args,0);
            return py::object(c.calendar_week_units(x_kwarg_utctime(args,kwargs,1,"t")));
        }
    };

	static void e_utctime() {
		class_<utctime>("time",
              doc_intro("time is represented as a number, in unit seconds.")
              doc_intro("For accuracy and performance, it's internally represented as 64bit integer, at micro-second resolution")
              doc_intro("It is usually used in two roles:")
              doc_intro("  1) time measured in seconds since epoch (1970.01.01UTC)")
              doc_intro("      often constructed using the Calendar class that takes calendar-coordinates (YMDhms etc) and ")
              doc_intro("      returns the corresponding time-point, taking time-zone and dst etc. into account.")
              doc_intro("      >>> utc = Calendar()")
              doc_intro("      >>> t1 = utc.time(2018,10,15,16,30,15)")
              doc_intro("      >>> t2 = time('2018-10-15T16:30:15Z')")
              doc_intro("      >>> t3 = time(1539621015)")
              doc_intro("  2) time measure, in unit of seconds (resolution up to 1us)")
              doc_intro("      often constructed from numbers, in SI-unit of seconds")
              doc_intro("      >>> dt1 = time(3600)")
              doc_intro("      >>> dt2 = time(1.5)")
              doc_intro("  ")
              doc_intro("It can be constructed supplying number of seconds, or a well defined iso8601 string\n")
              doc_intro("To convert it to a python number, use float() or int() to cast operations")
              doc_intro("If dealing with time-zones/calendars conversion, use the calendar.time(..)/.calendar_units functions.")
              doc_intro("If you want to use time-zone/calendar semantic add/diff/trim, use the corresponding Calendar methods.")
              doc_intro("")
              doc_see_also("Calendar,deltahours,deltaminutes")
              ,no_init
            )
			.def("__init__", make_constructor(&utctime_ext::create_default,
					default_call_policies()
				),
				doc_intro("construct a 0s")
			)
			.def("__init__", make_constructor(&utctime_ext::create_from_int,
				default_call_policies(),
				(py::arg("seconds"))
				),
				doc_intro("construct a from seconds as integer")
				doc_parameters()
				doc_parameter("seconds", "int", "seconds")
			)
			.def("__init__", make_constructor(&utctime_ext::create_from_double,
				default_call_policies(),
				(py::arg("seconds"))
				),
				doc_intro("construct a from seconds as decimal number, where fractions is fractions of second")
				doc_intro("- the resulting time preserves 1 micro-second digits")
				doc_parameters()
				doc_parameter("seconds", "float", "seconds")
			)
			.def("__init__", make_constructor(&utctime_ext::create_from_string,
				default_call_policies(),
				(py::arg("iso8601"))
				),
				doc_intro("construct a from iso 8601s string, e.g. '2014-11-12T19:12:14.505Z'")
				doc_intro("- the resulting time-span preserves 1 micro-second digits")
				doc_parameters()
				doc_parameter("iso8601", "str", "iso-formattet well defined time")
			)
			.def_pickle(utctime_picklers())
			.add_property("seconds",raw_function(utctime_ext::get_seconds,1),doc_intro("returns time in seconds"))
			.def("__abs__",&utctime_ext::abs_timespan,(py::arg("self")))
            .def("__round__",&utctime_ext::_round_,(py::arg("self")))
            .def("__float__",&utctime_ext::_float_,(py::arg("self")))
			.def("__int__",raw_function(utctime_ext::get_int_seconds,1),doc_intro("time as int seconds"))
			.def("__long__", raw_function(utctime_ext::get_int_seconds, 1), doc_intro("time as int seconds"))
			.def("__repr__",raw_function(utctime_ext::repr,1),doc_intro("repr of time"))
			.def("__str__", raw_function(utctime_ext::str, 1), doc_intro("str of time"))
			// rel operations
			.def("__lt__", raw_function(utctime_ext::_lt_, 1))
			.def("__le__", raw_function(utctime_ext::_le_, 1))
			.def("__gt__", raw_function(utctime_ext::_gt_, 1))
			.def("__ge__", raw_function(utctime_ext::_ge_, 1))
			.def("__eq__", raw_function(utctime_ext::_eq_, 1))
			.def("__nq__", raw_function(utctime_ext::_nq_, 1))

			// math
			.def("__add__",raw_function(utctime_ext::_add_, 1))
			.def("__radd__", raw_function(utctime_ext::_add_, 1))
			.def("__sub__", raw_function(utctime_ext::_sub_, 1))
			.def("__rsub__", raw_function(utctime_ext::_rsub_, 1))
			.def("__truediv__",raw_function(utctime_ext::_div_, 1))
			.def("__rtruediv__", raw_function(utctime_ext::_rdiv_, 1))
			.def("__floordiv__", raw_function(utctime_ext::_floordiv_, 1))
			.def("__rfloordiv__", raw_function(utctime_ext::_rfloordiv_, 1))
			.def("__mod__", raw_function(utctime_ext::_mod_, 1))

			.def("__mul__", raw_function(utctime_ext::_mult_, 1))
			.def("__rmul__", raw_function(utctime_ext::_mult_, 1))

			.def("sqrt",raw_function(utctime_ext::_sqrt_,1))
			//.def(self % self)
			.def(-self)
			;
			def("utctime_now",utctime_now,"returns time now as seconds since 1970s");
            def("deltahours",deltahours,args("n"),"returns time equal to specified n hours");
            def("deltaminutes",deltaminutes,args("n"),"returns time equal to specified n minutes");
            def("is_npos",is_npos,args("n"),"returns true if n is npos, - meaning no position");
            scope current;
            current.attr("max_utctime")= max_utctime;
            current.attr("min_utctime")= min_utctime;
            current.attr("no_utctime")=no_utctime;
            current.attr("npos")=string::npos;
            //TODO: require own type for utctime:
            //implicitly_convertible<utctime,int64_t>();
            //implicitly_convertible<int64_t,utctime>();
            //fx_implicitly_convertible<int64_t,utctime,ct_from_int64>();
            //fx_implicitly_convertible<utctime,int64_t,ct_from_utctime>();
            fx_implicitly_convertible<double,utctime,ct_utctime_from_double>();
            fx_implicitly_convertible<utctime,double,ct_double_from_utctime>();

	}

    typedef std::vector<utcperiod> UtcPeriodVector;


    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(calendar_time_overloads,calendar::time,1,7);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(calendar_time_overloads_week, calendar::time_from_week, 1, 7);

    static void e_calendar() {

        //std::string (shyft::core::calendar::*to_string_t)(shyft::core::utctime) const= &calendar::to_string;//selects correct ptr.
        std::string (calendar::*to_string_p)(utcperiod) const=&calendar::to_string;
        //int64_t (calendar::*diff_units)(utctime,utctime,utctimespan) const=&calendar::diff_units;
        utctime (calendar::*time_YMDhms)(YMDhms) const = &calendar::time;
        utctime (calendar::*time_YWdhms)(YWdhms) const = &calendar::time;
        utctime (calendar::*time_6)(int,int,int,int,int,int,int) const = &calendar::time;
        utctime (calendar::*time_from_week_6)(int, int, int, int, int, int,int) const = &calendar::time_from_week;
        class_<calendar, shared_ptr<calendar>>("Calendar",
            doc_intro("Calendar deals with the concept of human calendar")
            doc_intro(" In Shyft we practice the 'utc-time perimeter' principle,")
            doc_intro("  * so the core is utc-time only ")
            doc_intro("  * we deal with time-zone and calendars at the interfaces/perimeters")
            doc_intro(" in python, this corresponds to timestamp[64], or as the integer version of the time package representation")
            doc_intro(" e.g. the difference between time.time()- utctime_now() is in split-seconds")
            doc_intro("\n")
            doc_intro("Calendar functionality:")
            doc_intro(" -# Conversion between the calendar coordinates YMDhms or iso week YWdhms and utctime, taking  any timezone and DST into account\n")
            doc_intro(" -# Calendar constants, utctimespan like values for Year,Month,Week,Day,Hour,Minute,Second\n")
            doc_intro(" -# Calendar arithmetic, like adding calendar units, e.g. day,month,year etc.\n")
            doc_intro(" -# Calendar arithmetic, like trim/truncate a utctime down to nearest timespan/calendar unit. eg. day\n")
            doc_intro(" -# Calendar arithmetic, like calculate difference in calendar units(e.g days) between two utctime points\n")
            doc_intro(" -# Calendar Timezone and DST handling\n")
            doc_intro(" -# Converting time to string and vice-versa\n")
            doc_intro("\n")
            doc_notes()
            doc_note(" Please notice that although the calendar concept is complete")
            doc_note(" we only implement features as needed in the core and interfaces")
            doc_note(" Currently this includes most options, including olson time-zone handling")
            doc_note(" The time-zone support is currently a snapshot of rules ~2014")
            doc_note(" but we plan to use standard packages like Howard Hinnant's online approach for this later.")
            )
            .def(init<utctime>(args("tz_offset"),
                doc_intro("creates a calendar with constant tz-offset")
                doc_parameters()
                doc_parameter("tz_offset","time","specifies utc offset, time(3600) gives UTC+01 zone")
              )
            )
            .def(init<int>(args("tz_offset"),
                doc_intro("creates a calendar with constant tz-offset")
                doc_parameters()
                doc_parameter("tz_offset","int","seconds utc offset, 3600 gives UTC+01 zone")
              )
            )
            .def(init<string>(args("olson_tz_id"),
                doc_intro("create a Calendar from Olson timezone id, eg. 'Europe/Oslo'")
                doc_parameters()
                doc_parameter("olson_tz_id","str","Olson time-zone id, e.g. 'Europe/Oslo'")
                )
            )
            .def("to_string", raw_function(calendar_ext::to_string,2),//) to_string_t, (py::arg("self"),py::arg("t")),
                doc_intro("convert time t to readable iso standard string taking ")
                doc_intro(" the current calendar properties, including timezone into account")
                doc_parameters()
                doc_parameter("utctime","time","seconds utc since epoch")
                doc_returns("iso time string","str","iso standard formatted string,including tz info")
            )
            .def("to_string", to_string_p, (py::arg("self"),py::arg("utcperiod")),
                doc_intro("convert utcperiod p to readable string taking current calendar properties, including timezone into account")
                doc_parameters()
                doc_parameter("utcperiod", "UtcPeriod", "An UtcPeriod object")
                doc_returns("period-string", "str", "[start..end>, iso standard formatted string,including tz info")
            )
            .def("region_id_list", &calendar::region_id_list,
                doc_intro("Returns a list over predefined Olson time-zone identifiers")
                doc_notes()
                doc_note("the list is currently static and reflects tz-rules approx as of 2014")
            ).staticmethod("region_id_list")
            .def("calendar_units", raw_function(calendar_ext::calendar_units,2),
                doc_intro("returns YMDhms for specified t, in the time-zone as given by the calendar")
                doc_parameters()
                doc_parameter("t", "time", "timestamp utc seconds since epoch")
                doc_returns("calendar_units","YMDhms","calendar units as in year-month-day hour-minute-second")
            )
            .def("calendar_week_units", raw_function(calendar_ext::calendar_week_units,2),
                doc_intro("returns iso YWdhms for specified t, in the time-zone as given by the calendar")
                doc_parameters()
                doc_parameter("t", "time", "timestamp utc seconds since epoch")
                doc_returns("calendar_week_units", "YWdms", "calendar units as in iso year-week-week_day hour-minute-second")
            )
            .def("time", time_YMDhms,(py::arg("self"),py::arg("YMDhms")),
                doc_intro("convert calendar coordinates into time using the calendar time-zone")
                doc_parameters()
                doc_parameter("YMDhms","YMDhms","calendar cooordinate structure containg year,month,day, hour,minute,second")
                doc_returns("timestamp","int","timestamp as in seconds utc since epoch")
            )
            .def("time", time_YWdhms, (py::arg("self"),py::arg("YWdhms")),
                doc_intro("convert calendar iso week coordinates structure into time using the calendar time-zone")
                doc_parameters()
                doc_parameter("YWdhms", "YWdhms", "structure containg iso specification calendar coordinates")
                doc_returns("timestamp", "int", "timestamp as in seconds utc since epoch")
            )

            .def("time", time_6, calendar_time_overloads(
                doc_intro("convert calendar coordinates into time using the calendar time-zone")
                doc_parameters()
                doc_parameter("Y", "int", "Year")
                doc_parameter("M", "int", "Month  [1..12], default=1")
                doc_parameter("D", "int", "Day    [1..31], default=1")
                doc_parameter("h", "int", "hour   [0..23], default=0")
                doc_parameter("m", "int", "minute [0..59], default=0")
                doc_parameter("s", "int", "second [0..59], default=0")
                doc_parameter("us","int","micro second[0..999999], default=0")
                doc_returns("timestamp", "time", "timestamp as in seconds utc since epoch")
                ,
                (py::arg("self"),py::arg("Y"),py::arg( "M"),py::arg( "D"),py::arg( "h"),py::arg( "m"),py::arg( "s"),py::arg("us"))
            )
            )
            .def("time_from_week", time_from_week_6,
                calendar_time_overloads_week(
                    doc_intro("convert calendar iso week coordinates into time using the calendar time-zone")
                    doc_parameters()
                    doc_parameter("Y", "int", "ISO Year")
                    doc_parameter("W", "int", "ISO Week  [1..54], default=1")
                    doc_parameter("wd","int", "week_day  [1..7]=[mo..su], default=1")
                    doc_parameter("h", "int", "hour   [0..23], default=0")
                    doc_parameter("m", "int", "minute [0..59], default=0")
                    doc_parameter("s", "int", "second [0..59], default=0")
                    doc_parameter("us","int","micro second[0..999999], default=0")
                    doc_returns("timestamp", "time", "timestamp as in seconds utc since epoch")
                    ,
                    ( py::arg("self"), py::arg("Y"),py::arg("W"),py::arg("wd"),py::arg("h"),py::arg("m"),py::arg("s"),py::arg("us")  )
                )
            )

            .def("add", raw_function(calendar_ext::add,4), // args("t", "delta_t", "n"),
                doc_intro("calendar semantic add")
                doc_intro("conceptually this is similar to t + deltaT*n")
                doc_intro(" but with deltaT equal to calendar::DAY,WEEK,MONTH,YEAR")
                doc_intro(" and/or with dst enabled time-zone the variation of length due to dst")
                doc_intro(" or month/year length is taken into account")
                doc_intro(" e.g. add one day, and calendar have dst, could give 23,24 or 25 hours due to dst.")
                doc_intro(" similar for week or any other time steps.")
                doc_parameters()
                doc_parameter("t","time","timestamp utc seconds since epoch")
                doc_parameter("delta_t","time","timestep in seconds, with semantic interpretation of DAY,WEEK,MONTH,YEAR")
                doc_parameter("n","int","number of timesteps to add")
                doc_returns("t","time","new timestamp with the added time-steps, seconds utc since epoch")
                doc_see_also("diff_units(t1,t2,delta_t),trim(t,delta_t)")
            )
        .def("diff_units",raw_function(calendar_ext::diff_units,4),//diff_units,args("t1","t2","delta_t"),
             doc_intro("calculate the distance t1..t2 in specified units, taking dst into account if observed")
             doc_intro("The function takes calendar semantics when delta_t is calendar::DAY,WEEK,MONTH,YEAR,")
             doc_intro("and in addition also dst if observed.")
             doc_intro("e.g. the diff_units of calendar::DAY over summer->winter shift is 1,")
             doc_intro("even if the number of hours during those days are 23 and 25 summer and winter transition respectively")
             doc_intro("returns: calendar semantics of (t2-t1)/deltaT, where deltaT could be calendar units DAY,WEEK,MONTH,YEAR")
             doc_parameters()
             doc_parameter("t", "time", "timestamp utc seconds since epoch")
             doc_parameter("delta_t", "time", "timestep in seconds, with semantic interpretation of DAY,WEEK,MONTH,YEAR")
             doc_parameter("n", "int", "number of timesteps to add")
             doc_returns("n_units", "int", "number of units, so that t2 = c.add(t1,delta_t,n) + remainder(discarded)")
             doc_see_also("add(t,delta_t,n),trim(t,delta_t)")
        )
        .def("trim",raw_function(calendar_ext::trim,3),
            doc_intro("round time t down to the nearest calendar time-unit delta_t")
            doc_intro("taking the calendar time-zone and dst into account")
            doc_parameter("t", "time", "timestamp utc seconds since epoch")
            doc_parameter("delta_t", "time", "timestep in seconds, with semantic interpretation of Calendar.(DAY,WEEK,MONTH,YEAR)")
            doc_returns("t", "time", "new trimmed timestamp, seconds utc since epoch")
            doc_see_also("add(t,delta_t,n),diff_units(t1,t2,delta_t)")
        )

        .def("quarter",&calendar::quarter,args("t"),
            doc_intro("returns the quarter of the specified t, -1 if invalid t")
            doc_parameters()
            doc_parameter("t", "int", "timestamp utc seconds since epoch")
            doc_returns("quarter","int","in range[1..4], -1 if invalid time")
        )
        .def_readonly("YEAR",&calendar::YEAR, "defines a semantic year")
        .def_readonly("MONTH",&calendar::MONTH,"defines a semantic calendar month")
        .def_readonly("QUARTER",&calendar::QUARTER,"defines a semantic calendar quarter (3 months)")
        .def_readonly("DAY",&calendar::DAY,"defines a semantic calendar day")
        .def_readonly("WEEK",&calendar::WEEK,"defines a semantic calendar week")
        .def_readonly("HOUR",&calendar::HOUR,"hour, 3600 seconds")
        .def_readonly("MINUTE",&calendar::MINUTE,"minute, 60 seconds")
        .def_readonly("SECOND",&calendar::SECOND,"second, 1 second")
        .add_property("tz_info",&calendar::get_tz_info,"The TzInfo keeping the time-zone name, utc-offset and DST rules (if any)")//,return_value_policy<return_internal_reference>())
        ;

        class_<YMDhms>("YMDhms","Defines calendar coordinates as Year Month Day hour minute second")
        .def(init<int,py::optional<int,int,int,int,int,int>>( args("Y","M","D","h","m","s","us"),"Creates calendar coordinates specifying Y,M,D,h,m,s,us"))
        .def("is_valid",&YMDhms::is_valid,"returns true if YMDhms values are reasonable")
        .def("is_null",&YMDhms::is_null,"returns true if all values are 0, - the null definition")
        .def_readwrite("year",&YMDhms::year)
        .def_readwrite("month",&YMDhms::month)
        .def_readwrite("day",&YMDhms::day)
        .def_readwrite("hour",&YMDhms::hour)
        .def_readwrite("minute",&YMDhms::minute)
        .def_readwrite("second",&YMDhms::second)
        .def_readwrite("micro_second",&YMDhms::micro_second)
        .def("max",&YMDhms::max,"returns the maximum representation").staticmethod("max")
        .def("min",&YMDhms::max,"returns the minimum representation").staticmethod("min")
        .def(self==self)
        .def(self!=self)
           ;

        class_<YWdhms>("YWdhms", "Defines calendar coordinates as iso Year Week week-day hour minute second")
            .def(init<int, py::optional<int, int, int, int, int,int>>(args("Y", "W", "wd", "h", "m", "s","us"),
                doc_intro("Creates calendar coordinates specifying iso Y,W,wd,h,m,s")
                doc_parameters()
                doc_parameter("Y","int","iso-year")
                doc_parameter("W","int","iso week [1..53]")
                doc_parameter("wd","int","week_day [1..7]=[mo..sun]")
                doc_parameter("h","int","hour [0..23]")
                doc_parameter("m","int","minute [0..59]")
                doc_parameter("s","int","second [0..59]")
                doc_parameter("us","int","micro_second [0..999999]")
                )
            )
            .def("is_valid", &YWdhms::is_valid, "returns true if YWdhms values are reasonable")
            .def("is_null", &YWdhms::is_null, "returns true if all values are 0, - the null definition")
            .def_readwrite("iso_year", &YWdhms::iso_year)
            .def_readwrite("iso_week", &YWdhms::iso_week)
            .def_readwrite("week_day", &YWdhms::week_day,doc_intro("week_day,[1..7]=[mo..sun]"))
            .def_readwrite("hour", &YWdhms::hour)
            .def_readwrite("minute", &YWdhms::minute)
            .def_readwrite("second", &YWdhms::second)
            .def_readwrite("micro_second", &YWdhms::micro_second)
            .def("max", &YWdhms::max, "returns the maximum representation").staticmethod("max")
            .def("min", &YWdhms::max, "returns the minimum representation").staticmethod("min")
            .def(self == self)
            .def(self != self)
            ;

        class_<time_zone::tz_info_t,bases<>,time_zone::tz_info_t_,boost::noncopyable>("TzInfo",
            "TzInfo class is responsible for providing information about the\n"
            " time-zone of the calendar.\n"
            "  This include the\n"
            "   * name (olson identifier),\n"
            "   * base_offset\n"
            "   * utc_offset(t) time-dependent\n"
            "The Calendar class provides a shared pointer to it's TzInfo object \n",no_init
           )
        .def(init<utctimespan>(args("base_tz"),"creates a TzInfo with a fixed utc-offset(no dst-rules)"))
        .def("name",&time_zone::tz_info_t::name,"returns the olson time-zone identifier or name for the TzInfo")
        .def("base_offset",&time_zone::tz_info_t::base_offset,"returnes the time-invariant part of the utc-offset")
        .def("utc_offset",&time_zone::tz_info_t::utc_offset,args("t"),"returns the utc_offset at specified utc-time, takes DST into account if applicable")
        .def("is_dst",&time_zone::tz_info_t::is_dst,args("t"),"returns true if DST is observed at given utc-time t")
        ;
    }


    static void e_utcperiod() {
        enum_<shyft::core::trim_policy>("trim_policy",
                doc_intro("Enum to decide if to trim inwards or outwards where TRIM_IN means inwards")
            )
            .value("TRIM_IN", shyft::core::trim_policy::TRIM_IN)
            .value("TRIM_OUT", shyft::core::trim_policy::TRIM_OUT)
            .export_values()
            ;
        bool (utcperiod::*contains_t)(utctime) const = &utcperiod::contains;
        bool (utcperiod::*contains_i)(int64_t) const = &utcperiod::contains;
        bool (utcperiod::*contains_p)(const utcperiod&) const = &utcperiod::contains;
        utcperiod (utcperiod::*trim_i)(const calendar &, int64_t , trim_policy ) const=&utcperiod::trim;
        utcperiod (utcperiod::*trim_t)(const calendar &, utctimespan , trim_policy ) const=&utcperiod::trim;
        int64_t (utcperiod::*diff_units_t)(const calendar&, utctimespan) const = &utcperiod::diff_units;
        int64_t (utcperiod::*diff_units_i)(const calendar&, int64_t) const = &utcperiod::diff_units;

        class_<utcperiod>("UtcPeriod","UtcPeriod defines the open utc-time range [start..end> \nwhere end is required to be equal or greater than start")
        .def(init<utctime,utctime>(args("start,end"),"Create utcperiod given start and end"))
        .def(init<int64_t, int64_t>(args("start,end"), "Create utcperiod given start and end"))
        .def("valid",&utcperiod::valid,"returns true if start<=end otherwise false")
        .def("contains",contains_t,(py::arg("self"),py::arg("t")),"returns true if time t is contained in this utcperiod" )
        .def("contains",contains_i,(py::arg("self"),py::arg("t")),"returns true if time t is contained in this utcperiod" )
        .def("contains",contains_p,args("p"),"returns true if utcperiod p is contained in this utcperiod" )
        .def("overlaps",&utcperiod::overlaps,args("p"), "returns true if period p overlaps this utcperiod" )
        .def("trim",trim_t,(py::arg("self"),py::arg("calendar"),py::arg("delta_t"),py::arg("trim_policy")=shyft::core::trim_policy::TRIM_IN),
                doc_intro("Round UtcPeriod up or down to the nearest calendar time-unit delta_t")
                doc_intro("taking the calendar time-zone and dst into account")
                doc_parameters()
                doc_parameter("calendar", "calendar", "shyft calendar")
                doc_parameter("delta_t", "time", "timestep in seconds, with semantic interpretation of Calendar.(DAY,WEEK,MONTH,YEAR)")
                doc_parameter("trim_policy", "trim_policy", "TRIM_IN if rounding period inwards, else rounding outwards")
                doc_returns("trimmed_UtcPeriod", "UtcPeriod", "new trimmed UtcPeriod")
        )
        .def("trim",trim_i,(py::arg("self"),py::arg("calendar"),py::arg("delta_t"),py::arg("trim_policy")=shyft::core::trim_policy::TRIM_IN),
                doc_intro("Round UtcPeriod up or down to the nearest calendar time-unit delta_t")
                doc_intro("taking the calendar time-zone and dst into account")
                doc_parameters()
                doc_parameter("calendar", "calendar", "shyft calendar")
                doc_parameter("delta_t", "int", "timestep in seconds, with semantic interpretation of Calendar.(DAY,WEEK,MONTH,YEAR)")
                doc_parameter("trim_policy", "trim_policy", "TRIM_IN if rounding period inwards, else rounding outwards")
                doc_returns("trimmed_UtcPeriod", "UtcPeriod", "new trimmed UtcPeriod")
        )
        .def("diff_units",diff_units_t,(py::arg("self"),py::arg("calendar"),py::arg("delta_t")),
             doc_intro("Calculate the distance from start to end of UtcPeriod in specified units, taking dst into account if observed")
             doc_intro("The function takes calendar semantics when delta_t is calendar::DAY,WEEK,MONTH,YEAR,")
             doc_intro("and in addition also dst if observed.")
             doc_intro("e.g. the diff_units of calendar::DAY over summer->winter shift is 1,")
             doc_intro("even if the number of hours during those days are 23 and 25 summer and winter transition respectively")
             doc_parameters()
             doc_parameter("calendar", "calendar", "shyft calendar")
             doc_parameter("delta_t", "time", "timestep in seconds, with semantic interpretation of DAY,WEEK,MONTH,YEAR")
             doc_returns("n_units", "int", "number of units in UtcPeriod")
        )
        .def("diff_units",diff_units_i,(py::arg("self"),py::arg("calendar"),py::arg("delta_t")),
             doc_intro("Calculate the distance from start to end of UtcPeriod in specified units, taking dst into account if observed")
             doc_intro("The function takes calendar semantics when delta_t is calendar::DAY,WEEK,MONTH,YEAR,")
             doc_intro("and in addition also dst if observed.")
             doc_intro("e.g. the diff_units of calendar::DAY over summer->winter shift is 1,")
             doc_intro("even if the number of hours during those days are 23 and 25 summer and winter transition respectively")
             doc_parameters()
             doc_parameter("calendar", "calendar", "shyft calendar")
             doc_parameter("delta_t", "int", "timestep in seconds, with semantic interpretation of DAY,WEEK,MONTH,YEAR")
             doc_returns("n_units", "int", "number of units in UtcPeriod")
        )
        .def("__str__",&utcperiod::to_string,"returns the str using time-zone utc to convert to readable time")
        .def(self == self)
        .def(self != self)
        .def("timespan",&utcperiod::timespan,"returns end-start, the timespan of the period")
        .def_readwrite("start",&utcperiod::start,"Defines the start of the period, inclusive")
        .def_readwrite("end",&utcperiod::end,"Defines the end of the period, not inclusive")
        .def("intersection",&intersection,(py::arg("a"),py::arg("b")),
            doc_intro("Returns the intersection of two periods")
            doc_intro("if there is an intersection, the resulting period will be .valid() and .timespan()>0")
            doc_intro("If there is no intersection, an empty not .valid() period is returned")
            doc_parameters()
            doc_parameter("a","UtcPeriod","1st UtcPeriod argument")
            doc_parameter("b","UtcPeriod","2nd UtcPeriod argument")
            doc_returns("intersection","UtcPeriod","The computed intersection, or an empty not .valid() UtcPeriod")
        ).staticmethod("intersection");
        def("intersection",&intersection,(py::arg("a"),py::arg("b")),
            doc_intro("Returns the intersection of two periods")
            doc_intro("if there is an intersection, the resulting period will be .valid() and .timespan()>0")
            doc_intro("If there is no intersection, an empty not .valid() period is returned")
            doc_parameters()
            doc_parameter("a","UtcPeriod","1st UtcPeriod argument")
            doc_parameter("b","UtcPeriod","2nd UtcPeriod argument")
            doc_returns("intersection","UtcPeriod","The computed intersection, or an empty not .valid() UtcPeriod")
        );
    }
    void calendar_and_time() {
        e_utctime();
        e_utcperiod();
        e_calendar();
    }
}
