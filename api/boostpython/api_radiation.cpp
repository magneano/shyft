/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#include "boostpython_pch.h"

#include "core/radiation.h"
#include <armadillo>

namespace expose {

    void radiation() {
        using namespace shyft::core::radiation;
        using namespace boost::python;
        namespace py = boost::python;
        typedef arma::vec surface_normal;

        class_<parameter>("RadiationParameter")
            .def(init<double,double>(args("albedo","turbidity"),"a new object with specified parameters"))
            .def_readwrite("albedo",&parameter::albedo,"typical value 0.2")
            .def_readwrite("turbidity", &parameter::turbidity,"typical value 1.0")
            ;

        class_<response>("RadiationResponse")
            .def_readwrite("dir_radiation",&response::dir_radiation)
            .def_readwrite("dif_radiation",&response::dif_radiation)
            .def_readwrite("ref_radiation",&response::ref_radiation)
            .def_readwrite("psw_radiation",&response::psw_radiation)
            .def_readwrite("tsw_radiation",&response::tsw_radiation)
            ;

        typedef calculator<parameter,response> RadiationCalculator;
        class_<RadiationCalculator>("RadiationCalculator",
                "Radiation,R, (ref.: Allen, R. G.; Trezza, R. & Tasumi, M. Analytical integrated functions for daily solar radiation on slopes Agricultural and Forest Meteorology, 2006, 139, 55-73)\n"
                "primitive implementation for calculating predicted clear-sky short-wave solar radiation for inclined surfaces\n"
                "This function is plain and simple, taking albedo and turbidity\n"
                "into the constructor and provides 2 functions:\n"
                " psw_radiation calculates predicted solar radiation (if no measured data available);\n"
                " tsw_radiation translates measured horizontal radiation into sloped surface\n"
                "[mm/s] units.\n",no_init
            )
            .def(init<const parameter&>(args("param"),"create a calculator using supplied parameter"))
            .def("psw_radiation",&RadiationCalculator::psw_radiation<surface_normal>,(py::arg("self"),py::arg("response"), py::arg("latitude"), py::arg("t"), py::arg("surface_normal"), py::arg("temperature"), py::arg("rhumidity"), py::arg("elevation")),
                 doc_intro("calculates predicted short-wave radiation, updating response")
            )
            .def("tsw_radiation",&RadiationCalculator::tsw_radiation<surface_normal>,(py::arg("self"),py::arg("response"), py::arg("latitude"), py::arg("t"), py::arg("surface_normal"), py::arg("temperature"), py::arg("rhumidity"), py::arg("elevation"),  py::arg("rsm")),
                     doc_intro("calculates predicted short-wave radiation, updating response")
            )
            ;
    }
}

BOOST_PYTHON_MODULE(_radiation)
{

    boost::python::scope().attr("__doc__")="Shyft python api for the radiation model";
    boost::python::docstring_options doc_options(true, true, false);// all except c++ signatures
    expose::radiation();
}
