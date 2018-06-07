/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#include "boostpython_pch.h"

#include "py_convertible.h"

#include "core/utctime_utilities.h"
#include "core/time_axis.h"
#include "core/time_series.h"
#include "core/cell_model.h"

namespace expose {
    using namespace boost::python;
    namespace py=boost::python;
    void cell_environment() {

        typedef shyft::core::environment_t CellEnvironment;
        class_<CellEnvironment>("CellEnvironment","Contains all ts projected to a certain cell-model using interpolation step (if needed)")
            .def_readwrite("temperature",&CellEnvironment::temperature)
            .def_readwrite("precipitation",&CellEnvironment::precipitation)
            .def_readwrite("radiation",&CellEnvironment::radiation)
            .def_readwrite("wind_speed",&CellEnvironment::wind_speed)
            .def_readwrite("rel_hum",&CellEnvironment::rel_hum)
            .def("init",&CellEnvironment::init,(py::arg("self"),py::arg("ta")),
                 "zero all series, set time-axis ta"
            )
            .def("has_nan_values",&CellEnvironment::has_nan_values,(py::arg("self")),
                doc_intro("scans all time-series for nan-values")
                doc_returns("has_nan","bool","true if any nan is encounted, otherwise false")
            )
            ;

        typedef shyft::core::environment_const_rhum_and_wind_t CellEnvironmentConstRHumWind;
        class_<CellEnvironmentConstRHumWind>("CellEnvironmentConstRHumWind","Contains all ts projected to a certain cell-model using interpolation step (if needed)")
            .def_readwrite("temperature",&CellEnvironmentConstRHumWind::temperature)
            .def_readwrite("precipitation",&CellEnvironmentConstRHumWind::precipitation)
            .def_readwrite("radiation",&CellEnvironmentConstRHumWind::radiation)
            .def_readwrite("wind_speed",&CellEnvironmentConstRHumWind::wind_speed)
            .def_readwrite("rel_hum",&CellEnvironmentConstRHumWind::rel_hum)
            .def("init",&CellEnvironmentConstRHumWind::init,args("ta"),"zero all series, set time-axis ta")
            ;
            
        enum_<shyft::core::stat_scope>("stat_scope",
            doc_intro("Defines the scope of the indexes passed to statistics functions")
            doc_intro(".cell      : the indexes are considered as the i'th cell")
            doc_intro(".catchment : the indexes are considered catchment identifiers")
            doc_intro("")
            doc_intro("The statistics is then covering the cells that matches the selected criteria")
            )
            .value("cell",shyft::core::stat_scope::cell_ix)
            .value("catchment",shyft::core::stat_scope::catchment_ix)
            .export_values()
            ;
    }
}
