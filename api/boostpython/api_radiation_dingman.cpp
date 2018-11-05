/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#include "boostpython_pch.h"

#include "core/radiation_dingman.h"

namespace expose {

    void radiation_dingman() {
        using namespace shyft::core::radiation_dingman;
        using namespace boost::python;

        class_<parameter>("RadiationDingmanParameter")
            .def(init<optional<double,double>>(args("albedo","turbidity"),"a new object with specified parameters"))
            .def_readwrite("albedo",&parameter::albedo,"typical value 0.2")
            .def_readwrite("turbidity", &parameter::turbidity,"typical value 1.0 for clear air, 0.5 for dusty")
            ;

        class_<response>("RadiationDingmanResponse")
            .def_readwrite("rcs_radiation",&response::rcs_radiation)
            .def_readwrite("slope_factor",&response::slope_factor)
            ;
        typedef calculator<parameter> RadiationCalculator;
        class_<RadiationCalculator>("RadiationDingmanCalculator",
                "Radiation, (ref.: Dingman  )\n"
                "provide a function to calculate clear-sky radiation for inclined surface\n"
                "[W/m^2] units.\n",no_init
            )
            .def(init<const parameter&>(args("param"),"create a calculator using supplied parameter"))
            .def("rcs_radiation",&RadiationCalculator::compute_ra_radiation,args("latitude","rhumidity","temperature","utctime","slope","aspect","elevation"),
                "Calculate clear sky radiation, given specified parameters\n"
                "\n"
                "	 * param latitude [deg]\n"
                "	 * param utctime\n"
                "	 * param surface_normal [x,y,z]\n"
                "	 * param temperature in [degC]\n"
                "	 * param rhumidity in interval [0,100]\n"
                "	 * param elevation, [m]\n"
                "	 * return rcs_radiation in [W/m^2] units\n"
                "	 *\n"
                 )
            ;
    }
}
