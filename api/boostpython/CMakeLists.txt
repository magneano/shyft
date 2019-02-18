# Shyft python extension modules

set(core_api  "api")
set(core_api_sources
    api_actual_evapotranspiration.cpp  api_kirchner.cpp
    api_cell_environment.cpp           api_precipitation_correction.cpp
    api_dtss.cpp                       api_priestley_taylor.cpp
    api_gamma_snow.cpp                 api_pt_gs_k.cpp
    api_geo_cell_data.cpp              api_region_environment.cpp
    api_geo_point.cpp                  api_routing.cpp
    api_glacier_melt.cpp               api.cpp
    api_skaugen.cpp                    api_hbv_actual_evapotranspiration.cpp
    api_state.cpp                      api_kalman.cpp
    api_hbv_physical_snow.cpp          api_target_specification.cpp
    api_hbv_snow.cpp                   api_time_axis.cpp
    api_hbv_soil.cpp                   api_time_series.cpp
    api_hbv_tank.cpp                   api_utctime.cpp
    api_interpolation.cpp              api_vectors.cpp
)

# Force all Boost packages as direct dependencies. This avoids some linking
# issues with secondary dependencies (i.e. Boost::system)
set(api_link_libraries shyft_core shyft_api ${ARMADILLO_LIBRARIES} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES} dlib::dlib)

# Core
add_module(${core_api} ${core_api_sources})
target_link_libraries(${core_api} PRIVATE ${api_link_libraries})

# define each shyft method-stack that exposes complete stacks
foreach(shyft_stack  "pt_gs_k" "pt_ss_k" "pt_hs_k" "hbv_stack" "pt_hps_k")
    add_module(${shyft_stack} "${shyft_stack}.cpp")
    target_link_libraries(${shyft_stack} PRIVATE ${api_link_libraries})
endforeach()
