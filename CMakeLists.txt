# CMake build system for Shyft
# ============================
#
# This requires shyft_dependencies installed as well as a C++14/C++17 compliant compiler.
#
# Available options:
#
#   BUILD_TESTING: default ON
#       build test programs and generates the "test" target
#   BUILD_PYTHON_EXTENSIONS: default ON
#       build Python extensions for Shyft
#
# The next environment variables are honored:
#
#   SHYFT_DEPENDENCIES_DIR: default ${PROJECT_SOURCE_DIR}/../shyft_dependencies
#
cmake_minimum_required(VERSION 3.9.0)
project(shyft C CXX)

# Get the full version for Shyft
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/VERSION SHYFT_VERSION_STRING)

message("Configuring for Shyft version: " ${SHYFT_VERSION_STRING})

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

# options
option(BUILD_TESTING "Build test programs for SHYFT C++ core library" ON)
option(BUILD_PYTHON_EXTENSIONS "Build Python extensions for SHYFT" ON)
set(SHYFT_DEFAULT_BUILD_TYPE "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "No build type specified. Defaulting to '${SHYFT_DEFAULT_BUILD_TYPE}'.")
    set(CMAKE_BUILD_TYPE ${SHYFT_DEFAULT_BUILD_TYPE} CACHE STRING
        "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# The dependencies directory
if(DEFINED ENV{SHYFT_DEPENDENCIES_DIR})
  set(SHYFT_DEPENDENCIES_DIR $ENV{SHYFT_DEPENDENCIES_DIR})
else()
  set(SHYFT_DEPENDENCIES_DIR "${PROJECT_SOURCE_DIR}/../shyft_dependencies")
endif()
# Our code requires an absolute directory for the dependencies
get_filename_component(SHYFT_DEPENDENCIES ${SHYFT_DEPENDENCIES_DIR} ABSOLUTE)
message("SHYFT_DEPENDENCIES directory: " ${SHYFT_DEPENDENCIES})
message("You can override the above via the $SHYFT_DEPENDENCIES_DIR environment variable.")

# shyft uses boost, dlib, armadillo and doctest, and python +  numpy
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "${SHYFT_DEPENDENCIES}")

find_package(doctest REQUIRED)
find_package(dlib REQUIRED)
find_package(Armadillo REQUIRED)
#
get_filename_component(arma_libs  ${ARMADILLO_LIBRARIES} REALPATH CACHE)
if(${CMAKE_VERSION} VERSION_LESS "3.10.0") 
    # unfortunately, .travis does not have updated cmake, so find_package won't work as we need.
    message("old cmake, using workaround for boost libraries")
    set(boost_system_lib ${SHYFT_DEPENDENCIES}/lib/libboost_system.so)
    set(boost_filesystem_lib ${SHYFT_DEPENDENCIES}/lib/libboost_filesystem.so)
    set(boost_python3_lib ${SHYFT_DEPENDENCIES}/lib/libboost_python36.so)
    set(boost_serialization_lib ${SHYFT_DEPENDENCIES}/lib/libboost_serialization.so)
else()
    find_package(Boost COMPONENTS python3 serialization filesystem system REQUIRED)
    get_property(boost_system_lib TARGET Boost::system PROPERTY IMPORTED_LOCATION_RELEASE)
    get_property(boost_filesystem_lib TARGET Boost::filesystem PROPERTY IMPORTED_LOCATION_RELEASE)
    get_property(boost_serialization_lib  TARGET Boost::serialization PROPERTY IMPORTED_LOCATION_RELEASE)
    get_property(boost_python3_lib TARGET Boost::python3 PROPERTY IMPORTED_LOCATION_RELEASE)

endif()

set(boost_link_libraries ${boost_system_lib} ${boost_filesystem_lib} ${boost_serialization_lib})
set(boost_py_link_libraries ${boost_link_libraries} ${boost_python3_lib})

# The compiler requirement,flags definitions and include that applies for all
# Notice that we need MINIMIZE_SIZE to ensure g++ vs. ms c++ compat
set(CMAKE_CXX_STANDARD 17)
add_compile_options(-fPIC -fexceptions -pthread -Wno-deprecated-declarations)
add_definitions("-DARMA_USE_CXX11 -DARMA_NO_DEBUG -DBOOST_VARIANT_MINIMIZE_SIZE -DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS -DBOOST_MPL_LIMIT_LIST_SIZE=30 -DBOOST_MPL_LIMIT_VECTOR_SIZE=30")
include_directories(${CMAKE_SOURCE_DIR} ${SHYFT_DEPENDENCIES}/include) # we could add more explicit list boost,dlib,armadillo,doctest here


# C++ core and tests
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(core)
  add_subdirectory(api)
  add_subdirectory(test)
endif(BUILD_TESTING)

# Python extensions
if(BUILD_PYTHON_EXTENSIONS)
  add_subdirectory(api/boostpython)
endif(BUILD_PYTHON_EXTENSIONS)
