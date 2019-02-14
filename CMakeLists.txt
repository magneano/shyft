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

set(python_min_required_version 3.6)
set(boost_min_required_version 1.67)

# Get the full version for Shyft
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/VERSION SHYFT_VERSION_STRING)
message("Configuring for Shyft version: " ${SHYFT_VERSION_STRING})

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

# The dependencies directory if set
if(DEFINED ENV{SHYFT_DEPENDENCIES_DIR})
    set(SHYFT_DEPENDENCIES_DIR $ENV{SHYFT_DEPENDENCIES_DIR})
    # Our code requires an absolute directory for the dependencies
    get_filename_component(SHYFT_DEPENDENCIES ${SHYFT_DEPENDENCIES_DIR} ABSOLUTE)

    # Prevent conflict wit system boost if present
    set(Boost_NO_SYSTEM_PATHS ON)

    # Armadillo, dlib and doctest provide their own cmake config files
    # Setting CMAKE_PREFIX_PATH so that `find_module` will find them
    set(CMAKE_PREFIX_PATH "${SHYFT_DEPENDENCIES}" ${CMAKE_PREFIX_PATH})
    message("Using SHYFT_DEPENDENCIES directory: " ${SHYFT_DEPENDENCIES})
    message("You can change the above via the SHYFT_DEPENDENCIES_DIR environment variable.")
else()
    message("Environment variable SHYFT_DEPENDENCIES_DIR not set.")
endif()

# shyft uses boost, dlib, armadillo and doctest, and python +  numpy

# Find python version (Want matching Boost::python version)
find_package(PythonInterp ${python_min_required_version} REQUIRED)
find_package(PythonLibs ${python_min_required_version} REQUIRED)
set(boost_py "python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")

# Find required packages
find_package(doctest REQUIRED)
find_package(dlib REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Boost ${boost_min_required_version} COMPONENTS ${boost_py} serialization filesystem system REQUIRED)

# The compiler requirement,flags definitions and include that applies for all
# Notice that we need MINIMIZE_SIZE to ensure g++ vs. ms c++ compat
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
add_compile_options(-fPIC -fexceptions -pthread -Wno-deprecated-declarations -Wno-undef)

# Armadillo options
add_definitions("-DARMA_USE_CXX11 -DARMA_NO_DEBUG")

# Boost options
add_definitions("-DBOOST_ERROR_CODE_HEADER_ONLY -DBOOST_VARIANT_MINIMIZE_SIZE -DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS -DBOOST_MPL_LIMIT_LIST_SIZE=30 -DBOOST_MPL_LIMIT_VECTOR_SIZE=30")

# Include
include_directories(${CMAKE_SOURCE_DIR} ${ARMADILLO_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})

# Convenience function for adding python modules
function (add_module target_name)
    set(options SOURCE LIB)
    cmake_parse_arguments(ARG "" "" "${options}" ${ARGN})
    add_library(${target_name} MODULE ${ARG_SOURCE})
    target_link_libraries(${target_name} PRIVATE ${ARG_LIB})
    set_target_properties(${target_name} PROPERTIES PREFIX "_")
    set_target_properties(${target_name} PROPERTIES CXX_VISIBILITY_PRESET "hidden")
    set_property(TARGET ${target_name} APPEND PROPERTY COMPILE_DEFINITIONS SHYFT_EXTENSION)
endfunction()

# C++ core and tests
add_subdirectory(core)
add_subdirectory(api)

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

# Python extensions
if(BUILD_PYTHON_EXTENSIONS)
    add_subdirectory(api/boostpython)
endif()
