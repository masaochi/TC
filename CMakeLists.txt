cmake_minimum_required(VERSION 3.13)

project(TCpp
		VERSION 1.3
		LANGUAGES CXX Fortran)

## Project options
option(TCPP_WITH_MPI "Tcpp: Build with MPI support" ON)
option(TCPP_WITH_Boost "Tcpp: Build with Boost support" ON)
option(TCPP_WITH_Eigen3 "Tcpp: Build with Eigen3 support" ON)
option(TCPP_WITH_FFTW3 "Tcpp: Build with FFTW3 support" OFF)

## Set standard options
if (NOT DEFINED CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif ()
# Check C++ standards
if ("cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
	set(CMAKE_CXX_STANDARD 17)
else ()
	# TODO: Mark compiler standard for deprecation and announce release
	message(WARNING
			"Compiler ${CMAKE_CXX_COMPILER} does not support C++17 standard.\n"
			"We highly recommend upgrading the compiler to get better support and performance.")
	set(CMAKE_CXX_STANDARD 11)
endif ()
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(GNUInstallDirs)

## External libraries
if (TCPP_WITH_MPI)
	find_package(MPI REQUIRED)
else ()
	# TODO: Allow for non-mpi builds
	message(FATAL_ERROR "Currently Tcpp requires MPI")
endif ()

if (TCPP_WITH_Boost)
	# TODO: What components?
	find_package(Boost REQUIRED)
else ()
	# TODO: Allow for non-boost builds
	message(FATAL_ERROR "Currently Tcpp requires Boost")
endif ()
if (TCPP_WITH_Eigen3)
	find_package(Eigen3 REQUIRED)
else ()
	# TODO: Allow for non-boost builds
	message(FATAL_ERROR "Currently Tcpp requires Eigen3")
endif ()
if (TCPP_WITH_FFTW3)
	find_package(FFTW REQUIRED)
	set(FFTW_Target "FFTW::Double")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel.*")
	find_package(MKL REQUIRED)
	set(FFTW_Target "MKL::MKL")
else ()
	# TODO: Allow for non-boost builds
	message(FATAL_ERROR "Currently Tcpp requires FFTW3 or Intel MKL")
endif ()

# Main binary target
add_executable(TCpp_tc++)
# Target is namespaced for common practice. Add appropriate namespaced alias, export name and fix the executable name
add_executable(TCpp::tc++ ALIAS TCpp_tc++)
set_target_properties(TCpp_tc++ PROPERTIES
		EXPORT_NAME tc++
		OUTPUT_NAME tc++)
add_subdirectory(src)
