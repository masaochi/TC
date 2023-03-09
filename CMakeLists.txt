cmake_minimum_required(VERSION 3.13)

project(Tcpp
		VERSION 1.3
		LANGUAGES CXX Fortran)

## Project options
option(TCPP_WITH_MPI "TCpp: Build with MPI support" ON)
option(TCPP_WITH_Boost "TCpp: Build with Boost support" ON)
option(TCPP_WITH_Eigen3 "TCpp: Build with Eigen3 support" ON)
option(TCPP_WITH_FFTW3 "TCpp: Build with FFTW3 support" ON)
if (FFTW_INCLUDE AND FFTW_LIB)
   	set(TCPP_WITH_FFTW3 ON)
endif ()

## Set standard options
if (NOT DEFINED CMAKE_BUILD_TYPE)
   	set(CMAKE_BUILD_TYPE Release)
endif ()

# Check C++ standards
if ("cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
   	set(CMAKE_CXX_STANDARD 17)
else ()
	message(WARNING
			"Compiler ${CMAKE_CXX_COMPILER} does not support C++17 standard.\n"
			"We recommend upgrading the compiler to get better support and performance.")
	set(CMAKE_CXX_STANDARD 11)
endif ()
set(CMAKE_CXX_EXTENSIONS OFF)
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(GNUInstallDirs)

## External libraries
if (TCPP_WITH_MPI)
	find_package(MPI REQUIRED)
else ()
	message(FATAL_ERROR "Currently TC++ requires MPI")
endif ()

if (TCPP_WITH_Boost)
	if (BOOST_INCLUDE AND BOOST_LIB)
	   	add_library(Boost::boost INTERFACE IMPORTED)
    		set_target_properties(Boost::boost PROPERTIES
      						   IMPORTED_LOCATION ${BOOST_LIB}
						   INTERFACE_INCLUDE_DIRECTORIES ${BOOST_INCLUDE})
	else ()
		find_package(Boost)
      		if (NOT Boost_FOUND)
      	 	   	message(FATAL_ERROR "Boost was not found.\n"
         		   		    "Please specify the location of Boost include and library directories\n"
         		     		    "via -DBOOST_INCLUDE and -DBOOST_LIB options.")
      		endif ()
   	endif ()
else ()
	message(FATAL_ERROR "Currently TC++ requires Boost")
endif ()

if (TCPP_WITH_Eigen3)
	if (EIGEN3_INCLUDE)
		add_library(Eigen3::Eigen INTERFACE IMPORTED)
		target_include_directories(Eigen3::Eigen INTERFACE ${EIGEN3_INCLUDE})
	else ()
    		find_package(Eigen3)
		if (NOT Eigen3_FOUND)
        	   	message(FATAL_ERROR "Eigen3 was not found.\n"
        				    "Please specify the location of EIGEN3 include directories\n"
        				    "via -DEIGEN3_INCLUDE option.")
    		endif ()

		# for old versions of Eigen3 (before Eigen 3.3)
		if ((NOT TARGET Eigen3::Eigen) AND (DEFINED EIGEN3_INCLUDE_DIR))
       		   	 add_library(AliasEigen3 INTERFACE)
       			 target_include_directories(AliasEigen3 INTERFACE ${EIGEN3_INCLUDE_DIRS})
       			 add_library(Eigen3::Eigen ALIAS AliasEigen3)
		endif ()
    	endif ()
else ()
	message(FATAL_ERROR "Currently TC++ requires Eigen3")
endif ()

if (TCPP_WITH_FFTW3)
	if (FFTW_INCLUDE AND FFTW_LIB)
		add_library(FFTW_Target INTERFACE IMPORTED)
		set_target_properties(FFTW_Target PROPERTIES
						  INTERFACE_INCLUDE_DIRECTORIES ${FFTW_INCLUDE}
						  INTERFACE_LINK_LIBRARIES "${FFTW_LIB}/libfftw3.a")
	else ()
		find_package(FFTW)
    		if (FFTW_FOUND)
			set(FFTW_Target "FFTW::Double")
    		else ()
			message("FFTW3 was not found. Try to find INTEL MKL...")
			set(TCPP_WITH_FFTW3 OFF)
		endif ()
    	endif ()
endif ()

if (NOT TCPP_WITH_FFTW3)
	if (CMAKE_CXX_COMPILER_ID MATCHES "INTEL.*")
	   	find_package(MKL REQUIRED)
		set(FFTW_Target "MKL::MKL")
	else ()
		message(FATAL_ERROR "Currently TC++ requires FFTW3 or Intel MKL\n"
				    "To specify the location of FFTW3 manually, please use\n"
				    "-DFFTW_INCLUDE and -DFFTW_LIB options.")
	endif ()
endif ()

# Main binary target
add_executable(TCpp)
add_executable(TCpp::tc++ ALIAS TCpp)
set_target_properties(TCpp PROPERTIES
				     EXPORT_NAME tc++
				     OUTPUT_NAME tc++)
add_subdirectory(src)
