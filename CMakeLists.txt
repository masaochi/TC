cmake_minimum_required(VERSION 3.1)

# project
enable_language(Fortran)
set(ProjectName tc++)
project(${ProjectName} CXX Fortran)

# compile option
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED OFF) # C++11 is sufficient but C++17 will accelerate calculation
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_CXX_EXTENSIONS OFF)

# external libraries
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_DIRS} ${MPI_Fortran_INCLUDE_DIRS})

if (BOOST_INCLUDE AND BOOST_LIB)
    include_directories(${BOOST_INCLUDE})
    set(BOOST_LIBRARIES ${BOOST_LIB})
else()
    find_package(Boost)
    if (Boost_FOUND)
        include_directories(${Boost_INCLUDE_DIRS})
    else()
        message(FATAL_ERROR "Boost was not found.
        Please specify the location of Boost include and library directories
        via -DBOOST_INCLUDE and -DBOOST_LIB options.")
    endif()
endif()

if (EIGEN3_INCLUDE)
    include_directories(${EIGEN3_INCLUDE})
else()
    find_package(Eigen3)
    if (Eigen3_FOUND)
        include_directories(${EIGEN3_INCLUDE_DIR})
    else()
        message(FATAL_ERROR "Eigen3 was not found.
        Please specify the location of EIGEN3 include directories
        via -DEIGEN3_INCLUDE option.")
    endif()
endif()

if (FFTW_INCLUDE AND FFTW_LIB)
    include_directories(${FFTW_INCLUDE})
    set(FFTW_LIBRARIES ${FFTW_LIB})
else()
    find_package(FFTW)
    if (FFTW_FOUND)
        include_directories(${FFTW_INCLUDE_DIRS})
    else()
        message(FATAL_ERROR "FFTW was not found.
        Please specify the location of FFTW include and library directories
        via -DFFTW_INCLUDE and -DFFTW_LIB options.")
    endif()
endif()

include_directories("src/include")

add_executable(${ProjectName}
  src/my_clock.cpp
  src/error_messages.cpp
  src/parallelization.cpp
  src/method.cpp
  src/plane_wave_basis.cpp
  src/symmetry.cpp
  src/crystal_structure.cpp
  src/potentials_pseudopot.cpp
  src/potentials_coulomb.cpp
  src/kpoints.cpp
  src/file_names.cpp
  src/jastrow_initialize.cpp
  src/jastrow_functions.cpp
  src/io_tc_files_read_input_in.cpp
  src/io_tc_files_read_eigen.cpp
  src/io_tc_files_read_scfinfo.cpp
  src/io_tc_files_dump_eigen.cpp 
  src/io_tc_files_dump_scfinfo.cpp
  src/io_tc_files_dump_bandplot.cpp
  src/read_qe.cpp
  src/read_qe_xml.cpp 
  src/read_upf.cpp
  src/bloch_states_initialize.cpp
  src/bloch_states_scfloop_filling_density.cpp
  src/bloch_states_scfloop_phik.cpp
  src/bloch_states_scfloop_eigenvalues.cpp
  src/total_energy.cpp
  src/spin.cpp
  src/calc_hamiltonian_all.cpp
  src/calc_hamiltonian_kinetic.cpp
  src/calc_hamiltonian_pseudo.cpp
  src/calc_hamiltonian_hf2h.cpp
  src/calc_hamiltonian_hf2x.cpp
  src/calc_hamiltonian_tc2h.cpp
  src/calc_hamiltonian_tc2x.cpp
  src/calc_hamiltonian_tc3a1.cpp
  src/calc_hamiltonian_tc3b1.cpp
  src/calc_hamiltonian_tc3a2a4b2b5.cpp
  src/calc_hamiltonian_tc3a3b3b4b6.cpp
  src/diagonalization.cpp
  src/diagonalization_davidson.cpp
  src/diagonalization_utils.cpp
  src/main.cpp
  src/read_qe_wfc.f90
)
target_link_libraries(${ProjectName} ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})
target_link_libraries(${ProjectName} ${FFTW_LIBRARIES}/libfftw3.a)

install(TARGETS ${ProjectName}
	DESTINATION bin)
