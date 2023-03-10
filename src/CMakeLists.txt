## Configure main executable
target_sources(TCpp_tc++ PRIVATE
  			 my_clock.cpp
			 error_messages.cpp
			 parallelization.cpp
			 method.cpp
			 plane_wave_basis.cpp
			 symmetry.cpp
			 crystal_structure.cpp
			 potentials_pseudopot.cpp
			 potentials_coulomb.cpp
			 kpoints.cpp
			 file_names.cpp
			 jastrow_initialize.cpp
			 jastrow_functions.cpp
			 io_tc_files_read_input_in.cpp
			 io_tc_files_read_eigen.cpp
			 io_tc_files_read_scfinfo.cpp
			 io_tc_files_dump_eigen.cpp 
			 io_tc_files_dump_scfinfo.cpp
			 io_tc_files_dump_bandplot.cpp
			 read_qe.cpp
			 read_qe_xml.cpp
			 read_upf.cpp
			 bloch_states_initialize.cpp
			 bloch_states_scfloop_filling_density.cpp
			 bloch_states_scfloop_phik.cpp
			 bloch_states_scfloop_eigenvalues.cpp
			 total_energy.cpp
			 spin.cpp
			 calc_hamiltonian_all.cpp
			 calc_hamiltonian_kinetic.cpp
			 calc_hamiltonian_pseudo.cpp
			 calc_hamiltonian_hf2h.cpp
			 calc_hamiltonian_hf2x.cpp
			 calc_hamiltonian_tc2h.cpp
			 calc_hamiltonian_tc2x.cpp
			 calc_hamiltonian_tc3a1.cpp
			 calc_hamiltonian_tc3b1.cpp
			 calc_hamiltonian_tc3a2a4b2b5.cpp
			 calc_hamiltonian_tc3a3b3b4b6.cpp
			 diagonalization.cpp
			 diagonalization_davidson.cpp
			 diagonalization_utils.cpp
			 main.cpp
			 read_qe_wfc.f90
)
target_include_directories(TCpp_tc++ PRIVATE include)

## Link external libraries
target_link_libraries(TCpp_tc++ PRIVATE ${FFTW_Target})
if (TCPP_WITH_Boost)
	target_link_libraries(TCpp_tc++ PRIVATE Boost::boost)
endif ()
if (TCPP_WITH_MPI)
	target_link_libraries(TCpp_tc++ PRIVATE MPI::MPI_CXX MPI::MPI_Fortran)
endif ()
iF (TCPP_WITH_Eigen3)
	target_link_libraries(TCpp_tc++ PRIVATE Eigen3::Eigen)
endif ()

## Install
install(TARGETS TCpp_tc++
		DESTINATION ${CMAKE_INSTALL_BINDIR})
