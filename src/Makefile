## Intel compiler
CXX = mpiicpc
F90 = mpiifort
COMP_LINK = ${F90}
OPT_CXX = -O2 -std=c++17 -qopenmp #-std=c++11
OPT_F90 = -O2
OPT_LINK = -O2 -qopenmp -nofor-main -cxxlib
LIBS_LINK = 
OBJS_F90 =
SRC_F90 = io_qe_files_read_wfc.f90

## gcc + openmpi compiler
#CXX = mpicxx
#F90 = mpifort
#COMP_LINK = ${CXX}
#OPT_CXX = -O2 -fopenmp
#OPT_F90 = -O2
#OPT_LINK = -O2 -fopenmp
#LIBS_LINK = /usr/local/openmpi/lib/libmpi_mpifh.so -lquadmath /usr/local/gcc/10.1.0/lib64/libgfortran.a
#OBJS_F90 = io_qe_files_read_wfc.o
#SRC_F90 = 

# FFTW
FFTW_DIR = /home/user/fftw
LIBS_FFTW = ${FFTW_DIR}/lib/libfftw3.a
INCS_FFTW = -I${FFTW_DIR}/include

# Boost
BOOST_DIR = /home/user/boost
LIBS_BOOST = -L${BOOST_DIR}/lib
INCS_BOOST = -I${BOOST_DIR}/include

# Eigen
EIGEN_DIR = /home/user/eigen
INCS_EIGEN = -I${EIGEN_DIR}

LIBS = ${LIBS_FFTW} ${LIBS_BOOST} ${LIBS_LINK}
INCS = ${INCS_FFTW} ${INCS_BOOST} ${INCS_EIGEN}

OBJS_CXX = my_clock.o error_messages.o atomic_species.o parallelization.o \
method.o plane_wave_basis.o symmetry.o crystal_structure.o \
potentials_pseudopot.o potentials_coulomb.o potentials_derivative.o kpoints.o file_names.o \
jastrow_initialize.o jastrow_initialize_polynomial.o jastrow_functions.o \
io_tc_files_read_input_in.o io_tc_files_read_eigen.o io_tc_files_read_scfinfo.o \
io_tc_files_read_crystal_structure.o \
io_tc_files_dump_eigen.o io_tc_files_dump_scfinfo.o io_tc_files_dump_bandplot.o \
io_tc_files_dump_crystal_structure.o \
io_qmc_files_dump_jastrow_plt.o io_qmc_files_dump_pwfn.o \
io_qmc_files_read_casl.o io_qmc_files_dump_casl.o \
io_qe_files_read.o io_qe_files_read_xml.o io_qe_files_read_upf.o \
bloch_states_initialize.o bloch_states_scfloop_filling_density.o \
bloch_states_scfloop_phik.o bloch_states_scfloop_eigenvalues.o \
total_energy.o total_energy_ewald.o total_energy_force.o \
spin.o calc_hamiltonian_all.o \
calc_hamiltonian_kinetic.o calc_hamiltonian_pseudo.o \
calc_hamiltonian_hf2h.o calc_hamiltonian_hf2x.o \
calc_hamiltonian_tc2h.o calc_hamiltonian_tc2x.o \
calc_hamiltonian_tc3a1.o calc_hamiltonian_tc3b1.o \
calc_hamiltonian_tc3a2a4b2b5.o \
calc_hamiltonian_tc3a3b3b4b6.o \
calc_hamiltonian_force.o \
diagonalization.o diagonalization_davidson.o diagonalization_utils.o \
diagonalization_structural_optimization.o \
main.o

OBJS = ${OBJS_CXX} ${OBJS_F90}

PROGRAM = tc++

all: $(PROGRAM)
$(PROGRAM): ${OBJS} ${SRC_F90}
	${COMP_LINK} ${OPT_LINK} ${OBJS} ${LIBS} ${SRC_F90} -o tc++

.SUFFIXES: .o .cpp .f90
.cpp.o:
	${CXX} ${OPT_CXX} ${INCS} $< -c
.f90.o:
	${F90} ${OPT_F90} $< -c

.PHONY: clean
clean:
	rm -f *~ include/*~ *.o $(PROGRAM)
