#ifndef TC_HEADER_HPP
#define TC_HEADER_HPP

#define NDEBUG
// Eigen macros should be defined before including Eigen headers.
#define EIGEN_NO_DEBUG

// all quantities are in atomic unit in TC++

#include <cstdlib>
#include <cmath>

#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
// Eigen macros should be defined before including Eigen headers.
#define EIGEN_DONT_PARALLELIZE
#endif

#include <chrono>
#include <complex>
#include <vector>
#include <algorithm>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "fftw3.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/algorithm/string.hpp>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using Complex = std::complex<double>;
constexpr double PI = 3.14159265358979323846; // std::numbers::pi in C++20
constexpr double FourPI = 4*PI;
//constexpr Complex I {0.0,1.0}; // constexpr Complex does not work for some (old) compilers...
const Complex I {0.0,1.0};
constexpr double Ht_in_eV = 27.21138505;
constexpr double Bohr_in_ang = 0.529177210903;

// TC++ header files
#include "my_clock.hpp"
#include "error_messages.hpp"
#include "atomic_species.hpp"
#include "spin.hpp"
#include "file_names.hpp"
#include "method.hpp"

#include "crystal_structure.hpp"
#include "symmetry.hpp"
#include "kpoints.hpp"

#include "plane_wave_basis.hpp"
#include "bloch_states.hpp"
#include "total_energy.hpp"

#include "parallelization.hpp"

#include "jastrow.hpp"
#include "potentials.hpp"
#include "calc_hamiltonian.hpp"
#include "diagonalization.hpp"

#include "io_qe_files.hpp"
#include "io_qe_files_upf.hpp"
#include "io_tc_files.hpp"
#include "io_qmc_files.hpp"

#endif // TC_HEADER_HPP
