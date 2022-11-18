// [class Method]
// method (HF, TC, BITC, FREE) and mode (SCF, BAND)

#include "include/header.hpp"

namespace
{
    
void bcast_each_keyword(std::string &s, bool am_i_mpi_rank0)
{
    int length = s.size();
    MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!am_i_mpi_rank0) { s.resize(length); }
    MPI_Bcast(const_cast<char*>(s.data()), length, MPI_CHAR, 0, MPI_COMM_WORLD);
}
    
}

void Method::set_calc_method(const std::string &calc_method)
{
    const std::vector<std::string> list_method{"HF", "TC", "BITC", "FREE"};
    error_messages::compare_list("method", calc_method, list_method); // If not included in list_method, stop the program
    calc_method_ = calc_method;
}

void Method::set_calc_mode(const std::string &calc_mode)
{
    const std::vector<std::string> list_mode{"SCF", "BAND"};
    error_messages::compare_list("mode", calc_mode, list_mode);
    calc_mode_ = calc_mode;
}

void Method::check_consistency_of_calc_mode(const std::string &calc_mode_QE) const
{
    if (!((calc_mode_QE=="scf" && calc_mode_=="SCF") || 
          (calc_mode_QE=="bands" && calc_mode_=="BAND")))
    {
        error_messages::stop("calc_mode (input.in) and calculation in QE ('scf' or 'bands') are not consistent.");
    }
}

void Method::bcast(const bool am_i_mpi_rank0)
{
    std::string s1, s2;

    if (am_i_mpi_rank0) { s1 = calc_method_; }
    bcast_each_keyword(s1, am_i_mpi_rank0);
    if (!am_i_mpi_rank0) { set_calc_method(s1); }
    
    if (am_i_mpi_rank0) { s2 = calc_mode_; }
    bcast_each_keyword(s2, am_i_mpi_rank0);
    if (!am_i_mpi_rank0) { set_calc_mode(s2); }
}
