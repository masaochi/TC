// [class Spin]
// spin setting (no spin, noncollinear etc.), magnetization

#include "include/header.hpp"

void Spin::set(const std::string &noncolin, const std::string &lsda)
{
    if (noncolin == "true") // QE keyword
    {
        // nspin = 4 in QE (non-collinear)
        is_spinor_ = true;
        num_independent_spins_ = 1;
    }
    else
    {
        if (lsda == "true") // QE keyword
        {
            // nspin = 2 in QE (collinaer, spin-polarized)
            is_spinor_ = false;
            num_independent_spins_ = 2;
        }
        else
        {
            // nspin = 1 in QE (collinear, spin-non-polarized)
            is_spinor_ = false;
            num_independent_spins_ = 1;
        }
    }
}

void Spin::bcast(const bool am_i_mpi_rank0)
{
    MPI_Bcast(&is_spinor_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_independent_spins_, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
