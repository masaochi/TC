// [class Symmetry]
// symmetry operations

#include "include/header.hpp"

void Symmetry::set(const Spin &spin, const int num_symmetry_operations,
                   const std::vector<Eigen::Matrix3i> &rotation,
                   const std::vector<Eigen::Vector3d> &translation,
                   const std::string &no_t_rev, const std::string &noinv)
{
    assert(num_symmetry_operations >= 0);
    assert(rotation.size() == num_symmetry_operations);
    assert(translation.size() == num_symmetry_operations);

    num_symmetry_operations_ = num_symmetry_operations;
    rotation_ = rotation;
    translation_ = translation;

    // use time-reversal symmetry to make symmetrically-equivalent k-points
    if (!(no_t_rev=="true" && noinv=="true") && (spin.is_spinor() || spin.num_independent_spins()==2) )
    {
        error_messages::stop("Please specify noinv = .true. && no_t_rev = .true. in QE unless nspin=1...");
    }
    is_time_reversal_symmetric_ = (noinv=="true") ? false : true; // Note!
}

void Symmetry::bcast(const bool am_i_mpi_rank0)
{
    MPI_Bcast(&num_symmetry_operations_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!am_i_mpi_rank0)
    {
        rotation_.resize(num_symmetry_operations_);
        translation_.resize(num_symmetry_operations_);
    }
    for (int isym=0; isym<num_symmetry_operations_; isym++)
    {
        MPI_Bcast(rotation_[isym].data(), 9, MPI_INT, 0, MPI_COMM_WORLD);        
        MPI_Bcast(translation_[isym].data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    MPI_Bcast(&is_time_reversal_symmetric_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
}
