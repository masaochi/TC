// [class Crystal structure]
// Lattice, atomic positions, atomic species, etc.

#include "include/header.hpp"

void CrystalStructure::set(const Eigen::Matrix3d &lattice_vectors,
                           const Eigen::Matrix3d &reciprocal_vectors,
                           const int &num_atomic_species,
                           const std::vector<int> &index_of_atoms,
                           const std::vector<Eigen::Vector3d> &atomic_position_cartesian)
{
    lattice_vectors_ = lattice_vectors;
    reciprocal_vectors_ = reciprocal_vectors;
    unit_cell_volume_ = lattice_vectors_.determinant();

    num_atomic_species_ = num_atomic_species;

    num_atoms_ = index_of_atoms.size();
    assert(atomic_position_cartesian.size() == num_atoms_);

    index_of_atoms_.resize(num_atoms_);
    index_of_atoms_ = index_of_atoms;
    atomic_position_cartesian_.resize(num_atoms_);
    atomic_position_cartesian_ = atomic_position_cartesian;
}

void CrystalStructure::bcast(const bool am_i_mpi_rank0)
{
    MPI_Bcast(lattice_vectors_.data(), 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(reciprocal_vectors_.data(), 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (!am_i_mpi_rank0) { unit_cell_volume_ = lattice_vectors_.determinant(); }

    MPI_Bcast(&num_atomic_species_, 1, MPI_INT, 0, MPI_COMM_WORLD);

    num_atoms_ = index_of_atoms_.size();
    MPI_Bcast(&num_atoms_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!am_i_mpi_rank0)
    {
        index_of_atoms_.resize(num_atoms_);
        atomic_position_cartesian_.resize(num_atoms_);
    }
    MPI_Bcast(&index_of_atoms_[0], num_atoms_, MPI_INT, 0, MPI_COMM_WORLD);
    for (int iatom=0; iatom<num_atoms_; iatom++)
    {
        MPI_Bcast(atomic_position_cartesian_[iatom].data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

