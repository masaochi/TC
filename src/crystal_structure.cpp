// [class Crystal structure]
// Lattice, atomic positions, atomic species, etc.

#include "include/header.hpp"

void CrystalStructure::set(const Eigen::Matrix3d &lattice_vectors,
                           const Eigen::Matrix3d &reciprocal_vectors,
                           const int &num_atomic_species,
                           const std::vector<int> &index_of_atoms,
                           const std::vector<Eigen::Vector3d> &atomic_position_cartesian,
                           std::ostream *ost)
{
    lattice_vectors_ = lattice_vectors;
    reciprocal_vectors_ = reciprocal_vectors;
    unit_cell_volume_ = lattice_vectors_.determinant();

    num_atomic_species_ = num_atomic_species;

    num_atoms_ = index_of_atoms.size();
    assert(atomic_position_cartesian.size() == num_atoms_);

    index_of_atoms_ = index_of_atoms;
    atomic_position_cartesian_ = atomic_position_cartesian;

    set_atomic_position_crystal();

    print_atomic_positions(ost);
}


// set atomic_number_[index of atomic species] and Z_valence_[index of atomic species].
//  e.g., for bulk Si with num_atoms_ = 2 and num_atomic_species_ = 1:
//    atomic_number_[0] = 14 and Z_valence_[0] = 4 (for [Ne]-core pseudopot.)
void CrystalStructure::set_upf(const std::vector<int> &atomic_number,
                               const std::vector<double> &Z_valence)
{
    assert(atomic_number.size() == num_atomic_species_);
    atomic_number_ = atomic_number;

    assert(Z_valence.size() == num_atomic_species_);
    Z_valence_ = Z_valence;
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
        atomic_position_crystal_.resize(num_atoms_);
    }
    MPI_Bcast(&index_of_atoms_[0], num_atoms_, MPI_INT, 0, MPI_COMM_WORLD);
    for (int iatom=0; iatom<num_atoms_; iatom++)
    {
        MPI_Bcast(atomic_position_cartesian_[iatom].data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
        MPI_Bcast(atomic_position_crystal_[iatom].data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
    }
}

void CrystalStructure::bcast_upf(const bool am_i_mpi_rank0)
{
    if (!am_i_mpi_rank0)
    {
        atomic_number_.resize(num_atomic_species_);
        Z_valence_.resize(num_atomic_species_);
    }
    else // rank0
    {
        assert(atomic_number_.size() == num_atomic_species_);
        assert(Z_valence_.size() == num_atomic_species_);
    }
    MPI_Bcast(&atomic_number_[0], num_atomic_species_, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Z_valence_[0], num_atomic_species_, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void CrystalStructure::set_atomic_position_crystal()
{
    atomic_position_crystal_.resize(num_atoms_);
   for (int iatom=0; iatom<num_atoms_; iatom++)
    {
        atomic_position_crystal_[iatom] 
            = lattice_vectors_.transpose().inverse() * atomic_position_cartesian_[iatom];
    }
}

void CrystalStructure::update_atomic_positions(const Eigen::VectorXd &position,
                                               const bool am_i_mpi_rank0,
                                               std::ostream *ost)
{
    const int ndegrees = position.size();
    assert(ndegrees%3==0 && ndegrees/3==num_atoms_);

    // save old data
    atomic_position_cartesian_old_ = atomic_position_cartesian_;

    // update coordinates
    for (int idegree=0; idegree<ndegrees; idegree++)
    {
        int iatom = idegree/3; // 0, 0, 0, 1, 1, 1, ...
        int idim = idegree%3; // 0, 1, 2, 0, 1, 2, ...

        atomic_position_cartesian_[iatom](idim) = position(idegree);
    }

    // update atomic_position_crystal
    set_atomic_position_crystal();

    if (am_i_mpi_rank0) { *ost << "   Print NEW atomic coordinates." << std::endl; }
    if (am_i_mpi_rank0) { print_atomic_positions(ost); }
}

void CrystalStructure::print_atomic_positions(std::ostream *ost) const
{
    *ost << std::endl;
    *ost << "   Atomic coordinates (in Bohr, cartesian)" << std::endl;
    for (int iatom=0; iatom<num_atoms_; iatom++)
    {
        *ost << "    Atom " << iatom 
             << " ( " << atomic_position_cartesian_[iatom](0)
             << " , " << atomic_position_cartesian_[iatom](1)
             << " , " << atomic_position_cartesian_[iatom](2) << " )" << std::endl;;
    }
    *ost << std::endl;
    
    *ost << "   Atomic coordinates (crystal coordinate)" << std::endl;
    for (int iatom=0; iatom<num_atoms_; iatom++)
    {
        *ost << "    Atom " << iatom 
             << " ( " << atomic_position_crystal_[iatom](0)
             << " , " << atomic_position_crystal_[iatom](1)
             << " , " << atomic_position_crystal_[iatom](2) << " )" << std::endl;;
    }
    *ost << std::endl;
}

bool CrystalStructure::do_lattice_vectors_coincide(const Eigen::Matrix3d &lattice_vectors) const
{
    double diff = 0.0;
    for (int idim=0; idim<3; idim++)
    {
        for (int jdim=0; jdim<3; jdim++)
        {
            diff += std::abs(lattice_vectors(idim, jdim) - lattice_vectors_(idim, jdim));
        }
    }
    return (diff < 1e-8);
}

void CrystalStructure::update_atomic_positions_by_file(const std::vector<Eigen::Vector3d> &atomic_position_cartesian,
                                                       const bool am_i_mpi_rank0,
                                                       std::ostream *ost)
{
    assert(atomic_position_cartesian.size()==num_atoms_);
    atomic_position_cartesian_ = atomic_position_cartesian;

    // update atomic_position_crystal
    set_atomic_position_crystal();

    if (am_i_mpi_rank0) { *ost << "   Print NEW atomic coordinates read from file." << std::endl; }
    if (am_i_mpi_rank0) { print_atomic_positions(ost); }
}
