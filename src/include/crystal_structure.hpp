// [class CrystalStructure]
// Lattice, atomic positions, etc.

// atomic unit (bohr etc.)

#ifndef TC_CRYSTAL_STRUCTURE_HPP
#define TC_CRYSTAL_STRUCTURE_HPP

class CrystalStructure
{
private:
    // (0,0:2) = a1, (1,0:2) = a2, (2,0:2) = a3 in cartesian coordinate
    // (in bohr. For example, a1 = (-5.13, 0.0, 5.13))
    Eigen::Matrix3d lattice_vectors_;
    // (0,0:2) = b1, (1,0:2) = b2, (2,0:2) = b3 in cartesian coordinate
    // (in bohr^{-1}. e.g., lattice_vectors_ * reciprocal_vectors_.transpose() = 2pi*identity
    Eigen::Matrix3d reciprocal_vectors_;

    double unit_cell_volume_; // = lattice_vectors_.determinant()

    // e.g. index_of_atoms_[0] = 0 (Ti), [1] = 1 (S), [2] = 1 (S) for TiS2
    //      atomic_number_[0] = 22, [1] = 16, [2] = 16 (see also atomic_species.hpp)
    // In this case, num_atoms_ = 3, num_atomic_species_ = 2
    int num_atoms_;
    int num_atomic_species_;
    std::vector<int> index_of_atoms_;
    std::vector<int> atomic_number_; // atomic_number_[index of atomic species (< num_atomic_species_)]
    std::vector<double> Z_valence_; // Z_valence_[index of atomic species (< num_atomic_species_)]
    std::vector<Eigen::Vector3d> atomic_position_cartesian_; // cartesian coordinate of each atom (in bohr)
    std::vector<Eigen::Vector3d> atomic_position_cartesian_old_; // for structural opt. Atomic positions in a previous loop (in bohr)
    std::vector<Eigen::Vector3d> atomic_position_crystal_; // crystal coordinate 

public:
    const Eigen::Matrix3d &lattice_vectors() const { return lattice_vectors_; }
    const Eigen::Matrix3d &reciprocal_vectors() const { return reciprocal_vectors_; }
    double unit_cell_volume() const { return unit_cell_volume_; }
    int num_atoms() const { return num_atoms_; }
    int num_atomic_species() const { return num_atomic_species_; }
    const std::vector<int> &index_of_atoms() const { return index_of_atoms_; }
    const std::vector<int> &atomic_number() const { return atomic_number_; }
    const std::vector<double> &Z_valence() const { return Z_valence_; }
    const std::vector<Eigen::Vector3d> &atomic_position_cartesian() const { return atomic_position_cartesian_; }
    const std::vector<Eigen::Vector3d> &atomic_position_cartesian_old() const { return atomic_position_cartesian_old_; }
    const std::vector<Eigen::Vector3d> &atomic_position_crystal() const { return atomic_position_crystal_; }

    CrystalStructure() : unit_cell_volume_(0.0), num_atoms_(0), num_atomic_species_(0) {}

    // set(): called from MPI rank=0
    void set(const Eigen::Matrix3d &lattice_vectors,
             const Eigen::Matrix3d &reciprocal_vectors,
             const int &num_atomic_species,
             const std::vector<int> &index_of_atoms,
             const std::vector<Eigen::Vector3d> &atomic_position_cartesian,
             std::ostream *ost);
    void set_upf(const std::vector<int> &atomic_number,
                 const std::vector<double> &Z_valence);
    void bcast(const bool am_i_mpi_rank0);
    void bcast_upf(const bool am_i_mpi_rank0);

    void set_atomic_position_crystal();

    void update_atomic_positions(const Eigen::VectorXd &position,
                                 const bool am_i_mpi_rank0,
                                 std::ostream *ost);
    void print_atomic_positions(std::ostream *ost) const;

    // for io_tc_files::read_crystal_structure()
    bool do_lattice_vectors_coincide(const Eigen::Matrix3d &lattice_vectors) const;
    void update_atomic_positions_by_file(const std::vector<Eigen::Vector3d> &atomic_position_cartesian,
                                         const bool am_i_mpi_rank0,
                                         std::ostream *ost);

};


#endif // TC_CRYSTAL_STRUCTURE_HPP

