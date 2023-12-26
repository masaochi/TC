// [class Symmetry]
// symmetry operations

#ifndef TC_SYMMETRY_HPP
#define TC_SYMMETRY_HPP

class Symmetry
{
private:
    int num_symmetry_operations_; // num. of symmetry operations
    // r -> (rotation_)*(r + translation_) (r: crystal coordinate)
    std::vector<Eigen::Matrix3i> rotation_; // rotation[nsym](0:2,0:2), in crystal coordinate
    std::vector<Eigen::Vector3d> translation_; // translation[nsym](0:2), in crystal coordinate

    // atom_index_after_sym[nsym][natoms]
    // e.g., "atom_index_after_sym[3][0] = 1" means that "Atom 0 will be moved to Atom 1 by symmetry operation 3"
    std::vector<std::vector<int> > atom_index_after_sym_; 

    std::vector<int> order_; // order_[nsym]. order_ = n such that S^n = identity for each symmetry operation S
    const int num_order_max_;

    bool is_time_reversal_symmetric_; // use time reversal symmetry or not (at present we can use it only for nospin)

public:
    int num_symmetry_operations() const { return num_symmetry_operations_; }
    const std::vector<Eigen::Matrix3i> &rotation() const { return rotation_; }
    const std::vector<Eigen::Vector3d> &translation() const { return translation_; }
    const std::vector<std::vector<int> > &atom_index_after_sym() const { return atom_index_after_sym_; }
    const std::vector<int> &order() const { return order_; }
    bool is_time_reversal_symmetric() const { return is_time_reversal_symmetric_; }

    Symmetry() : num_symmetry_operations_(0), num_order_max_(96) {}

    // set(): called from MPI rank=0
    void set(const Spin &spin, const CrystalStructure &crystal_structure,
             const int num_symmetry_operations,
             const std::vector<Eigen::Matrix3i> &rotation,
             const std::vector<Eigen::Vector3d> &translation,
             const std::string &no_t_rev, const std::string &noinv);
    void bcast(const bool am_i_mpi_rank0);
};


#endif // TC_SYMMETRY_HPP

