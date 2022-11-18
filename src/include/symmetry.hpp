// [class Symmetry]
// symmetry operations

#ifndef TC_SYMMETRY_HPP
#define TC_SYMMETRY_HPP

class Symmetry
{
private:
    int num_symmetry_operations_; // num. of symmetry operations
    // r -> (rotation_)*r + translation_ (r: crystal coordinate)
    std::vector<Eigen::Matrix3i> rotation_; // rotation[nsym](0:2,0:2), in crystal coordinate
    std::vector<Eigen::Vector3d> translation_; // translation[nsym](0:2), in crystal coordinate

    bool is_time_reversal_symmetric_; // use time reversal symmetry or not (at present we can use it only for nospin)

public:
    int num_symmetry_operations() const { return num_symmetry_operations_; }
    const std::vector<Eigen::Matrix3i> &rotation() const { return rotation_; }
    const std::vector<Eigen::Vector3d> &translation() const { return translation_; }
    bool is_time_reversal_symmetric() const { return is_time_reversal_symmetric_; }

    Symmetry() : num_symmetry_operations_(0) {}

    // set(): called from MPI rank=0
    void set(const Spin &spin, const int num_symmetry_operations,
             const std::vector<Eigen::Matrix3i> &rotation,
             const std::vector<Eigen::Vector3d> &translation,
             const std::string &no_t_rev, const std::string &noinv);
    void bcast(const bool am_i_mpi_rank0);
};


#endif // TC_SYMMETRY_HPP

