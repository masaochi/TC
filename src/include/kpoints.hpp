// [class Kpoints]
// list of k-point, indicies of symmetry operations for k-points

#ifndef TC_KPOINTS_HPP
#define TC_KPOINTS_HPP

class Kpoints
{
private:
    int num_irreducible_kpoints_scf_; // num. of irreducible k-points for the SCF k-mesh
    int num_irreducible_kpoints_band_; // since band k-points do not use symmetry, this is equivalent to num. kpoints for bands
    int num_kpoints_; // total num. of non-zero-weight k-points (i.e. k-points on the SCF k-mesh. 4x4x4 then 64)

    int num_kpoints_all_scf_; // total num. of all scf k-points including zero-weight ones (cf. fakescf calculation with band k-points)
    std::vector<std::vector<int> > index_all_kscf_; // index_all_kscf_[iq_isymq][0] = iq, index_all_kscf_[iq_isymq][1] = isymq for 0 <= iq_isymq < num_kpoints_all_scf_.

    // k-point weight for irreducible k-points (QE variable)
    // summation of kweight_ = 2 for no-spin, 1 for spin-polarized or non-collinear calculations.
    // Information of kweight is taken into account by bloch_states.filling_ in TC++ calculation.
    std::vector<double> kweight_scf_; // for SCF k-mesh. Note: kweight_band_ is not necessary. kweight read from the QE xml file is not used for BAND calc.

    // kvec[ik][is](0:2) = k-vector. ik = irreducible k-index, is = symmetry index (0 = irreducible k-point)
    // crystal coordinate. e.g., (1.0, 0.0, 0.0) is on the G-grid and equivalent to the Gamma point
    std::vector<std::vector<Eigen::Vector3d> > kvectors_scf_;
    std::vector<std::vector<Eigen::Vector3d> > kvectors_band_;

    // index of symmetry operation for making (ik,is) k-vector from (ik,0) irreducible k-vector
    // k'(ik,is) = U^+ k (ik,0) (is_time_reversal_used_at_k_==false)
    // k'(ik,is) = -(U^+ k(ik,0)) (is_time_reversal_used_at_k_==true)
    // Note: k = SCF k-points, since no symmetry operation is performed for band k-points.
    std::vector<std::vector<int> > index_of_rotation_at_k_; // index of rotation operator in symmetries.rotation
    std::vector<std::vector<bool> > is_time_reversal_used_at_k_; // time reveral symmetry is used ("true") or not ("false")

    // variables for smearing
    std::string smearing_mode_; // fixed or gaussian (default: gaussian)
    double smearing_width_; // in Hartree. ot used when smearing_mode_ = fixed. (default: 0.01)

    // called in set() 
    void set_kvectors_by_symmetry(const Spin &spin, const Symmetry &symmetry, std::ostream *ost);
    void show_kvectors(const Symmetry &symmetry, std::ostream *ost) const;

public:
    int num_irreducible_kpoints_scf() const { return num_irreducible_kpoints_scf_; }
    int num_irreducible_kpoints_band() const { return num_irreducible_kpoints_band_; }
    int num_kpoints() const { return num_kpoints_; }
    int num_kpoints_all_scf() const { return num_kpoints_all_scf_; }
    const std::vector<std::vector<int> > &index_all_kscf() const { return index_all_kscf_; }
    const std::vector<double> &kweight_scf() const { return kweight_scf_; }
    const std::vector<std::vector<Eigen::Vector3d> > &kvectors_scf() const { return kvectors_scf_; }
    const std::vector<std::vector<Eigen::Vector3d> > &kvectors_band() const { return kvectors_band_; }
    const std::vector<std::vector<int> > &index_of_rotation_at_k() const { return index_of_rotation_at_k_; }
    const std::vector<std::vector<bool> > &is_time_reversal_used_at_k() const { return is_time_reversal_used_at_k_; }
    const std::string &smearing_mode() const { return smearing_mode_; }
    double smearing_width() const { return smearing_width_; }

    Kpoints() : num_irreducible_kpoints_scf_(-1), num_irreducible_kpoints_band_(-1),
                num_kpoints_(-1), num_kpoints_all_scf_(-1), smearing_mode_("gaussian"), smearing_width_(0.01) {}

    void set_tc_input(const std::string &smearing_mode, const double &smearing_width);
    void bcast_tc_input(const bool am_i_mpi_rank0);

    void set_qe_input(const Spin &spin, const Symmetry &symmetry,
                      const int num_irreducible_kpoints,
                      const std::vector<double> kweight,
                      const std::vector<Eigen::Vector3d> &kvectors_irred,
                      const std::string &calc_mode,
                      std::ostream *ost);
    void bcast_qe_input(const std::string &calc_mode,
                        const bool am_i_mpi_rank0);

    // initialize SCF variables after reading tc_scfinfo.dat (called in BAND calculation)
    // Bcast is done by bcast_qe_input().
    void set_scfinfo(const std::vector<double> &kweight_scf,
                     const std::vector<std::vector<Eigen::Vector3d> > &kvectors_scf,
                     const std::vector<std::vector<int> > &index_of_rotation_at_k,
                     const std::vector<std::vector<bool> > &is_time_reversal_used_at_k);

    // return band filling including the spin factor
    double return_band_filling(const Spin &spin, const double &orbital_energy, const double &fermi_energy,
                               const int &num_electrons, const int &iband) const;
};


#endif // TC_KPOINTS_HPP

