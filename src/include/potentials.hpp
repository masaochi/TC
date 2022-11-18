// [class Potentials]
// pseudopotential, Jastrow, Ewald sum

// given in atomic unit (Hartree)

#ifndef TC_POTENTIALS_HPP
#define TC_POTENTIALS_HPP

class Potentials
{
private:
    bool is_heg_; // (true) homogeneous-electron-gas mode, where pseudopot. and ewald_energy_ are neglected

    // pseudopotential (pp) data
    Eigen::VectorXd Vpp_local_; // local part in R-space (on FFT grid)

    // pseudopotential (pp) data for each atomic species using upf-radial grid
    std::vector<std::vector<double> > pseudo_rmesh_; // radial grid (r value) [iatomic_species][imesh]
    std::vector<std::vector<double> > pseudo_weight_rmesh_; // weight for integration on the radial grid [iatomic_species][imesh]
    std::vector<std::vector<std::vector<double> > > pseudo_beta_; // Kleinman-Bylander projectors [iatomic_species][iproj][imesh]

    // pseudopotential (pp) data for each atomic species
    std::vector<std::vector<int> > pseudo_lbeta_; // angular momentum for pseudo_beta_ (see below)
    std::vector<Eigen::MatrixXd> pseudo_dij_; // coefficients of Kleinman-Bylander projectors. |i>dij[iatomic_species]<j|

    // for divegence correction (Coulomb potential)
    bool includes_div_correction_; // Gygi-Baldereschi treatment of the potential divergencies in q-space
    double alpha_Vaux_; // 4pi*exp[-alpha_Vaux_*q^2]/q^2 is used
    std::vector<double> sum_of_Vaux_scf_; // [ik: irreducible kpoints on the SCF k-mesh]
    std::vector<double> sum_of_Vaux_band_; // [ik: on the band k-mesh]
    std::vector<Eigen::Vector3d> sum_of_kVaux_scf_; // [ik: irreducible kpoints on the SCF k-mesh](idim). sum of \nabla Vaux (set only for zero-weight k-points)
    std::vector<Eigen::Vector3d> sum_of_kVaux_band_; // [ik: on the band k-mesh](idim). sum of \nabla Vaux

    // called in set_upf()
    void set_unit(std::vector<std::vector<double> > &pseudo_local,
                  std::vector<Eigen::MatrixXd> &pseudo_dij);
    void set_Vpp_local(const std::vector<double> &Z_valence,
                       const std::vector<std::vector<double> > &pseudo_local,
                       const CrystalStructure &crystal_structure,
                       PlaneWaveBasis &plane_wave_basis);

    // Simpson's rule integration on the pseudo-potential radial mesh. Integrate "function".
    double simpson(const std::vector<double> &function, const int &iatomic_species) const;

    // called in set_sum_of_Vaux(). set sum_of_Vaux_scf/band_
    void set_sum_of_Vaux_each(const Parallelization &parallelization,
                              const CrystalStructure &crystal_structure,
                              const Kpoints &kpoints,
                              const std::vector<std::vector<Eigen::Vector3d> > &kvectors,
                              std::vector<double> &sum_of_Vaux_);

    // called in set_sum_of_Vaux(). set sum_of_kVaux_scf_/band_
    void set_sum_of_kVaux_each(const Parallelization &parallelization,
                               const CrystalStructure &crystal_structure,
                               const Kpoints &kpoints,
                               const std::vector<std::vector<Eigen::Vector3d> > &kvectors,
                               std::vector<Eigen::Vector3d> &sum_of_kVaux_);

    // set parallelization index (used only in set_sum_of_Vaux_each() and set_sum_of_kVaux_band)
    // cf. Parallelization::set_orbital_parallelization()
    void set_parallelization_in_sum_Vaux(const Parallelization &parallelization,
                                         const Kpoints &kpoints,
                                         const std::vector<std::vector<Eigen::Vector3d> > &kvectors,
                                         std::vector<bool> &if_match,
                                         int &ibegin, int &iend);

public:
    Jastrow jastrow; // Jatrow parameters (class Jastrow)

    bool is_heg() const { return is_heg_; }
    const Eigen::VectorXd &Vpp_local() const { return Vpp_local_; }
    const std::vector<std::vector<double> > &pseudo_rmesh() const { return pseudo_rmesh_; }
    const std::vector<std::vector<double> > &pseudo_weight_rmesh() const { return pseudo_weight_rmesh_; }
    const std::vector<std::vector<std::vector<double> > > &pseudo_beta() const { return pseudo_beta_; }
    const std::vector<std::vector<int> > &pseudo_lbeta() const { return pseudo_lbeta_; }
    const std::vector<Eigen::MatrixXd> &pseudo_dij() const { return pseudo_dij_; }
    bool includes_div_correction() const { return includes_div_correction_; }
    double alpha_Vaux() const { return alpha_Vaux_; }
    const std::vector<double> &sum_of_Vaux_scf() const { return sum_of_Vaux_scf_; }
    const std::vector<double> &sum_of_Vaux_band() const { return sum_of_Vaux_band_; }
    const std::vector<Eigen::Vector3d> &sum_of_kVaux_band() const { return sum_of_kVaux_band_; }
    const std::vector<Eigen::Vector3d> &sum_of_kVaux_scf() const { return sum_of_kVaux_scf_; }

    Potentials() :
        is_heg_(false), 
        includes_div_correction_(true),
        alpha_Vaux_(0.0) {}

    void set_is_heg(const bool is_heg) { is_heg_ = is_heg; }
    void set_includes_div_correction(const bool includes_div_correction) { includes_div_correction_ = includes_div_correction; }
    void bcast_tc_input() 
    {
        MPI_Bcast(&is_heg_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
        MPI_Bcast(&includes_div_correction_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    }

    void set_upf(const std::vector<double> &Z_valence,
                 const std::vector<std::vector<double> > &rmesh,
                 const std::vector<std::vector<double> > &weight_rmesh,
                 std::vector<std::vector<double> > &pseudo_local,
                 const std::vector<std::vector<std::vector<double> > > &pseudo_beta,
                 const std::vector<std::vector<int> > &pseudo_lbeta,
                 std::vector<Eigen::MatrixXd> &pseudo_dij,
                 const CrystalStructure &crystal_structure,
                 PlaneWaveBasis &plane_wave_basis);
    void bcast_upf(const bool am_i_mpi_rank0);

    std::vector<Eigen::VectorXcd> projector_nonlocal(const CrystalStructure &crystal_structure,
                                                     const PlaneWaveBasis &plane_wave_basis,
                                                     const Eigen::VectorXi &Gindex_at_k,
                                                     const Eigen::Vector3d &kvector,
                                                     const std::vector<std::vector<Eigen::VectorXcd> > &ylm,
                                                     const int &ispin, const int &ik,
                                                     const int &iatomic_species, const int &iproj) const;

    // for divergence correction (Coulomb potential)
    void set_sum_of_Vaux(const Parallelization &parallelization,
                         const Method &method,
                         const CrystalStructure &crystal_structure,
                         const Kpoints &kpoints);
};

#endif // TC_POTENTIALS_HPP

