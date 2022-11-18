// [class BlochStates]
// wave functions, eigenvalues, filling, num. of bands, etc.

#ifndef TC_BLOCH_STATES_HPP
#define TC_BLOCH_STATES_HPP

class BlochStates
{
private:
    double num_electrons_; // including spin degrees of freedom
    double fermi_energy_; // in Hartree.

    // num_bands is array: (spin-polarized && non-collinear) num_bands[0]=up, num_bands[1]=down. (otherwise) num_bands[0].
    std::vector<int> num_bands_qe_; // num. of bands in Quantum Espresso
    // num. of bands in TC++ for the SCF or BAND orbitals, should satisfy "num_bands_scf/band_ <= num_bands_qe_" 
    std::vector<int> num_bands_scf_; 
    std::vector<int> num_bands_band_; 

    // num_occupied_bands_[spin][kpoint]
    std::vector<std::vector<int> > num_occupied_bands_; // num. of bands where filling_ > 1e-8

    // eigenvalues[spin][kpoint][band]. Num. of spin indicies = 2 for spin-polarized && non-collinear, otherwise 1. (same as num_bands)
    std::vector<std::vector<std::vector<double> > > eigenvalues_qe_;
    std::vector<std::vector<std::vector<Complex> > > eigenvalues_scf_; // Note: TC Hamiltonian is non-Hermitian
    std::vector<std::vector<std::vector<Complex> > > eigenvalues_band_; // Note: TC Hamiltonian is non-Hermitian
    std::vector<std::vector<std::vector<double> > > filling_; // 0<=filling<=2 for spin-unpolarized && collinear. otherwise, 0<=filling<=1.

    double density_difference_; // L1 norm of the density difference between this and previous SCF loops (to check convergence)

    // phik[spin][k-point][band][num_spinor][npw]: wave function in G-space at each irreducible k-point
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > phik_scf_;
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > phik_left_scf_; // for method.calc_method_=="BITC"
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > phik_band_;
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > phik_left_band_; // for method.calc_method_=="BITC"
    std::vector<Eigen::VectorXcd> density_; // 0: up, 1: down. in R-space. Note: complex for BITC method.

    // called in set_qe_xml()
    void set_num_bands_qe(const std::vector<int> &num_bands_qe,
                          const std::string &calc_mode);
    void set_eigenvalues_qe(const std::vector<std::vector<double> > &eigenvalues,
                            const std::string &calc_mode);

    void resize_filling(const int &num_independent_spins,
                        const int &num_irreducible_kpoints_scf);
    void resize_num_occupied_bands(const int &num_independent_spins,
                                   const int &num_irreducible_kpoints_scf);

    // called in set_density()
    void integrate_local_density(const CrystalStructure &crystal_structure,
                                 const Kpoints &kpoints, const PlaneWaveBasis &plane_wave_basis,
                                 std::ostream *ost);

public:
    double num_electrons() const { return num_electrons_; }
    double fermi_energy() const { return fermi_energy_; }
    const std::vector<int> &num_bands_qe() const { return num_bands_qe_; }
    const std::vector<int> &num_bands_scf() const { return num_bands_scf_; }
    const std::vector<int> &num_bands_band() const { return num_bands_band_; }
    const std::vector<std::vector<int> > &num_occupied_bands() const { return num_occupied_bands_; }
    const std::vector<std::vector<std::vector<double> > > &eigenvalues_qe() const { return eigenvalues_qe_; }
    const std::vector<std::vector<std::vector<Complex> > > &eigenvalues_scf() const { return eigenvalues_scf_; }
    const std::vector<std::vector<std::vector<Complex> > > &eigenvalues_band() const { return eigenvalues_band_; }
    const std::vector<std::vector<std::vector<double> > > &filling() const { return filling_; }
    double density_difference() const { return density_difference_; }
    const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_scf() const { return phik_scf_; }
    const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_left_scf() const { return phik_left_scf_; }
    const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_band() const { return phik_band_; }
    const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_left_band() const { return phik_left_band_; }
    const std::vector<Eigen::VectorXcd> &density() const { return density_; }

    BlochStates() : num_electrons_(-1.0), num_bands_qe_{-1}, num_bands_scf_{-1}, num_bands_band_{-1}, 
        density_difference_(0.0) {}

    // [bloch_states_initialize.cpp]
    // set_num_bands_tc() &
    // set_qe_xml() & resize_phik_qe() & set_phik_qe() & bcast_qe_input() &
    // bcast_scfinfo() & set_phik_scf_from_tc_wfc() & set_eigenvalues_scf_from_tc_energy()

    void set_num_bands_tc(const std::vector<int> &num_bands_tc,
                          const std::string &calc_mode);
    void set_qe_xml(const std::vector<int> &num_bands_qe,
                    const std::vector<std::vector<double> > &eigenvalues,
                    const std::string &calc_mode,
                    const double &num_electrons);
    void resize_phik_qe(const PlaneWaveBasis &plane_wave_basis,
                        const bool is_spinor, const int &num_independent_spins,
                        const std::string &calc_mode);
    void set_phik_qe(const int ispin, const int ik, 
                     const std::vector<Complex> &evc,
                     const std::string &calc_mode);
    void bcast_qe_input(const PlaneWaveBasis &plane_wave_basis,
                        const bool is_spinor,
                        const std::string &calc_mode,
                        const std::string &calc_method,
                        const bool am_i_mpi_rank0);

    // called in BAND calculation
    void bcast_scfinfo(const PlaneWaveBasis &plane_wave_basis,
                       const bool is_spinor,
                       const std::string &calc_method,
                       const bool am_i_mpi_rank0);

    // set phik(_left)_"calc_mode"_ by reading tc_wfc.dat
    void set_phik_from_tc_wfc(const int &ispin, const int &ik, const int &iband,
                              const Eigen::VectorXcd &wfc_temp,
                              const std::string &right_or_left,
                              const std::string &calc_mode);
    // set eigenvalues_"calc_mode"_ by reading tc_energy.dat
    void set_eigenvalues_from_tc_energy(const int &ispin, const int &ik,
                                        const Eigen::VectorXcd &eigenvalues_temp,
                                        const std::string &calc_mode);

    // [bloch_states_scfloop_filling_density.cpp]
    // set_filling() & set_density()

    // set fermi_energy_ & filling_ & num_occupied_bands_
    // including print eigenvalues_scf_
    void set_filling(const Spin &spin, const Kpoints &kpoints, 
                     const bool am_i_mpi_rank0, std::ostream *ost);

    // (is_first_iter==false) density mixing is performed.
    // non-collinear calculation is not supported
    void set_density(const CrystalStructure &crystal_structure,
                     const Kpoints &kpoints, PlaneWaveBasis &plane_wave_basis,
                     const double &mixing_beta,
                     const bool is_first_iter, const bool is_bitc,
                     const bool am_i_mpi_rank0, std::ostream *ost);

    // [bloch_states_scfloop_phik.cpp]
    // reset_phik() & bcast_phik()

    // reset phik_scf/band in diagonalization
    // non-collinear calculation is not supported
    void reset_phik(const int &ispin, const int &ik,
                    const std::vector<Eigen::VectorXcd> &V,
                    const std::string &right_or_left,
                    const std::string &calc_mode);
    void bcast_phik(const bool bcast_phik_left,
                    const std::string &calc_mode);

    // [bloch_states_scfloop_eigenvalues.cpp] 
    // print_band_energies() &
    // reset_eigenvalues() & bcast_eigenvalues()

    // print eigenvalues_band_
    void print_band_energies(const Kpoints &kpoints, std::ostream *ost) const;

    // reset eigenvalues_scf/band in diagonalization
    void reset_eigenvalues(const int &ispin, const int &ik,
                           const Eigen::VectorXcd &eigenvalues,
                           const std::string &calc_mode);
    void bcast_eigenvalues(const std::string &calc_mode);

};


#endif // TC_BLOCH_STATES_HPP

