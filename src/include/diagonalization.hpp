// [class Diagonalization]
// scf_loop, convergence check, control paramters, etc.

#ifndef TC_DIAGONALIZATION_HPP
#define TC_DIAGONALIZATION_HPP

class Diagonalization
{
private:
    bool restarts_; // restart from the previous run
    double energy_tolerance_; // Total energy for scf convergence, sum of eigenvalues for band convergence (in Hartree)
    double charge_tolerance_; // for scf convergence (in e-)
    double force_tolerance_; // for structural optimization (in eV/ang) 
    int max_num_iterations_; // for scf
    int max_num_ionic_steps_; // for structural optimization

    // mixing parameters
    bool mixes_density_matrix_; // for scf. [true] mixes the density matrix [false] mixes the density
    double mixing_beta_; // for scf

    // Parameters for block-Davidson algorithm
    int num_refresh_david_; // update trial vectors by "num_refresh_david"-times for each update of the Fock operator // should >=1
    int max_num_blocks_david_; // determines a size of subspace dimension (see "diago_david_ndim" in QE)
    bool biortho_david_; // takes biorthogonal orbitals for Davidson diagonalization

    // Parameters for BFGS structural opt.
    const double size_ionic_displacement_;

    // for disk-io saving
    bool dumps_pwfn_; // (true) dumping pwfn.data (filenames.tc_pwfn_) (default: false)

    // for reading the crystal structure (mainly for restarting structural-optimization calculations, but can be used for normal SCF & BAND)
    bool reads_crystal_structure_; // (true) reading tc_crystal_structure.dat (filenames.tc_crystal_structure_) (default: false)

    // "return value = false" means that the orthonormalized new vector is zero.
    bool Gram_Schmidt(const std::vector<Eigen::VectorXcd> &vectors, const int &target_band,
                      Eigen::VectorXcd &target_vector) const;
    bool Gram_Schmidt_biortho(const std::vector<Eigen::VectorXcd> &vectors_right, 
                              const std::vector<Eigen::VectorXcd> &vectors_left,
                              const int &target_band,
                              Eigen::VectorXcd &target_vector_right,
                              Eigen::VectorXcd &target_vector_left) const;
    void sort_eigen(Eigen::VectorXcd &eigenvalues, Eigen::MatrixXcd &eigenvectors) const;
    void sort_eigen_overlap(Eigen::VectorXcd &eigenvalues, Eigen::MatrixXcd &eigenvectors) const;
    void make_a_new_trial_vector(const Eigen::VectorXd &kGvect2,
                                 Eigen::VectorXcd &phi,
                                 const Eigen::VectorXcd &H1phi,
                                 const Eigen::VectorXcd &H2phi,
                                 const Eigen::VectorXcd &H3phi,
                                 const Complex &eigenvalue,
                                 const bool uses_3body) const;
    void structural_optimization(const Method &method,
                                 CrystalStructure &crystal_structure,
                                 const Symmetry &symmetry,
                                 TotalEnergy &total_energy,
                                 const int &istep_ionic,
                                 const bool am_i_mpi_rank0,
                                 std::ostream *ost) const;

public:
    bool restarts() const { return restarts_; }
    double energy_tolerance() const { return energy_tolerance_; }
    double charge_tolerance() const { return charge_tolerance_; }
    double force_tolerance() const { return force_tolerance_; }
    int max_num_iterations() const { return max_num_iterations_; }
    int max_num_ionic_steps() const { return max_num_ionic_steps_; }
    bool mixes_density_matrix() const {return mixes_density_matrix_; }
    double mixing_beta() const { return mixing_beta_; }
    int num_refresh_david() const { return num_refresh_david_; }
    int max_num_blocks_david() const { return max_num_blocks_david_; }
    bool biortho_david() const { return biortho_david_; }
    bool dumps_pwfn() const { return dumps_pwfn_; }
    bool reads_crystal_structure() const { return reads_crystal_structure_; }

    Diagonalization() :
        restarts_(false),
        charge_tolerance_(1e-5),
        energy_tolerance_(1e-6),
        force_tolerance_(1e-2),
        max_num_iterations_(-1), // (SCF) 30, (BAND) 15, set in "void set()".
        max_num_ionic_steps_(0),
        mixes_density_matrix_(false),
        mixing_beta_(0.7),
        num_refresh_david_(1),
        max_num_blocks_david_(2),
        biortho_david_(false),
        size_ionic_displacement_(1e-3), // empirical parameter (in Bohr)
        dumps_pwfn_(false),
        reads_crystal_structure_(false) {}

    void set(const bool restarts,
             const double &energy_tolerance, const double &charge_tolerance,
             const double &force_tolerance,
             const int &max_num_iterations, const int &max_num_ionic_steps,
             const bool &mixes_density_matrix, const double &mixing_beta,
             const int &num_refresh_david, const int &max_num_blocks_david,
             const bool &biortho_david, const bool &dumps_pwfn, 
             const bool &reads_crystal_structure,
             const std::string &calc_mode, const std::string &calc_method, const bool is_heg);
    void bcast(const bool am_i_mpi_rank0);

    void scf(Parallelization &parallelization,
             const FileNames &file_names, MyClock &my_clock,
             const Method &method, CrystalStructure &crystal_structure,
             const Symmetry &symmetry, Potentials &potentials,
             const Spin &spin, const Kpoints &kpoints,
             PlaneWaveBasis &plane_wave_basis, 
             BlochStates &bloch_states, TotalEnergy &total_energy,
             std::ostream *ost) const;
    void band(Parallelization &parallelization,
              const FileNames &file_names, MyClock &my_clock,
              const Method &method, const CrystalStructure &crystal_structure,
              const Symmetry &symmetry, Potentials &potentials,
              const Spin &spin, const Kpoints &kpoints,
              PlaneWaveBasis &plane_wave_basis, 
              BlochStates &bloch_states, TotalEnergy &total_energy,
              std::ostream *ost) const;
    void block_davidson(Parallelization &parallelization, 
                        MyClock &my_clock, const Method &method,
                        const CrystalStructure &crystal_structure,
                        const Symmetry &symmetry, Potentials &potentials,
                        const Spin &spin, const Kpoints &kpoints,
                        PlaneWaveBasis &plane_wave_basis, 
                        BlochStates &bloch_states, TotalEnergy &total_energy,
                        std::ostream *ost) const;
};


#endif // TC_DIAGONALIZATION_HPP

