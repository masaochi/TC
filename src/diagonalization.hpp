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
    int max_num_iterations_; // for scf

    bool mixes_density_matrix_; // for scf. [true] mixes the density matrix [false] mixes the density
    double mixing_beta_; // for scf

    // Parameters for block-Davidson algorithm
    int num_refresh_david_; // update trial vectors by "num_refresh_david"-times for each update of the Fock operator // should >=1
    int max_num_blocks_david_; // determines a size of subspace dimension (see "diago_david_ndim" in QE)

    // "return value = false" means that the orthonormalized new vector is zero.
    bool Gram_Schmidt(std::vector<Eigen::VectorXcd> &vectors, const int &target_band,
                      Eigen::VectorXcd &target_vector);
    void sort_eigen(Eigen::VectorXcd &eigenvalues, Eigen::MatrixXcd &eigenvectors);
    void make_a_new_trial_vector(const Eigen::VectorXd &kGvect2,
                                 Eigen::VectorXcd &phi,
                                 const Eigen::VectorXcd &H1phi,
                                 const Eigen::VectorXcd &H2phi,
                                 const Eigen::VectorXcd &H3phi,
                                 const Complex &eigenvalue,
                                 const bool uses_3body);

public:
    bool restarts() const { return restarts_; }
    double energy_tolerance() const { return energy_tolerance_; }
    double charge_tolerance() const { return charge_tolerance_; }
    int max_num_iterations() const { return max_num_iterations_; }
    bool mixes_density_matrix() const {return mixes_density_matrix_; }
    double mixing_beta() const { return mixing_beta_; }
    int num_refresh_david() const { return num_refresh_david_; }
    int max_num_blocks_david() const { return max_num_blocks_david_; }

    Diagonalization() :
        restarts_(false),
        charge_tolerance_(1e-4),
        energy_tolerance_(1e-5),
        max_num_iterations_(-1), // (SCF) 30, (BAND) 15, set in "void set()".
        mixes_density_matrix_(false),
        mixing_beta_(0.7),
        num_refresh_david_(1),
        max_num_blocks_david_(2) {}

    void set(const bool restarts,
             const double &energy_tolerance, const double &charge_tolerance,
             const int &max_num_iterations,
             const bool &mixes_density_matrix, const double &mixing_beta,
             const int &num_refresh_david, const int &max_num_blocks_david,
             const std::string &calc_mode, const bool is_heg);
    void bcast(const bool am_i_mpi_rank0);

    void scf(Parallelization &parallelization,
             const FileNames &file_names, MyClock &my_clock,
             const Method &method, const CrystalStructure &crystal_structure,
             const Symmetry &symmetry, const Potentials &potentials,
             const Spin &spin, const Kpoints &kpoints,
             PlaneWaveBasis &plane_wave_basis, 
             BlochStates &bloch_states, TotalEnergy &total_energy,
             std::ostream *ost);
    void band(Parallelization &parallelization,
              const FileNames &file_names, MyClock &my_clock,
              const Method &method, const CrystalStructure &crystal_structure,
              const Symmetry &symmetry, const Potentials &potentials,
              const Spin &spin, const Kpoints &kpoints,
              PlaneWaveBasis &plane_wave_basis, 
              BlochStates &bloch_states, TotalEnergy &total_energy,
              std::ostream *ost);
    void block_davidson(Parallelization &parallelization, 
                        MyClock &my_clock, const Method &method,
                        const CrystalStructure &crystal_structure,
                        const Symmetry &symmetry, const Potentials &potentials,
                        const Spin &spin, const Kpoints &kpoints,
                        PlaneWaveBasis &plane_wave_basis, 
                        BlochStates &bloch_states, TotalEnergy &total_energy,
                        std::ostream *ost);
};


#endif // TC_DIAGONALIZATION_HPP

