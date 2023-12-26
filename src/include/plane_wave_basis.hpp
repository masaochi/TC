// [class PlaneWaveBasis]
// plane-wave (PW) basis, list of PWs, FFT

#ifndef TC_PLANE_WAVE_BASIS_HPP
#define TC_PLANE_WAVE_BASIS_HPP

class PlaneWaveBasis
{
private:
    std::vector<int> num_G_at_k_scf_; // num_G_at_k_[irreducible k-point]: num. of plane waves (G) at each (irreducible) k-point
    std::vector<int> num_G_at_k_band_; // num_G_at_k_[irreducible k-point]: num. of plane waves (G) at each (irreducible) k-point

    // arrays[num_independent_spins][irreducible k-point][num. of equivalent k-points][num_G_at_k_]
    std::vector<std::vector<std::vector<Eigen::VectorXi> > > Gindex_at_k_scf_; // G-vector index (on the FFT grid) at each k-point
    std::vector<std::vector<std::vector<Eigen::VectorXi> > > Gindex_at_k_band_; // G-vector index (on the FFT grid) at each k-point
    std::vector<std::vector<std::vector<Eigen::VectorXcd> > > phase_factor_at_k_; // phase factor by translational operation (for SCF k-points)
    // phase_factor_at_k_band = 1.0 thus not required, since no symmetry operation is performed for band k-points. 

    std::vector<int> size_FFT_grid_vec_; // size_FFT_grid_vec_[0:2]: size of the plane-wave grid for FFT e.g. 24 x 24 x 12 then {24,24,12}
    int size_FFT_grid_; // size_FFT_grid_ = size_FFT_grid_vec_[0]*size_FFT_grid_vec_[1]*size_FFT_grid_vec_[2]

    bool gamma_only_; // memory saving in QE. need to make f(-G) from f(G).

    bool are_FFTarrays_initialized_;
    fftw_plan plan_forward_;
    fftw_plan plan_backward_;
    Complex* in_forward_;
    Complex* out_forward_;
    Complex* in_backward_;
    Complex* out_backward_;

    // called from set_Gvector_at_k()
    void set_G_from_mill(const int ipw_at_k, const int npw_at_k_org, const std::vector<int> &mill,
                         Eigen::Vector3i &Gvector, int &icount, int &zero_index);

public:
    const std::vector<int> &num_G_at_k_scf() const { return num_G_at_k_scf_; }
    const std::vector<int> &num_G_at_k_band() const { return num_G_at_k_band_; }
    const std::vector<std::vector<std::vector<Eigen::VectorXi> > > &Gindex_at_k_scf() const { return Gindex_at_k_scf_; }
    const std::vector<std::vector<std::vector<Eigen::VectorXi> > > &Gindex_at_k_band() const { return Gindex_at_k_band_; }
    const std::vector<std::vector<std::vector<Eigen::VectorXcd> > > &phase_factor_at_k() const { return phase_factor_at_k_; }
    const std::vector<int> &size_FFT_grid_vec() const { return size_FFT_grid_vec_; }
    int size_FFT_grid() const { return size_FFT_grid_; }
    bool gamma_only() const { return gamma_only_; }
    bool are_FFTarrays_initialized() const { return are_FFTarrays_initialized_; }

    PlaneWaveBasis() : size_FFT_grid_(0), gamma_only_(false), are_FFTarrays_initialized_(false) {}
    PlaneWaveBasis(const PlaneWaveBasis &pw)
    {
        num_G_at_k_scf_ = pw.num_G_at_k_scf();
        num_G_at_k_band_ = pw.num_G_at_k_band();
        Gindex_at_k_scf_ = pw.Gindex_at_k_scf();
        Gindex_at_k_band_ = pw.Gindex_at_k_band();
        phase_factor_at_k_ = pw.phase_factor_at_k();
        gamma_only_ = pw.gamma_only();

        are_FFTarrays_initialized_ = false;
        if (pw.are_FFTarrays_initialized())
        {
            setup_FFT(pw.size_FFT_grid_vec()[0], pw.size_FFT_grid_vec()[1], pw.size_FFT_grid_vec()[2]);
        }
    }

    PlaneWaveBasis& operator=(const PlaneWaveBasis &pw)
    {
        num_G_at_k_scf_ = pw.num_G_at_k_scf();
        num_G_at_k_band_ = pw.num_G_at_k_band();
        Gindex_at_k_scf_ = pw.Gindex_at_k_scf();
        Gindex_at_k_band_ = pw.Gindex_at_k_band();
        phase_factor_at_k_ = pw.phase_factor_at_k();
        gamma_only_ = pw.gamma_only();

        are_FFTarrays_initialized_ = false;
        if (pw.are_FFTarrays_initialized())
        {
            setup_FFT(pw.size_FFT_grid_vec()[0], pw.size_FFT_grid_vec()[1], pw.size_FFT_grid_vec()[2]);
        }
        return *this;
    }

    // set the size of num_G_at_k_ and Gvector_at_k_
    void resize_G_at_k(const int num_independent_spins, const int num_irreducible_kpoints, const std::string &calc_mode);
    void set_num_G_at_k(const std::vector<int> num_G_at_k, const std::string &calc_mode);
    void set_Gvector_at_k(const Symmetry &symmetry, const Kpoints &kpoints,
                          const int ispin, const int ik, const std::vector<int> &mill,
                          int &zero_index, const std::string &calc_mode); // mill = QE variable
    void set_scfinfo(const std::vector<int> &num_G_at_k_scf,
                     const std::vector<std::vector<std::vector<Eigen::VectorXi> > > &Gindex_at_k_scf,
                     const std::vector<std::vector<std::vector<Eigen::VectorXcd> > > &phase_factor_at_k);
        
    // set size_FFT_grid_vec_, size_FFT_grid_, and FFT variables
    void setup_FFT(const int nr1, const int nr2, const int nr3);
    void set_gamma_only(const bool gamma_only);

    // ipw -> (i1, i2, i3)
    Eigen::Vector3i get_Gvector(const int &ipw) const;

    // convert phi on "reduced G-vector grid on each k-point" to orbital on "FFT grid"
    // non-collinear calc. not supported
    void get_orbital_FFTgrid(const int &ispin, const int &ik, const int &isym,
                             const bool is_time_reversal_used_at_k,
                             const Eigen::VectorXcd &phi,
                             Eigen::VectorXcd &orbital_FFTgrid,
                             const std::string &calc_mode) const;

    void FFT_forward(const Eigen::VectorXcd &in, Eigen::VectorXcd &out);
    void FFT_backward(const Eigen::VectorXcd &in, Eigen::VectorXcd &out);

    void bcast(const std::string &calc_mode, const bool am_i_mpi_rank0,
               const bool setup_FFTgrid);
    ~PlaneWaveBasis();
};


#endif // TC_PLANE_WAVE_BASIS_HPP

