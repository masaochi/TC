// [class PlaneWaveBasis]
// plane-wave (PW) basis, list of PWs, FFT

#include "include/header.hpp"

void PlaneWaveBasis::resize_G_at_k(const int num_independent_spins, const int num_irreducible_kpoints,
                                   const std::string &calc_mode)
{ 
    assert(num_independent_spins > 0);
    assert(num_irreducible_kpoints > 0);
    assert(calc_mode=="SCF" || calc_mode=="BAND");

    std::vector<int> &num_G_at_k_ref_ =
        calc_mode=="SCF" ? num_G_at_k_scf_ : num_G_at_k_band_;
    num_G_at_k_ref_.resize(num_irreducible_kpoints);

    std::vector<std::vector<std::vector<Eigen::VectorXi> > > &Gindex_at_k_ref_ =
        calc_mode=="SCF" ? Gindex_at_k_scf_ : Gindex_at_k_band_;
    Gindex_at_k_ref_.resize(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++) { Gindex_at_k_ref_[ispin].resize(num_irreducible_kpoints); }

    if (calc_mode=="SCF") // related to the symmetry operation, not used for band k-points
    {
        phase_factor_at_k_.resize(num_independent_spins);
        for (int ispin=0; ispin<num_independent_spins; ispin++) { phase_factor_at_k_[ispin].resize(num_irreducible_kpoints); }
    }
}

void PlaneWaveBasis::set_num_G_at_k(const std::vector<int> num_G_at_k, const std::string &calc_mode) 
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    std::vector<int> &num_G_at_k_ref_ =
        calc_mode=="SCF" ? num_G_at_k_scf_ : num_G_at_k_band_;

    assert(num_G_at_k_ref_.size() == num_G_at_k.size());
    num_G_at_k_ref_ = num_G_at_k;
}

void PlaneWaveBasis::set_Gvector_at_k(const Symmetry &symmetry, const Kpoints &kpoints,
                                      const int ispin, const int ik, const std::vector<int> &mill,
                                      const std::string &calc_mode)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    const std::vector<int> &num_G_at_k_ref_ =
        calc_mode=="SCF" ? num_G_at_k_scf_ : num_G_at_k_band_;
    std::vector<std::vector<std::vector<Eigen::VectorXi> > > &Gindex_at_k_ref_ =
        calc_mode=="SCF" ? Gindex_at_k_scf_ : Gindex_at_k_band_;

    assert(ispin >= 0);
    assert(ispin < Gindex_at_k_ref_.size());

    assert(ik >= 0);
    assert(ik < Gindex_at_k_ref_[ispin].size());

    const int nsym = calc_mode=="SCF" ? kpoints.kvectors_scf()[ik].size() : 1; // num. of equivalent k-points
    const int npw_at_k = num_G_at_k_ref_[ik]; // num. of plane waves

    Gindex_at_k_ref_[ispin][ik].resize(nsym);    
    for (int isym=0; isym<nsym; isym++) { Gindex_at_k_ref_[ispin][ik][isym].resize(npw_at_k); }

    if (calc_mode=="SCF") 
    {
        phase_factor_at_k_[ispin][ik].resize(nsym);
        for (int isym=0; isym<nsym; isym++) { phase_factor_at_k_[ispin][ik][isym].resize(npw_at_k); }
    }

    Eigen::Vector3i Gvector, Gvector_sym;
    if (calc_mode=="SCF")
    {
        int icount = 0;
        for (int ipw_at_k=0; ipw_at_k<npw_at_k ; ipw_at_k++) // G-grid at k (smaller than the FFT grid)
        {
            for (int idim=0; idim<3; idim++)
            {
                Gvector(idim) = mill[icount]; // mill(1:3,1:npw_at_k) in the fortran code (see read_qe_wfc.f90)
                ++icount;

                while (Gvector(idim) <= -(size_FFT_grid_vec_[idim]+1)/2) { Gvector(idim) += size_FFT_grid_vec_[idim]; }
                while (Gvector(idim) > size_FFT_grid_vec_[idim]/2) { Gvector(idim) -= size_FFT_grid_vec_[idim]; }
            }
            
            for (int isym=0; isym<nsym; isym++)
            {
                int ind_sym = kpoints.index_of_rotation_at_k()[ik][isym];
                
                Gvector_sym = symmetry.rotation()[ind_sym].transpose() * Gvector; // U^+ G
                if (kpoints.is_time_reversal_used_at_k()[ik][isym]) { Gvector_sym = -Gvector_sym; }
                double Gr = symmetry.translation()[ind_sym].dot(Gvector_sym.cast<double>()); // G_sym*r_0
                
                for (int idim=0; idim<3; idim++)
                {
                    while (Gvector_sym(idim)<0) { Gvector_sym(idim) += size_FFT_grid_vec_[idim]; }
                    while (Gvector_sym(idim)>=size_FFT_grid_vec_[idim]) { Gvector_sym(idim) -= size_FFT_grid_vec_[idim]; }
                }
                Gindex_at_k_scf_[ispin][ik][isym](ipw_at_k) =
                    Gvector_sym(0) 
                    + Gvector_sym(1)*size_FFT_grid_vec_[0]
                    + Gvector_sym(2)*size_FFT_grid_vec_[0]*size_FFT_grid_vec_[1];
                phase_factor_at_k_[ispin][ik][isym](ipw_at_k) = std::cos(2*PI*Gr) + I*std::sin(2*PI*Gr);
            } // isym
        } // ipw_at_k
    }
    else // BAND
    {
        int icount = 0;
        for (int ipw_at_k=0; ipw_at_k<npw_at_k ; ipw_at_k++) // G-grid at k (smaller than the FFT grid)
        {
            for (int idim=0; idim<3; idim++)
            {
                Gvector(idim) = mill[icount]; // mill(1:3,1:npw_at_k) in the fortran code (see read_qe_wfc.f90)
                ++icount;

                while (Gvector(idim)<0) { Gvector(idim) += size_FFT_grid_vec_[idim]; }
                while (Gvector(idim)>=size_FFT_grid_vec_[idim]) { Gvector(idim) -= size_FFT_grid_vec_[idim]; }
            }
            Gindex_at_k_band_[ispin][ik][0](ipw_at_k) =
                Gvector(0) 
                + Gvector(1)*size_FFT_grid_vec_[0]
                + Gvector(2)*size_FFT_grid_vec_[0]*size_FFT_grid_vec_[1];
        } // ipw_at_k
    } // calc_mode
}

void PlaneWaveBasis::set_scfinfo(const std::vector<int> &num_G_at_k_scf,
                                 const std::vector<std::vector<std::vector<Eigen::VectorXi> > > &Gindex_at_k_scf,
                                 const std::vector<std::vector<std::vector<Eigen::VectorXcd> > > &phase_factor_at_k)
{
    num_G_at_k_scf_ = num_G_at_k_scf;
    Gindex_at_k_scf_ = Gindex_at_k_scf;
    phase_factor_at_k_ = phase_factor_at_k;
}

void PlaneWaveBasis::setup_FFT(const int nr1, const int nr2, const int nr3)
{
    assert(are_FFTarrays_initialized_ == false); // not initialized twice!
    assert(nr1 > 0 && nr2 > 0 && nr3 > 0);

    size_FFT_grid_vec_.resize(3);
    size_FFT_grid_vec_ = {nr1, nr2, nr3};
    size_FFT_grid_ = nr1*nr2*nr3;
    
    in_forward = (Complex*) fftw_malloc(sizeof(Complex)*size_FFT_grid_);
    out_forward = (Complex*) fftw_malloc(sizeof(Complex)*size_FFT_grid_);
    in_backward = (Complex*) fftw_malloc(sizeof(Complex)*size_FFT_grid_);
    out_backward = (Complex*) fftw_malloc(sizeof(Complex)*size_FFT_grid_);
    
    plan_forward = fftw_plan_dft_3d(size_FFT_grid_vec_[2],size_FFT_grid_vec_[1],size_FFT_grid_vec_[0],
                                    reinterpret_cast<fftw_complex*>(in_forward),
                                    reinterpret_cast<fftw_complex*>(out_forward),
                                    FFTW_FORWARD,FFTW_PATIENT);
    plan_backward = fftw_plan_dft_3d(size_FFT_grid_vec_[2],size_FFT_grid_vec_[1],size_FFT_grid_vec_[0],
                                     reinterpret_cast<fftw_complex*>(in_backward),
                                     reinterpret_cast<fftw_complex*>(out_backward),
                                     FFTW_BACKWARD,FFTW_PATIENT);
    are_FFTarrays_initialized_ = true;
}

Eigen::Vector3i PlaneWaveBasis::get_Gvector(const int &ipw) const
{
    int i1 = ipw%size_FFT_grid_vec_[0];
    int i2 = (ipw/size_FFT_grid_vec_[0])%size_FFT_grid_vec_[1];
    int i3 = ipw/(size_FFT_grid_vec_[0]*size_FFT_grid_vec_[1]);

    i1 = (i1<=size_FFT_grid_vec_[0]/2) ? i1 : i1 - size_FFT_grid_vec_[0];
    i2 = (i2<=size_FFT_grid_vec_[1]/2) ? i2 : i2 - size_FFT_grid_vec_[1];
    i3 = (i3<=size_FFT_grid_vec_[2]/2) ? i3 : i3 - size_FFT_grid_vec_[2];

    return {i1, i2, i3};
}

// convert phi to orbital on FFT grid
// non-collinear calc. not supported
// Note: even in the BAND mode, calc_mode=="SCF" can be used to get the SCF orbitals
void PlaneWaveBasis::get_orbital_FFTgrid(const int &ispin, const int &ik, const int &isym,
                                         const bool is_time_reversal_used_at_k,
                                         const Eigen::VectorXcd &phi,
                                         Eigen::VectorXcd &orbital_FFTgrid,
                                         const std::string &calc_mode) const
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    const std::vector<std::vector<std::vector<Eigen::VectorXi> > > &Gindex_at_k_ref_ =
        calc_mode=="SCF" ? Gindex_at_k_scf_ : Gindex_at_k_band_;
    assert(calc_mode=="SCF" || isym==0); // isym=0 is allowed for calc_mode=="BAND"

    assert(ispin >= 0);
    assert(ispin < Gindex_at_k_ref_.size());

    assert(ik >= 0);
    assert(ik < Gindex_at_k_ref_[ispin].size());

    assert(isym >= 0);
    assert(isym < Gindex_at_k_ref_[ispin][ik].size());

    orbital_FFTgrid = Eigen::VectorXcd::Zero(size_FFT_grid_);
    if (calc_mode=="BAND") // no symmetry is used
    {
        for (int ipw_at_k=0; ipw_at_k<num_G_at_k_band_[ik]; ipw_at_k++)
        {
            orbital_FFTgrid(Gindex_at_k_band_[ispin][ik][0](ipw_at_k)) 
                = phi(ipw_at_k);
        }
    }
    else
    {
        if (is_time_reversal_used_at_k)
        {
            for (int ipw_at_k=0; ipw_at_k<num_G_at_k_scf_[ik]; ipw_at_k++)
            {
                orbital_FFTgrid(Gindex_at_k_scf_[ispin][ik][isym](ipw_at_k)) 
                    = std::conj(phi(ipw_at_k)) // time-reversal sym.
                    * phase_factor_at_k_[ispin][ik][isym](ipw_at_k);
            }
        }
        else
        {
            for (int ipw_at_k=0; ipw_at_k<num_G_at_k_scf_[ik]; ipw_at_k++)
            {
                orbital_FFTgrid(Gindex_at_k_scf_[ispin][ik][isym](ipw_at_k)) 
                    = phi(ipw_at_k)
                    * phase_factor_at_k_[ispin][ik][isym](ipw_at_k);
            }
        }
    }
}

// R-space -> G-space
void PlaneWaveBasis::FFT_forward(const Eigen::VectorXcd &in, Eigen::VectorXcd &out)
{
    assert(in.rows() == size_FFT_grid_ && out.rows() == size_FFT_grid_);
    for (int i=0;i<size_FFT_grid_;++i) { in_forward[i] = in(i); }
    fftw_execute(plan_forward);
    for (int i=0;i<size_FFT_grid_;++i) { out(i) = out_forward[i]/static_cast<double>(size_FFT_grid_); }
}

// G-space -> R-space
void PlaneWaveBasis::FFT_backward(const Eigen::VectorXcd &in, Eigen::VectorXcd &out)
{
    assert(in.rows() == size_FFT_grid_ && out.rows() == size_FFT_grid_);
    for (int i=0;i<size_FFT_grid_;++i) { in_backward[i] = in(i); }
    fftw_execute(plan_backward);
    for (int i=0;i<size_FFT_grid_;++i) { out(i) = out_backward[i]; }
}

void PlaneWaveBasis::bcast(const std::string &calc_mode, const bool am_i_mpi_rank0,
                           const bool setup_FFTgrid)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    std::vector<int> &num_G_at_k_ref_ =
        calc_mode=="SCF" ? num_G_at_k_scf_ : num_G_at_k_band_;
    std::vector<std::vector<std::vector<Eigen::VectorXi> > > &Gindex_at_k_ref_ =
        calc_mode=="SCF" ? Gindex_at_k_scf_ : Gindex_at_k_band_;

    int num_irreducible_kpoints = num_G_at_k_ref_.size();
    MPI_Bcast(&num_irreducible_kpoints, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int num_independent_spins = Gindex_at_k_ref_.size();
    MPI_Bcast(&num_independent_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // set array size of num_G_at_k_ref_ and Gvector_at_k_ref_
    if (!am_i_mpi_rank0) { resize_G_at_k(num_independent_spins, num_irreducible_kpoints, calc_mode); }

    MPI_Bcast(&num_G_at_k_ref_[0], num_irreducible_kpoints, MPI_INT, 0, MPI_COMM_WORLD);

    // bcast Gindex_at_k_ref_ & phase_factor_at_k_ (for SCF)
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            int nsym = Gindex_at_k_ref_[ispin][ik].size();
            MPI_Bcast(&nsym, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (!am_i_mpi_rank0) { Gindex_at_k_ref_[ispin][ik].resize(nsym); }
            if (!am_i_mpi_rank0 && calc_mode=="SCF") { phase_factor_at_k_[ispin][ik].resize(nsym); }

            int npw_at_k = num_G_at_k_ref_[ik];
            for (int isym=0; isym<nsym; isym++)
            {
                if (!am_i_mpi_rank0) { Gindex_at_k_ref_[ispin][ik][isym].resize(npw_at_k); }
                if (!am_i_mpi_rank0 && calc_mode=="SCF") { phase_factor_at_k_[ispin][ik][isym].resize(npw_at_k); }

                MPI_Bcast(Gindex_at_k_ref_[ispin][ik][isym].data(), npw_at_k, MPI_INT, 0, MPI_COMM_WORLD);
                if (calc_mode=="SCF") { MPI_Bcast(phase_factor_at_k_[ispin][ik][isym].data(), npw_at_k, MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD); }
            }
        } // ik
    } // ispin

    if (setup_FFTgrid)
    {
        // bcast FFT grid and setup_FFT()
        std::vector<int> ivec(3);
        if (am_i_mpi_rank0) { ivec = size_FFT_grid_vec_; }
        MPI_Bcast(&ivec[0], 3, MPI_INT, 0, MPI_COMM_WORLD);
        if (!am_i_mpi_rank0) { setup_FFT(ivec[0], ivec[1], ivec[2]); }
    }
}

PlaneWaveBasis::~PlaneWaveBasis()
{
    if (are_FFTarrays_initialized_)
    {
        fftw_destroy_plan(plan_forward);
        fftw_destroy_plan(plan_backward);
        fftw_free(in_forward);
        fftw_free(out_forward);
        fftw_free(in_backward);
        fftw_free(out_backward);
    }
}
