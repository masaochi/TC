// [class BlochStates]
// wave functions, eigenvalues, filling, num. of bands, etc.

#include "include/header.hpp"

// private

// called in set_qe_xml()
void BlochStates::set_num_bands_qe(const std::vector<int> &num_bands_qe,
                                   const std::string &calc_mode) 
{ 
    const int num_independent_spins = num_bands_qe.size();

    assert(calc_mode=="SCF" || calc_mode=="BAND");
    std::vector<int> &num_bands_tc_ 
        = calc_mode=="SCF" ? num_bands_scf_ : num_bands_band_;

    // set num_bands_qe_
    num_bands_qe_.resize(num_independent_spins);
    for (auto nb : num_bands_qe)
    {
        assert(nb > 0);
    }
    num_bands_qe_ = num_bands_qe;

    // set num_bands_tc_
    if (num_bands_tc_[0]<0) // num_bands_tc_ is not specified in input.in
    {
        num_bands_tc_ = num_bands_qe_; // default
    }
    else // check consistency
    {
        const int nbtmp = num_bands_tc_[0];
        for (int i=0; i<num_independent_spins; i++)
        {
            if (nbtmp > num_bands_qe_[i]) 
            {
                error_messages::inappropriate_argument("num_bands_tc", nbtmp,
                                                       "should satisfy num_bands_tc <= num_bands_qe");
            }
        }
        num_bands_tc_.assign(num_independent_spins, nbtmp); // change the size of the array and filled with num_bands_tc_[0]
    }
}

// called in set_qe_xml(). also set eigenvalues_scf/band = eigenvalues_qe
void BlochStates::set_eigenvalues_qe(const std::vector<std::vector<double> > &eigenvalues,
                                     const std::string &calc_mode)
{
    const int num_irreducible_kpoints = eigenvalues.size();
    const int num_independent_spins = num_bands_qe_.size();
    const int num_bands_qe_total = num_independent_spins==2 ? num_bands_qe_[0]+num_bands_qe_[1] : num_bands_qe_[0];

    assert(calc_mode=="SCF" || calc_mode=="BAND");
    std::vector<std::vector<std::vector<Complex> > > &eigenvalues_tc_
        = calc_mode=="SCF" ? eigenvalues_scf_ : eigenvalues_band_;
    const std::vector<int> &num_bands_tc_ 
        = calc_mode=="SCF" ? num_bands_scf_ : num_bands_band_;

    for (int ik=0; ik<num_irreducible_kpoints; ik++)
    {
        assert(eigenvalues[ik].size() == num_bands_qe_total);
    }

    eigenvalues_qe_.resize(num_independent_spins);
    eigenvalues_tc_.resize(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        eigenvalues_qe_[ispin].resize(num_irreducible_kpoints);
        eigenvalues_tc_[ispin].resize(num_irreducible_kpoints);
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            eigenvalues_qe_[ispin][ik].resize(num_bands_qe_[ispin]);
            eigenvalues_tc_[ispin][ik].resize(num_bands_tc_[ispin]); // Note: num_band_tc_ (not num_bands_qe_)

            int i = ispin==0 ? 0 : num_bands_qe_[0];
            for (int iband=0; iband<num_bands_qe_[ispin]; iband++)
            {
                eigenvalues_qe_[ispin][ik][iband] = eigenvalues[ik][i + iband];
                if (iband<num_bands_tc_[ispin]) { eigenvalues_tc_[ispin][ik][iband] = eigenvalues_qe_[ispin][ik][iband]; }
            }
        }
    }
}

void BlochStates::resize_filling(const int &num_independent_spins,
                                 const int &num_irreducible_kpoints_scf)
{
    assert(num_independent_spins > 0);
    assert(num_irreducible_kpoints_scf > 0);

    filling_.resize(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        filling_[ispin].resize(num_irreducible_kpoints_scf);
        for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
        {
            filling_[ispin][ik].resize(num_bands_scf_[ispin]);
        }
    }

    // initialize (while not used in BAND calculation...)
    filling_old_.resize(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        filling_old_[ispin].resize(num_irreducible_kpoints_scf);
        for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
        {
            filling_old_[ispin][ik].resize(0); // resized in the SCF loop. Only occupied bands are allocated for filling_old_.
        }
    }
}

void BlochStates::resize_num_occupied_bands(const int &num_independent_spins,
                                            const int &num_irreducible_kpoints_scf)
{
    assert(num_independent_spins > 0);
    assert(num_irreducible_kpoints_scf > 0);

    num_occupied_bands_.resize(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        num_occupied_bands_[ispin].resize(num_irreducible_kpoints_scf);
    }
}

// public

// In the SCF calculation, this function is called once when reading input.in with the "SCF" switch.
// In the BAND calculation, this function is called twice.
//   First call = "BAND" switch when reading input.in to initialize num_bands_band
//   Second call = "SCF" switch when reading tc_scfinfo.dat to initialize num_bands_scf
// For both cases, the first call is skipped when num_bands_tc is not specified in input.in.
// In that case, num_bands_tc = num_bands_qe is used as a default setting in set_num_bands_qe().
void BlochStates::set_num_bands_tc(const std::vector<int> &num_bands_tc,
                                   const std::string &calc_mode) 
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    std::vector<int> &num_bands_tc_ 
        = calc_mode=="SCF" ? num_bands_scf_ : num_bands_band_;

    for (auto nb : num_bands_tc)
    {
        if (nb <= 0) { error_messages::inappropriate_argument("num_bands_tc", nb, "should be positive"); }
    }
    num_bands_tc_ = num_bands_tc;
}

void BlochStates::set_qe_xml(const std::vector<int> &num_bands_qe,
                             const std::vector<std::vector<double> > &eigenvalues,
                             const std::string &calc_mode,
                             const double &num_electrons)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");

    set_num_bands_qe(num_bands_qe, calc_mode);
    set_eigenvalues_qe(eigenvalues, calc_mode); // calc_mode==SCF or BAND
    assert(num_electrons > 0);
    num_electrons_ = num_electrons;
}

// phik[spin][k-point][band][num_spinor][npw_at_k]: wave function in G-space at each irreducible k-point
void BlochStates::resize_phik_qe(const PlaneWaveBasis &plane_wave_basis,
                                 const bool is_spinor, const int &num_independent_spins,
                                 const std::string &calc_mode)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    assert(num_independent_spins > 0);

    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_ref_
        = calc_mode=="SCF" ? phik_scf_ : phik_band_;
    const std::vector<int> &num_G_at_k = calc_mode=="SCF" ?
        plane_wave_basis.num_G_at_k_scf() : plane_wave_basis.num_G_at_k_band();
    const std::vector<int> &num_bands_tc_ 
        = calc_mode=="SCF" ? num_bands_scf_ : num_bands_band_;

    const int num_irreducible_kpoints = num_G_at_k.size();
    assert(num_irreducible_kpoints > 0);
    const int num_spinor = (!is_spinor) ? 1 : 2;

    phik_ref_.resize(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        phik_ref_[ispin].resize(num_irreducible_kpoints);
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            assert(num_bands_tc_[ispin] > 0);
            phik_ref_[ispin][ik].resize(num_bands_tc_[ispin]);
            for (int iband=0; iband<num_bands_tc_[ispin]; iband++)
            {
                phik_ref_[ispin][ik][iband].resize(num_spinor);
                for (int ispinor=0; ispinor<num_spinor; ispinor++)
                {
                    assert(num_G_at_k[ik] > 0);
                    phik_ref_[ispin][ik][iband][ispinor].resize(num_G_at_k[ik]);
                    // phik_ref_[ispin][ik][iband][ispinor].shrink_to_fit(); // memory saving
                }
            }
        }
    }

    if (calc_mode=="SCF")
    {
        phik_scf_old_.resize(num_independent_spins);
        for (int ispin=0; ispin<num_independent_spins; ispin++)
        {
            phik_scf_old_[ispin].resize(num_irreducible_kpoints);
            for (int ik=0; ik<num_irreducible_kpoints; ik++)
            {
                phik_scf_old_[ispin][ik].resize(0); // resized in the SCF loop
            }
        }
    }
}

// see also plane_wave_basis.set_Gvector_at_k() for gamma_only treatment
void BlochStates::set_phik_qe(const int ispin, const int ik, const std::vector<Complex> &evc,
                              const std::string &calc_mode, const bool gamma_only, const int zero_index)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    assert(!gamma_only || zero_index>=0); // zero_index should be set (>=0) for gamma_only==true
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_ref_
        = calc_mode=="SCF" ? phik_scf_ : phik_band_;
    const std::vector<int> &num_bands_tc_ 
        = calc_mode=="SCF" ? num_bands_scf_ : num_bands_band_;

    assert(ispin >= 0);
    assert(ispin < phik_ref_.size());

    assert(ik >= 0);
    assert(ik < phik_ref_[ispin].size());

    int icount = 0;
    for (int iband=0; iband<num_bands_tc_[ispin]; iband++)
    {
        int num_spinor = phik_ref_[ispin][ik][iband].size();
        double norm = 0.0;
        for (int ispinor=0; ispinor<num_spinor; ispinor++)
        {
            int npw_at_k_org = phik_ref_[ispin][ik][iband][ispinor].size();
            int npw_at_k;
            if (gamma_only) // u(-G)=u(G)^*
            {
                npw_at_k = 2*npw_at_k_org - 1;
                phik_ref_[ispin][ik][iband][ispinor].resize(npw_at_k);
            }
            else
            {
                npw_at_k = npw_at_k_org;
            }
            for (int ipw_at_k=0; ipw_at_k<npw_at_k_org; ipw_at_k++)
            {
                phik_ref_[ispin][ik][iband][ispinor](ipw_at_k) = evc[icount]; // evc(npol*igwx, nbnd[is]) in fortran
                norm += std::norm(evc[icount]);
                if (gamma_only && ipw_at_k!=zero_index) // gamma_only && G!=(0,0,0)
                {
                    int ipw_at_k_minus = ipw_at_k<zero_index ? ipw_at_k + npw_at_k_org : ipw_at_k + npw_at_k_org - 1; // -G index
                    phik_ref_[ispin][ik][iband][ispinor](ipw_at_k_minus) = std::conj(evc[icount]); // u(-G)=u(G)^*
                    norm += std::norm(evc[icount]);
                }
                ++icount;
            }
        }
        assert(std::abs(norm - 1.0) < 1e-8); // check normalization    
    }
}

void BlochStates::bcast_qe_input(const PlaneWaveBasis &plane_wave_basis,
                                 const bool is_spinor,
                                 const std::string &calc_mode, const std::string &calc_method,
                                 const bool am_i_mpi_rank0)
{
    int num_irreducible_kpoints;
    if (am_i_mpi_rank0) { num_irreducible_kpoints = eigenvalues_qe_[0].size(); }
    MPI_Bcast(&num_irreducible_kpoints, 1, MPI_INT, 0, MPI_COMM_WORLD);

    assert(calc_mode=="SCF" || calc_mode=="BAND");
    std::vector<std::vector<std::vector<Complex> > > &eigenvalues_tc_
        = calc_mode=="SCF" ? eigenvalues_scf_ : eigenvalues_band_;
    std::vector<int> &num_bands_tc_ 
        = calc_mode=="SCF" ? num_bands_scf_ : num_bands_band_;

    // num_bands_tc & num_bands_qe ... the former can be specified by input.in or initialized by io_qe_files::read_qe().
    // Thus, Bcast is performed here.
    int num_independent_spins = num_bands_qe_.size();
    MPI_Bcast(&num_independent_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> ivec(num_independent_spins);
    if (am_i_mpi_rank0) { ivec = num_bands_tc_; }
    MPI_Bcast(&ivec[0], ivec.size(), MPI_INT, 0, MPI_COMM_WORLD);
    if (!am_i_mpi_rank0) { set_num_bands_tc(ivec, calc_mode); }
    if (am_i_mpi_rank0) { ivec = num_bands_qe_; }
    MPI_Bcast(&ivec[0], ivec.size(), MPI_INT, 0, MPI_COMM_WORLD);
    if (!am_i_mpi_rank0) { set_num_bands_qe(ivec, calc_mode); }

    // eigenvalues_qe_ and eigenvalues_tc_
    if (!am_i_mpi_rank0) { eigenvalues_qe_.resize(num_independent_spins); }
    if (!am_i_mpi_rank0) { eigenvalues_tc_.resize(num_independent_spins); }
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        if (!am_i_mpi_rank0) { eigenvalues_qe_[ispin].resize(num_irreducible_kpoints); }
        if (!am_i_mpi_rank0) { eigenvalues_tc_[ispin].resize(num_irreducible_kpoints); }
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            if (!am_i_mpi_rank0) { eigenvalues_qe_[ispin][ik].resize(num_bands_qe_[ispin]); }
            if (!am_i_mpi_rank0) { eigenvalues_tc_[ispin][ik].resize(num_bands_tc_[ispin]); }
            MPI_Bcast(&eigenvalues_qe_[ispin][ik][0], eigenvalues_qe_[ispin][ik].size(),
                      MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&eigenvalues_tc_[ispin][ik][0], eigenvalues_tc_[ispin][ik].size(),
                      MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
        }
    }

    // num_electrons_
    MPI_Bcast(&num_electrons_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // phik
    if (!am_i_mpi_rank0) { resize_phik_qe(plane_wave_basis, is_spinor, num_independent_spins, calc_mode); }
    bcast_phik(false, false, calc_mode, am_i_mpi_rank0); // false(1st): do not bcast phik_left. false(2nd): do not bcast phik_scf_old
    if (calc_method=="BITC" && calc_mode=="SCF") // phik_left
    {
        phik_left_scf_ = phik_scf_;
        phik_left_scf_old_ = phik_scf_old_;
    }
    if (calc_method=="BITC" && calc_mode=="BAND") { phik_left_band_ = phik_band_; } // phik_left

    // not bcast but required initialization
    // In BAND calculation, resize_filling & resize_num_occupied_bands will be called in bcast_scfinfo().
    if (calc_mode=="SCF")
    {
        resize_filling(num_independent_spins, num_irreducible_kpoints);
        resize_num_occupied_bands(num_independent_spins, num_irreducible_kpoints);
    }
}

// called in BAND calculation
void BlochStates::bcast_scfinfo(const PlaneWaveBasis &plane_wave_basis,
                                const bool is_spinor,
                                const std::string &calc_method,
                                const bool am_i_mpi_rank0)
{
    int num_irreducible_kpoints_scf;
    if (am_i_mpi_rank0) { num_irreducible_kpoints_scf = plane_wave_basis.num_G_at_k_scf().size(); }
    MPI_Bcast(&num_irreducible_kpoints_scf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int num_independent_spins;
    if (am_i_mpi_rank0) { num_independent_spins = num_bands_scf_.size(); }
    MPI_Bcast(&num_independent_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // num_bands_scf_
    if (!am_i_mpi_rank0) { num_bands_scf_.resize(num_independent_spins); }
    MPI_Bcast(&num_bands_scf_[0], num_independent_spins, MPI_INT, 0, MPI_COMM_WORLD);

    // allocate eigenvalues_scf_
    eigenvalues_scf_.resize(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        eigenvalues_scf_[ispin].resize(num_irreducible_kpoints_scf);
        for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
        {
            eigenvalues_scf_[ispin][ik].resize(num_bands_scf_[ispin]);
        }
    }

    // allocate phik_scf (& phik_left_scf)
    resize_phik_qe(plane_wave_basis, is_spinor, num_independent_spins, "SCF");
    if (calc_method=="BITC") { phik_left_scf_ = phik_scf_; } // phik_left

    // allocate the following variables
    resize_filling(num_independent_spins, num_irreducible_kpoints_scf);
    resize_num_occupied_bands(num_independent_spins, num_irreducible_kpoints_scf);
}

// set phik_scf by reading tc_wfc.dat
// non-collinear calculation not supported
void BlochStates::set_phik_from_tc_wfc(const int &ispin, const int &ik, const int &iband,
                                       const Eigen::VectorXcd &wfc_temp,
                                       const std::string &right_or_left,
                                       const std::string &calc_mode)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    assert(right_or_left=="right" || right_or_left=="left");
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_ref_ = 
        calc_mode=="SCF" ?
        (right_or_left=="right" ? phik_scf_ : phik_left_scf_) :
        (right_or_left=="right" ? phik_band_ : phik_left_band_);

    assert(ispin >= 0 && ispin < phik_ref_.size());
    assert(ik >= 0 && ik < phik_ref_[ispin].size());
    assert(iband >=0 && iband < phik_ref_[ispin][0].size());

    const int num_spinor = phik_ref_[ispin][ik][iband].size();
    assert(num_spinor==1);

    const int num_G_at_k = phik_ref_[ispin][ik][iband][0].size();
    assert(wfc_temp.size() == num_G_at_k);

    phik_ref_[ispin][ik][iband][0] = wfc_temp;
}

// set eigenvalues_scf by reading tc_energy.dat
void BlochStates::set_eigenvalues_from_tc_energy(const int &ispin, const int &ik,
                                                 const Eigen::VectorXcd &eigenvalues_temp,
                                                 const std::string &calc_mode)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    std::vector<std::vector<std::vector<Complex> > > &eigenvalues_ref_ = calc_mode=="SCF" ?
        eigenvalues_scf_ : eigenvalues_band_;
    const std::vector<int> num_bands_ref_ = calc_mode=="SCF" ? 
        num_bands_scf_ : num_bands_band_;

    assert(ispin >= 0 && ispin < eigenvalues_ref_.size());
    assert(ik >= 0 && ik < eigenvalues_ref_[ispin].size());
    assert(eigenvalues_temp.size() == num_bands_ref_[ispin]);
    assert(eigenvalues_temp.size() == eigenvalues_ref_[ispin][ik].size());

    for (int iband=0; iband<num_bands_ref_[ispin]; iband++)
    {
        eigenvalues_ref_[ispin][ik][iband] = eigenvalues_temp(iband);
    }
}
