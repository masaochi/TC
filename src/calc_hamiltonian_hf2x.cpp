// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

// two-body, exchange term (Hartree-Fock)
// non-collinear spin not supported
void calc_hamiltonian::hf2x(const Parallelization &parallelization, 
                            const Method &method,
                            const CrystalStructure &crystal_structure,
                            const Potentials &potentials,
                            const Spin &spin, const Kpoints &kpoints,
                            PlaneWaveBasis &plane_wave_basis, 
                            const BlochStates &bloch_states, 
                            const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
                            std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H2phi,
                            std::ostream *ost)
{
    assert(!spin.is_spinor());
    const bool is_bitc = method.calc_method()=="BITC" ? true : false;
    const int num_independent_spins = spin.num_independent_spins();
    const int num_irreducible_kpoints_scf = kpoints.num_irreducible_kpoints_scf();
    const int num_irreducible_kpoints_ref = method.calc_mode()=="SCF" ? // kpoints for "phi"
        kpoints.num_irreducible_kpoints_scf() : kpoints.num_irreducible_kpoints_band();
    const int num_kpoints_all_scf = kpoints.num_kpoints_all_scf();
    const std::vector<int> num_bands_tc = method.calc_mode()=="SCF" ?
        bloch_states.num_bands_scf() : bloch_states.num_bands_band();

    double two_body_factor = -FourPI/(kpoints.num_kpoints()*crystal_structure.unit_cell_volume()); // "-": exchange, 4pi: Coulomb
    if (num_independent_spins==1) { two_body_factor /= 2.0; } // exchange: parallel spins only, "filling" should be divided by 2

    std::vector<Eigen::Vector3d> Gvect(plane_wave_basis.size_FFT_grid());
    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
    {
        Gvect[ipw] = crystal_structure.reciprocal_vectors().transpose()
            * plane_wave_basis.get_Gvector(ipw).cast<double>();
    }

    Eigen::Vector3d kqvect;
    Eigen::VectorXcd phij(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd phiq(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd phiqj(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd Hphij_sub(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXd V_coulomb(plane_wave_basis.size_FFT_grid());

    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<num_irreducible_kpoints_ref; ik++)
        {
            const Eigen::Vector3d kvector_ref = method.calc_mode()=="SCF" ? 
                kpoints.kvectors_scf()[ik][0] : kpoints.kvectors_band()[ik][0];
            const int num_G_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.num_G_at_k_scf()[ik] : plane_wave_basis.num_G_at_k_band()[ik];
            const Eigen::VectorXi Gindex_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.Gindex_at_k_scf()[ispin][ik][0] :  plane_wave_basis.Gindex_at_k_band()[ispin][ik][0];

            for (int jband=0; jband<H2phi[ispin][ik].size(); jband++)
            {
                if (!parallelization.is_assigned_irreducible_kpoints_all_bands()[ispin][ik][jband]) { continue; }
                plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, // isym = 0
                                                     false, // time-rersal not used for isym=0
                                                     phi[ispin][ik][jband][0], phij,
                                                     method.calc_mode());
                plane_wave_basis.FFT_backward(phij, phij);
  
                Hphij_sub = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
#ifdef _OPENMP
                #pragma omp parallel firstprivate(kqvect, V_coulomb, phiq, phiqj)
                {
                    PlaneWaveBasis plane_wave_basis_thread;
                    #pragma omp critical // making a plan in FFTW is not thread-safe.
                    {
                        plane_wave_basis_thread = plane_wave_basis; 
                    }
                    
                    #pragma omp for
#else
                    PlaneWaveBasis& plane_wave_basis_thread = plane_wave_basis;
#endif
                    for (int iq_isymq=0; iq_isymq<num_kpoints_all_scf; iq_isymq++)
                    {
                        // We use a "iq_isymq" loop instead of "iq" & "isymq" double loops
                        // for OpenMP parallelization with better efficiency. 
                        int iq = kpoints.index_all_kscf()[iq_isymq][0];
                        int isymq = kpoints.index_all_kscf()[iq_isymq][1];

                        assert(iq>=0 && iq<num_irreducible_kpoints_scf);
                        assert(isymq>=0 && isymq<kpoints.kvectors_scf()[iq].size());

                        if (bloch_states.num_occupied_bands()[ispin][iq]==0) { continue; }
                        kqvect = crystal_structure.reciprocal_vectors().transpose()
                            * (kvector_ref - kpoints.kvectors_scf()[iq][isymq]);
                        for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                        {
                            double tmp = (kqvect + Gvect[ipw]).squaredNorm();
                            if (tmp<1e-8) 
                            {
                                V_coulomb(ipw) = 0.0;
                            }
                            else
                            {
                                V_coulomb(ipw) = 1.0/tmp;
                            }
                        } // ipw
                        
                        const int nbands_old = bloch_states.filling_old()[ispin][iq].size();
                        for (int ibandq=-nbands_old; ibandq<bloch_states.num_occupied_bands()[ispin][iq]; ibandq++)
                        {
                            if (ibandq>=0)
                            {
                                plane_wave_basis.get_orbital_FFTgrid(ispin, iq, isymq,
                                                                     kpoints.is_time_reversal_used_at_k()[iq][isymq],
                                                                     bloch_states.phik_scf()[ispin][iq][ibandq][0], phiq,
                                                                     "SCF");
                            }
                            else
                            {
                                plane_wave_basis.get_orbital_FFTgrid(ispin, iq, isymq,
                                                                     kpoints.is_time_reversal_used_at_k()[iq][isymq],
                                                                     bloch_states.phik_scf_old()[ispin][iq][-1-ibandq][0], phiq,
                                                                     "SCF");
                            }
                            plane_wave_basis_thread.FFT_backward(phiq, phiq);
                            phiqj = phiq.conjugate().array() * phij.array(); // in R-space
                            plane_wave_basis_thread.FFT_forward(phiqj, phiqj);
                            phiqj = phiqj.array() * V_coulomb.array(); // in G-space
                            plane_wave_basis_thread.FFT_backward(phiqj, phiqj);
                            
                            double filq = ibandq>=0 ? bloch_states.filling()[ispin][iq][ibandq] : bloch_states.filling_old()[ispin][iq][-1-ibandq];
#ifdef _OPENMP
                            #pragma omp critical
#endif
                            {
                                Hphij_sub = Hphij_sub.array() + filq * phiqj.array() * phiq.array(); // in R-space
                            }
                        } // iband
                    } // iq_isymq
                
#ifdef _OPENMP
                } // pragma omp parallel
#endif

                plane_wave_basis.FFT_forward(Hphij_sub, Hphij_sub);
                for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                {
                    H2phi[ispin][ik][jband][0](ipw_at_k) += two_body_factor*Hphij_sub(Gindex_at_k(ipw_at_k));
                }

                // Divergence correction
                if (potentials.includes_div_correction())
                {
                    // \int d^3q/(2pi)^3 4pi*exp(-alpha q^2)/q^2 = 2/pi \int dq_0^infty exp(-alpha q^2)/q^2 = 1/sqrt(pi*alpha)
                    double div_corr = 1.0/std::sqrt(PI*potentials.alpha_Vaux());
                    if (method.calc_mode()=="SCF") 
                    {
                        div_corr -= potentials.sum_of_Vaux_scf()[ik]/(kpoints.num_kpoints()*crystal_structure.unit_cell_volume());
                    }
                    else if (method.calc_mode()=="BAND") 
                    {
                        div_corr -= potentials.sum_of_Vaux_band()[ik]/(kpoints.num_kpoints()*crystal_structure.unit_cell_volume());
                    }
                    div_corr *= (-1); // "-": exchange (same as two_body_factor)
                    if (num_independent_spins==1) { div_corr /= 2.0; } // same as two_body_factor
                    
                    const int nbands_old = 
                        method.calc_mode()=="BAND" ? 0 : bloch_states.filling_old()[ispin][ik].size();
                    for (int ibandk=-nbands_old; ibandk<num_bands_tc[ispin]; ibandk++)
                    {
                        double filling_ibandk = 
                            ibandk>=0 ? (method.calc_mode()=="BAND" ?
                                         kpoints.return_band_filling(spin, bloch_states.eigenvalues_band()[ispin][ik][ibandk].real(),
                                                                     bloch_states.fermi_energy(),
                                                                     bloch_states.num_electrons(), ibandk) :
                                         bloch_states.filling()[ispin][ik][ibandk]) :
                            bloch_states.filling_old()[ispin][ik][-1-ibandk];

                        if (ibandk>=0 && method.calc_mode()=="SCF")
                        {
                            // above filling_ibandk (=filling[ispin][ik][ibandk]) is not correct for zero-weight k-points in SCF calc.
                            if (kpoints.kweight_scf()[ik] < 1e-8)
                            {
                                filling_ibandk = 
                                    kpoints.return_band_filling(spin, bloch_states.eigenvalues_scf()[ispin][ik][ibandk].real(),
                                                                bloch_states.fermi_energy(),
                                                                bloch_states.num_electrons(), ibandk);
                            }
                        }

                        if (filling_ibandk < 1e-8) { continue; } // no filling then no correction
                            
                        const Eigen::VectorXcd &phik_ref =
                            ibandk>=0 ? (method.calc_mode()=="BAND" ?
                                         bloch_states.phik_band()[ispin][ik][ibandk][0] :
                                         bloch_states.phik_scf()[ispin][ik][ibandk][0]) :
                            bloch_states.phik_scf_old()[ispin][ik][-1-ibandk][0];
 
                        Complex prod = phik_ref.dot(phi[ispin][ik][jband][0]); // conj(phik)*phi
                        prod *= filling_ibandk * div_corr;
                        H2phi[ispin][ik][jband][0] += prod * phik_ref;
                    } // ibandk
                } // if (includes_div_correction)
            } // jband
        } // ik
    } // ispin
}

