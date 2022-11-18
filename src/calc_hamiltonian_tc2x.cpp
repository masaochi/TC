// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

// two-body, exchange term (TC or BITC)
// calculate -sum_q <*,q| 1/r + \nabla^2 u - (\nabla u)^2 |q,j> (2a_x)
// calculate -sum_q <*,q| \nabla_1 u_12 \nabla_1 |q,j> (2b_x1)
// calculate -sum_q <*,q| \nabla_2 u_21 \nabla_2 |q,j> (2b_x2)
// non-collinear spin not supported
void calc_hamiltonian::tc2x(const Parallelization &parallelization, 
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
    const std::vector<int> num_bands_tc = method.calc_mode()=="SCF" ?
        bloch_states.num_bands_scf() : bloch_states.num_bands_band();

    double two_body_factor = -1.0/(kpoints.num_kpoints()*crystal_structure.unit_cell_volume()); // "-": exchange
    if (num_independent_spins==1) { two_body_factor /= 2.0; } // exchange: parallel spins only, "filling" should be divided by 2
    // coefficient of [2b_x1]: FFT I^2 = -1 is additionally multiplied 
    // coefficient of [2b_x2]: FFT I^2 * convolution (\nabla_2 u_21 = "-1" * \nabla_1 u_12) = +1.

    std::vector<Eigen::Vector3d> Gvect(plane_wave_basis.size_FFT_grid()); // Gvect[ipw](idim) for 3x3 matrix operation
    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
    {
        Gvect[ipw] = crystal_structure.reciprocal_vectors().transpose()
            * plane_wave_basis.get_Gvector(ipw).cast<double>();
    }
    std::vector<Eigen::VectorXd> qGvect(3); // qGvect[idim](ipw)
    for (int idim=0; idim<3; idim++)
    {
        qGvect[idim].resize(plane_wave_basis.size_FFT_grid()); // kqGvect[idim](ipw)
    }
    std::vector<Eigen::VectorXd> kqGvect(3);
    for (int idim=0; idim<3; idim++)
    {
        kqGvect[idim].resize(plane_wave_basis.size_FFT_grid());
    }

    Eigen::Vector3d kqvect, kvect, qvect;
    Eigen::VectorXcd phij(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd phiq(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd chiq;
    if (is_bitc) { chiq.resize(plane_wave_basis.size_FFT_grid()); }
    Eigen::VectorXcd phiqj(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd Hphij_sub(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXd V_2body(plane_wave_basis.size_FFT_grid());

    Eigen::VectorXcd temp(plane_wave_basis.size_FFT_grid());

    std::vector<Eigen::VectorXcd> grad_phij(3);
    for (int idim=0; idim<3; idim++)
    {
        grad_phij[idim].resize(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> grad_phiq(3);
    for (int idim=0; idim<3; idim++)
    {
        grad_phiq[idim].resize(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> phi_grad_phi(3);
    for (int idim=0; idim<3; idim++)
    {
        phi_grad_phi[idim].resize(plane_wave_basis.size_FFT_grid());
    }

    std::vector<Eigen::VectorXd> duk(3);
    for (int idim=0; idim<3; idim++)
    {
        duk[idim].resize(plane_wave_basis.size_FFT_grid());
    }

    // divergence correction for \nabla u
    Eigen::Vector3d kVaux;
    std::vector<Eigen::VectorXcd> div_corr_2bx1;
    std::vector<std::vector<Eigen::VectorXcd> > div_corr_2bx2;

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

            kvect = crystal_structure.reciprocal_vectors().transpose() * kvector_ref;

            // divergence correction for \nabla u
            if (potentials.includes_div_correction()) 
            {
                kVaux = method.calc_mode()=="BAND" ?
                    potentials.sum_of_kVaux_band()[ik] / (kpoints.num_kpoints()*crystal_structure.unit_cell_volume()) :
                    potentials.sum_of_kVaux_scf()[ik] / (kpoints.num_kpoints()*crystal_structure.unit_cell_volume());
                // (-1) of exchange interaction is cancelled with the other "-1":
                //  \int \nabla Vaux - \sum \nabla Vaux = 0 - sum_of_kVaux = "-1" * sum_of_kVaux.
                if (num_independent_spins==1) { kVaux /= 2.0; } // same as two_body_factor
                kVaux *= potentials.jastrow.A_long()[ispin][ispin]; // since we consider not 1/r but A/r-like divergence

                if (kVaux.squaredNorm() > 1e-8) // if = 0 then div_corr_2bx1/x2 is not needed
                {
                    div_corr_2bx1.resize(num_bands_tc[ispin]);
                    div_corr_2bx2.resize(num_bands_tc[ispin]);
                    for (int ibandk=0; ibandk<num_bands_tc[ispin]; ibandk++)
                    {
                        div_corr_2bx2[ibandk].resize(3);
                    }

                    for (int ibandk=0; ibandk<num_bands_tc[ispin]; ibandk++)
                    {
                        double filling_ibandk = method.calc_mode()=="BAND" ?
                            kpoints.return_band_filling(spin, bloch_states.eigenvalues_band()[ispin][ik][ibandk].real(),
                                                        bloch_states.fermi_energy(),
                                                        bloch_states.num_electrons(), ibandk) :
                            kpoints.return_band_filling(spin, bloch_states.eigenvalues_scf()[ispin][ik][ibandk].real(),
                                                        bloch_states.fermi_energy(),
                                                        bloch_states.num_electrons(), ibandk);
                        if (filling_ibandk < 1e-8) { continue; } // no filling                                                                    
                        // NOTE! filling_ibandk is non-zero even for zero-weight k-points (in fake scf) 

                        // [2b_x1]
                        div_corr_2bx1[ibandk] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());

                        const Eigen::VectorXcd &phik_ref = method.calc_mode()=="BAND" ?
                            bloch_states.phik_band()[ispin][ik][ibandk][0] :
                            bloch_states.phik_scf()[ispin][ik][ibandk][0];

                        // Note: names "phiq" and "grad_phiq" are meaningless here.
                        plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, // isym = 0
                                                             ibandk, false, // time-rersal not used for isym=0
                                                             phik_ref, phiq, method.calc_mode());
                        for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                        {
                            for (int idim=0; idim<3; idim++)
                            {
                                grad_phiq[idim](ipw) = (kvect(idim) + Gvect[ipw](idim)) * phiq(ipw);
                            }
                        }
                        Eigen::Vector3d coeff_vec_2bx1 = filling_ibandk * kVaux;
                        for (int idim=0; idim<3; idim++)
                        {
                            div_corr_2bx1[ibandk] += coeff_vec_2bx1(idim) * grad_phiq[idim]; // NOTE! function in k-space.
                        }

                        // [2b_x2]
                        Eigen::VectorXcd &chi_ref_div = is_bitc ? chiq : phiq; // bra orbital
                        if (is_bitc)
                        {
                            const Eigen::VectorXcd &chik_ref = method.calc_mode()=="BAND" ?
                                bloch_states.phik_left_band()[ispin][ik][ibandk][0] :
                                bloch_states.phik_left_scf()[ispin][ik][ibandk][0];

                            plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, ibandk, false,
                                                                 chik_ref, chiq, method.calc_mode());
                        }                        
                        for (int idim=0; idim<3; idim++)
                        {
                            div_corr_2bx2[ibandk][idim] = 
                                kVaux(idim) * chi_ref_div; // NOTE! function in k-space. conjugate() will be later taken.
                        }
                    } // ibandk
                } // if (kVaux != 0 )
            } // if (includes_div_correction)

            for (int jband=0; jband<H2phi[ispin][ik].size(); jband++)
            {
                if (!parallelization.is_assigned_irreducible_kpoints_all_bands()[ispin][ik][jband]) { continue; }
                // calculate phij and grad_phij
                plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, // isym = 0
                                                     jband, false, // time-rersal not used for isym=0
                                                     phi[ispin][ik][jband][0], phij,
                                                     method.calc_mode());
                for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                {
                    for (int idim=0; idim<3; idim++)
                    {
                        grad_phij[idim](ipw) = (kvect(idim) + Gvect[ipw](idim)) * phij(ipw);
                    }
                }
                for (int idim=0; idim<3; idim++)
                {
                    plane_wave_basis.FFT_backward(grad_phij[idim], grad_phij[idim]);
                }
                plane_wave_basis.FFT_backward(phij, phij);
                
                Hphij_sub = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                for (int iq=0; iq<num_irreducible_kpoints_scf; iq++)
                {
                    for (int isymq=0; isymq<kpoints.kvectors_scf()[iq].size(); isymq++)
                    {
                        if (bloch_states.num_occupied_bands()[ispin][iq]==0) { continue; }

                        qvect = crystal_structure.reciprocal_vectors().transpose() *
                            kpoints.kvectors_scf()[iq][isymq];
                        kqvect = kvect - qvect;
                        for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                        {
                            for (int idim=0; idim<3; idim++)
                            {
                                qGvect[idim](ipw) = qvect(idim) + Gvect[ipw](idim); // NOTE: order of (ipw,idim)
                                kqGvect[idim](ipw) = kqvect(idim) + Gvect[ipw](idim);
                            }
                        }

                        // calculate Jastrow
                        for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                        {
                            V_2body(ipw)
                                = potentials.jastrow.tc_2body(kqvect + Gvect[ipw], ispin, ispin);
                        
                            double uk = potentials.jastrow.uk(kqvect + Gvect[ipw], ispin, ispin);
                            for (int idim=0; idim<3; idim++)
                            {
                                duk[idim](ipw) = kqGvect[idim](ipw) * uk;
                            }
                        }
                        for (int ibandq=0; ibandq<bloch_states.num_occupied_bands()[ispin][iq]; ibandq++)
                        {
                            // calculate phiq and grad_phiq
                            Eigen::VectorXcd &chiq_ref = is_bitc ? chiq : phiq; // bra orbital

                            plane_wave_basis.get_orbital_FFTgrid(ispin, iq, isymq, ibandq,
                                                                 kpoints.is_time_reversal_used_at_k()[iq][isymq],
                                                                 bloch_states.phik_scf()[ispin][iq][ibandq][0], phiq,
                                                                 "SCF");
                            if (is_bitc)
                            {
                                plane_wave_basis.get_orbital_FFTgrid(ispin, iq, isymq, ibandq,
                                                                     kpoints.is_time_reversal_used_at_k()[iq][isymq],
                                                                     bloch_states.phik_left_scf()[ispin][iq][ibandq][0], chiq,
                                                                     "SCF");
                                plane_wave_basis.FFT_backward(chiq, chiq);
                            }

                            for (int idim=0; idim<3; idim++)
                            {
                                grad_phiq[idim] = qGvect[idim].array() * phiq.array();
                                plane_wave_basis.FFT_backward(grad_phiq[idim], grad_phiq[idim]);
                            }
                            plane_wave_basis.FFT_backward(phiq, phiq);

                            // [2a_x & 2b_x1]
                            phiqj = chiq_ref.conjugate().array() * phij.array(); // in R-space
                            plane_wave_basis.FFT_forward(phiqj, phiqj); // -> phiqj(G)

                            // [2b_x1] <*,q| \nabla_1 u_12 |*,j>
                           for (int idim=0; idim<3; idim++)
                            {
                                phi_grad_phi[idim] = phiqj.array() * duk[idim].array(); // note: meaningless name "phi_grad_phi"
                                plane_wave_basis.FFT_backward(phi_grad_phi[idim], phi_grad_phi[idim]);
                            }

                           // [2b_x1] <*,q| \nabla_1 u_12 \nabla_1 |q,j>
                           temp = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                           for (int idim=0; idim<3; idim++)
                            {
                                temp = temp.array()
                                    + phi_grad_phi[idim].array() * grad_phiq[idim].array();
                            }
                           Hphij_sub -= bloch_states.filling()[ispin][iq][ibandq] * temp; // -1 = I^2

                            // [2b_x2]
                            for (int idim=0; idim<3; idim++)
                            {
                                phi_grad_phi[idim] = chiq_ref.conjugate().array() * grad_phij[idim].array();
                                plane_wave_basis.FFT_forward(phi_grad_phi[idim], phi_grad_phi[idim]);
                            }

                            // [2a_x] <*,q| V_2body |j,*>
                            phiqj = phiqj.array() * V_2body.array(); // in G-space
                            // [2b_x2] <*,q| \nablba_2 u_12 \nabla_2 |*,j>
                            for (int idim=0; idim<3; idim++) 
                            {
                                phiqj = phiqj.array() +
                                    phi_grad_phi[idim].array() * duk[idim].array(); // in G-space
                            }
                            plane_wave_basis.FFT_backward(phiqj, phiqj); // -> (R)
                            Hphij_sub = Hphij_sub.array() +
                                bloch_states.filling()[ispin][iq][ibandq] * phiqj.array() * phiq.array(); // in R-space
                        } // ibandq
                    } // isymq
                } // iq
                plane_wave_basis.FFT_forward(Hphij_sub, Hphij_sub);
                for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                {
                    H2phi[ispin][ik][jband][0](ipw_at_k) += two_body_factor * Hphij_sub(Gindex_at_k(ipw_at_k));
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
                    
                    for (int ibandk=0; ibandk<num_bands_tc[ispin]; ibandk++)
                    {
                        double filling_ibandk = method.calc_mode()=="BAND" ?
                            kpoints.return_band_filling(spin, bloch_states.eigenvalues_band()[ispin][ik][ibandk].real(),
                                                        bloch_states.fermi_energy(),
                                                        bloch_states.num_electrons(), ibandk) :
                            kpoints.return_band_filling(spin, bloch_states.eigenvalues_scf()[ispin][ik][ibandk].real(),
                                                        bloch_states.fermi_energy(),
                                                        bloch_states.num_electrons(), ibandk);
                        if (filling_ibandk < 1e-8) { continue; } // no filling then no correction
                        
                        const Eigen::VectorXcd &phik_ref = method.calc_mode()=="BAND" ?
                            bloch_states.phik_band()[ispin][ik][ibandk][0] :
                            bloch_states.phik_scf()[ispin][ik][ibandk][0];
                        
                        const Eigen::VectorXcd &chik_ref = !is_bitc ? phik_ref :
                            (method.calc_mode()=="BAND" ? bloch_states.phik_left_band()[ispin][ik][ibandk][0] :
                             bloch_states.phik_left_scf()[ispin][ik][ibandk][0]);
                        
                        Complex prod = chik_ref.dot(phi[ispin][ik][jband][0]); // conj(phii)*phij
                        prod *= div_corr; // [2a_x]
                        
                        if (kVaux.squaredNorm() > 1e-8) // \sum \nabla u correction (e.g. not required for non-zero k-weight in SCF)
                        {
                            for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                            {
                                // NOTE! function in k-space
                                H2phi[ispin][ik][jband][0](ipw_at_k) 
                                    -= div_corr_2bx1[ibandk](Gindex_at_k(ipw_at_k)); // [2b_x1] I^2 = -1
                            }
                            for (int idim=0; idim<3; idim++)
                            {
                                // NOTE! inner product in k-space
                                plane_wave_basis.FFT_forward(grad_phij[idim], grad_phiq[idim]); // save in "grad_phiq" (meaningless name)
                                prod += div_corr_2bx2[ibandk][idim].dot(grad_phiq[idim]); // [2b_x2] I^2 * convolution = +1. conj(div)*gphij
                            }
                        } // if (||kVaux||^2 > 1e-8)
                        
                        prod *= filling_ibandk;
                        H2phi[ispin][ik][jband][0] += prod * phik_ref;
                    } // ibandk
                } // divergence correction
            } // jband
        } // ik
    } // ispin
}

