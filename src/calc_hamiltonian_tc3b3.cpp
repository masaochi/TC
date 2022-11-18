// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

// *** tc3b3 ***
// calculate 0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |q1,j,q2>
void calc_hamiltonian::tc3b3(const Parallelization &parallelization, 
                             const Method &method,
                             const CrystalStructure &crystal_structure,
                             const Potentials &potentials,
                             const Spin &spin, const Kpoints &kpoints,
                             PlaneWaveBasis &plane_wave_basis, 
                             const BlochStates &bloch_states, 
                             const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
                             std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H3phi,
                             std::ostream *ost)
{
    assert(!spin.is_spinor());
    const bool is_bitc = method.calc_method()=="BITC" ? true : false;
    const int num_independent_spins = spin.num_independent_spins();
    const int num_irreducible_kpoints_scf = kpoints.num_irreducible_kpoints_scf();
    const int num_irreducible_kpoints = method.calc_mode()=="SCF" ? 
        kpoints.num_irreducible_kpoints_scf() : kpoints.num_irreducible_kpoints_band();
    const int num_spinor = (spin.is_spinor()==false) ? 1 : 2;
    const std::vector<int> num_bands_tc = method.calc_mode()=="SCF" ?
        bloch_states.num_bands_scf() : bloch_states.num_bands_band();

    // (0.5*2) * I^2 * (-1) (convolution \nabla_2 u_21 = -\nabla_1 u_12) = 1.0
    double three_body_factor = 1.0/(kpoints.num_kpoints()*kpoints.num_kpoints()*
                                    crystal_structure.unit_cell_volume()*
                                    crystal_structure.unit_cell_volume());
    // spin indices of x1 and x2 should be the same: filling should be divided with 2
    if (num_independent_spins==1) { three_body_factor /= 2.0; } 

    std::vector<Eigen::Vector3d> Gvect(plane_wave_basis.size_FFT_grid());
    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
    {
        Gvect[ipw] = crystal_structure.reciprocal_vectors().transpose()
            * plane_wave_basis.get_Gvector(ipw).cast<double>();
    }

    Eigen::Vector3d kqvect, kvect, qvect;
    Eigen::VectorXcd orbital(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd phij(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd phiq(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd chiq;
    if (is_bitc) { chiq.resize(plane_wave_basis.size_FFT_grid()); }
    Eigen::VectorXcd Hphij_sub(plane_wave_basis.size_FFT_grid());

    std::vector<std::vector<Eigen::VectorXcd> > dnu(num_independent_spins); // spin index of x2
    for (int ispin=0; ispin<num_independent_spins; ispin++) 
    {
        dnu[ispin].resize(3);
        for (int idim=0; idim<3; idim++)
        {
            dnu[ispin][idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
        }
    }
    std::vector<Eigen::VectorXcd> dnu_phij(3);
    for (int idim=0; idim<3; idim++)
    {
        dnu_phij[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> dnu_phiphi(3);
    for (int idim=0; idim<3; idim++)
    {
        dnu_phiphi[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXd> du(3);
    for (int idim=0; idim<3; idim++)
    {
        du[idim] = Eigen::VectorXd::Zero(plane_wave_basis.size_FFT_grid());
    }

    // divergence correction
    Eigen::Vector3d kVaux;

    for (int ispin2=0; ispin2<num_independent_spins; ispin2++) // spin index of x2
    {
        for (int ispin3=0; ispin3<2; ispin3++) // spin index of x3
        {
            orbital = bloch_states.density()[ispin3]; // density(R). Name "orbital" is meaningless here
            plane_wave_basis.FFT_forward(orbital, orbital); // -> density(G)

            // set \sum_q2 <*,*,q2| \nabla_2 u_23 |*,*,q2>
            for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
            {
                double uk = potentials.jastrow.uk(Gvect[ipw], ispin2, ispin3);
                for (int idim=0; idim<3; idim++)
                {
                    dnu[ispin2][idim](ipw) +=
                        uk * Gvect[ipw](idim) * orbital(ipw);
                }
            }
        }
        for (int idim=0; idim<3; idim++)
        {
            plane_wave_basis.FFT_backward(dnu[ispin2][idim], dnu[ispin2][idim]); // -> dnu(R)
        }
    } // ispin2

    // Since spin (x1) = spin (x2) for tc3b3, we denote it as "ispin" hereafter.
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            const Eigen::Vector3d kvector_ref = method.calc_mode()=="SCF" ? 
                kpoints.kvectors_scf()[ik][0] : kpoints.kvectors_band()[ik][0];
            const int num_G_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.num_G_at_k_scf()[ik] : plane_wave_basis.num_G_at_k_band()[ik];
            const Eigen::VectorXi Gindex_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.Gindex_at_k_scf()[ispin][ik][0] :  plane_wave_basis.Gindex_at_k_band()[ispin][ik][0];

            kvect = crystal_structure.reciprocal_vectors().transpose() * kvector_ref;

            for (int jband=0; jband<H3phi[ispin][ik].size(); jband++)
            {
                if (!parallelization.is_assigned_irreducible_kpoints_all_bands()[ispin][ik][jband]) { continue; }
                for (int jspinor=0; jspinor<num_spinor; jspinor++)
                {
                    // phi -> phij on the FFT-grid
                    plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, // isym = 0
                                                         jband, false, // time-rersal not used for isym=0
                                                         phi[ispin][ik][jband][jspinor], phij,
                                                         method.calc_mode());
                    plane_wave_basis.FFT_backward(phij, phij); // -> phij(R)

                    // \sum_q2 <*,*,q2| \nabla_2 u_23 |*,j,q2>
                    for (int idim=0; idim<3; idim++)
                    {
                        dnu_phij[idim] = dnu[ispin][idim].array() * phij.array();
                    }

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
                                double uk = potentials.jastrow.uk(kqvect + Gvect[ipw], ispin, ispin);
                                for (int idim=0; idim<3; idim++)
                                {
                                    du[idim](ipw) = (kqvect(idim) + Gvect[ipw](idim)) * uk;
                                }
                            }

                            for (int ibandq=0; ibandq<bloch_states.num_occupied_bands()[ispin][iq]; ibandq++)
                            {
                                // calculate phiq
                                Eigen::VectorXcd &chiq_ref = is_bitc ? chiq : phiq; // bra orbital

                                plane_wave_basis.get_orbital_FFTgrid(ispin, iq, isymq, ibandq,
                                                                     kpoints.is_time_reversal_used_at_k()[iq][isymq],
                                                                     bloch_states.phik_scf()[ispin][iq][ibandq][0], phiq,
                                                                     "SCF");
                                plane_wave_basis.FFT_backward(phiq, phiq);

                                if (is_bitc)
                                {
                                    plane_wave_basis.get_orbital_FFTgrid(ispin, iq, isymq, ibandq,
                                                                         kpoints.is_time_reversal_used_at_k()[iq][isymq],
                                                                         bloch_states.phik_left_scf()[ispin][iq][ibandq][0], chiq,
                                                                         "SCF");
                                    plane_wave_basis.FFT_backward(chiq, chiq);
                                }

                                // \sum_q2 <*,q1,q2| \nabla_2 u_23 |*,j,q2>
                                for (int idim=0; idim<3; idim++)
                                {
                                    dnu_phiphi[idim] = dnu_phij[idim].array() * chiq_ref.conjugate().array();
                                    plane_wave_basis.FFT_forward(dnu_phiphi[idim], dnu_phiphi[idim]);
                                }

                                // \sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |*,j,q2>
                                orbital = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid()); // "orbital" = meaningless name
                                for (int idim=0; idim<3; idim++)
                                {
                                    orbital = orbital.array() + 
                                        du[idim].array() * dnu_phiphi[idim].array();
                                }
                                plane_wave_basis.FFT_backward(orbital, orbital);

                                // \sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |q1,j,q2>
                                Hphij_sub = Hphij_sub.array() +
                                    bloch_states.filling()[ispin][iq][ibandq] * orbital.array() * phiq.array();
                            } // ibandq
                        } // isymq
                    } // iq
                    plane_wave_basis.FFT_forward(Hphij_sub, Hphij_sub);
                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                    {
                        H3phi[ispin][ik][jband][0](ipw_at_k) += 
                            three_body_factor * Hphij_sub(Gindex_at_k(ipw_at_k));
                    }

                    // divergence correction for \nabla u
                    if (potentials.includes_div_correction()) 
                    {
                        kVaux = method.calc_mode()=="BAND" ?
                            -potentials.sum_of_kVaux_band()[ik] :
                            -potentials.sum_of_kVaux_scf()[ik];
                        //  Here, "-1" comes from: \int \nabla Vaux - \sum \nabla Vaux = 0 - sum_of_kVaux = "-1" * sum_of_kVaux.
                        kVaux *= potentials.jastrow.A_long()[ispin][ispin]; // since we consider not 1/r but A/r-like divergence

                        if (kVaux.squaredNorm() > 1e-8) // if = 0 then no correction is needed
                        {
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

                                const Eigen::VectorXcd &phik_ref = method.calc_mode()=="BAND" ?
                                    bloch_states.phik_band()[ispin][ik][ibandk][0] :
                                    bloch_states.phik_scf()[ispin][ik][ibandk][0];
                                if (is_bitc)
                                {
                                    const Eigen::VectorXcd &chik_ref = method.calc_mode()=="BAND" ?
                                        bloch_states.phik_left_band()[ispin][ik][ibandk][0] :
                                        bloch_states.phik_left_scf()[ispin][ik][ibandk][0];

                                    plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, ibandk, false,
                                                                         chik_ref, phiq, method.calc_mode());
                                }
                                else
                                {
                                    plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, // isym = 0
                                                                         ibandk, false, // time-rersal not used for isym=0
                                                                         phik_ref, phiq, method.calc_mode());
                                }
                                plane_wave_basis.FFT_backward(phiq, phiq);

                                Eigen::Vector3cd coeff_vec_3b3 = {0.0, 0.0, 0.0};
                                for (int idim=0; idim<3; idim++)
                                {
                                    coeff_vec_3b3[idim] = phiq.dot(dnu_phij[idim]); // conj(phiq)*dnu_phij
                                }
                                coeff_vec_3b3 *= filling_ibandk * three_body_factor /
                                    plane_wave_basis.size_FFT_grid(); // integral in r-space: normalized const.

                                Complex temp_div_3b3 = 0.0;
                                for (int idim=0; idim<3; idim++)
                                {
                                    temp_div_3b3 += coeff_vec_3b3[idim] * kVaux[idim];
                                }

                                H3phi[ispin][ik][jband][0] += temp_div_3b3 * phik_ref;
                            } // ibandk
                        } // if (kVaux != 0)
                    } // if (includes_div_correction)
                } // jspinor
            } // jband
        } // ik
    } // ispin
}
