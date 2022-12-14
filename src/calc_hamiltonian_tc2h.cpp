// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

// two-body, Hartree term (TC or BiTC)
// calculate sum_q <*,q| 1/r + \nabla^2 u - (\nabla u)^2 |j,q> (2a_h)
// calculate sum_q <*,q| \nabla_1 u_12 \nabla_1 |j,q> (2b_h1)
// calculate sum_q <*,q| \nabla_2 u_21 \nabla_2 |j,q> (2b_h2)
void calc_hamiltonian::tc2h(const Parallelization &parallelization, 
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
    const bool is_bitc = method.calc_method()=="BITC" ? true : false;
    const int num_independent_spins = spin.num_independent_spins();
    const int num_irreducible_kpoints_scf = kpoints.num_irreducible_kpoints_scf();
    const int num_irreducible_kpoints = method.calc_mode()=="SCF" ? 
        kpoints.num_irreducible_kpoints_scf() : kpoints.num_irreducible_kpoints_band();
    const int num_spinor = (spin.is_spinor()==false) ? 1 : 2;

    const double two_body_factor = 1.0/(kpoints.num_kpoints()*crystal_structure.unit_cell_volume());
    // coefficient of [2b_h1]: FFT I^2 = -1 will be is additionally multiplied 
    // coefficient of [2b_h2]: FFT I^2 * convolution (\nabla_2 u_21 = "-1" * \nabla_1 u_12) = +1.
    //                         But since we consider both u_para and u_anti, will be divided with 2 when num_independent_spins==1.
    //                         *Note* density in 2a_h and 2b_h1 are already divided with 2 (to say, density(up) or density(down)).

    std::vector<Eigen::Vector3d> Gvect(plane_wave_basis.size_FFT_grid());
    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
    {
        Gvect[ipw] = crystal_structure.reciprocal_vectors().transpose()
            * plane_wave_basis.get_Gvector(ipw).cast<double>();
    }

    Eigen::Vector3d kvect, qvect;
    Eigen::VectorXcd orbital(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd V_orbital(plane_wave_basis.size_FFT_grid());
    std::vector<Eigen::VectorXcd> grad_orbital(3);
    for (int idim=0; idim<3; idim++) 
    {
        grad_orbital[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }

    // [2b_h2] phiq_gphiq = sum_q <q|\nabla|q>
    std::vector<std::vector<Eigen::VectorXcd> > phiq_gphiq(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++) 
    {
        phiq_gphiq[ispin].resize(3);
        for (int idim=0; idim<3; idim++) 
        {
            phiq_gphiq[ispin][idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid()); 
        }
    }
    for (int ispin=0; ispin<num_independent_spins; ispin++) 
    {
        for (int iq=0; iq<num_irreducible_kpoints_scf; iq++)
        {
            for (int isymq=0; isymq<kpoints.kvectors_scf()[iq].size(); isymq++)
            {
                if (bloch_states.num_occupied_bands()[ispin][iq]==0) { continue; }
                qvect = crystal_structure.reciprocal_vectors().transpose() *
                    kpoints.kvectors_scf()[iq][isymq];

                const int nbands_old = bloch_states.filling_old()[ispin][iq].size();
                for (int ibandq=-nbands_old; ibandq<bloch_states.num_occupied_bands()[ispin][iq]; ibandq++)
                {
                    if (ibandq>=0)
                    {
                        plane_wave_basis.get_orbital_FFTgrid(ispin, iq, isymq,
                                                             kpoints.is_time_reversal_used_at_k()[iq][isymq],
                                                             bloch_states.phik_scf()[ispin][iq][ibandq][0],
                                                             orbital, "SCF");
                    }
                    else
                    {
                        plane_wave_basis.get_orbital_FFTgrid(ispin, iq, isymq,
                                                             kpoints.is_time_reversal_used_at_k()[iq][isymq],
                                                             bloch_states.phik_scf_old()[ispin][iq][-1-ibandq][0],
                                                             orbital, "SCF");
                    }
                    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                    {
                        for (int idim=0; idim<3; idim++)
                        {
                            grad_orbital[idim](ipw) = (qvect(idim) + Gvect[ipw](idim)) * orbital(ipw);
                        }
                    }
                    if (is_bitc)
                    {
                        if (ibandq>=0)
                        {
                            plane_wave_basis.get_orbital_FFTgrid(ispin, iq, isymq,
                                                                 kpoints.is_time_reversal_used_at_k()[iq][isymq],
                                                                 bloch_states.phik_left_scf()[ispin][iq][ibandq][0],
                                                                 orbital, "SCF");
                        }
                        else
                        {
                            plane_wave_basis.get_orbital_FFTgrid(ispin, iq, isymq,
                                                                 kpoints.is_time_reversal_used_at_k()[iq][isymq],
                                                                 bloch_states.phik_left_scf_old()[ispin][iq][-1-ibandq][0],
                                                                 orbital, "SCF");
                        }
                    }
                    plane_wave_basis.FFT_backward(orbital, orbital); // -> phi_q(R)
                    double filq = ibandq>=0 ? bloch_states.filling()[ispin][iq][ibandq] : bloch_states.filling_old()[ispin][iq][-1-ibandq];
                    orbital = filq * orbital.conjugate();
                    for (int idim=0; idim<3; idim++)
                    {
                        plane_wave_basis.FFT_backward(grad_orbital[idim], grad_orbital[idim]);
                        phiq_gphiq[ispin][idim] = phiq_gphiq[ispin][idim].array() +
                            orbital.array() * grad_orbital[idim].array(); // <q|\nabla|q>
                    }
                } // ibandq
            } // isymq
        } // iq
        for (int idim=0; idim<3; idim++)
        {
            plane_wave_basis.FFT_forward(phiq_gphiq[ispin][idim], phiq_gphiq[ispin][idim]); // -> (G)
        }
    } // ispin

    // [2b_h1] dnu = sum_q <*,q| \nabla_1 u_12 | *,q>
    std::vector<std::vector<Eigen::VectorXcd> > dnu(num_independent_spins); // spin index of x1
    for (int ispin=0; ispin<num_independent_spins; ispin++) 
    {
        dnu[ispin].resize(3);
        for (int idim=0; idim<3; idim++)
        {
            dnu[ispin][idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
        }
    }
    // [2a_h + 2b_h2] V_hartree = sum_q <*,q| tc_2body(1,2) + \nabla_2 u_21 \nabla_2 |*,q>
    std::vector<Eigen::VectorXcd> V_2body(num_independent_spins); 
    for (int ispin=0; ispin<num_independent_spins; ispin++) 
    {
        V_2body[ispin] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    double spin_factor = num_independent_spins==1 ? 2.0 : 1.0; // canceled with filling=2.0 for num_independent_spins==1
    for (int ispin=0; ispin<num_independent_spins; ispin++) 
    {
        for (int jspin=0; jspin<2; jspin++)
        {
            int jspin_ref = num_independent_spins==1 ? 0 : jspin;

            orbital = bloch_states.density()[jspin]; // density(R). Name "orbital" is meaningless here
            plane_wave_basis.FFT_forward(orbital, orbital); // -> density(G)

            for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
            {
                // [2a_h] sum_q <*,q| tc_2body(1,2) |*,q>
                V_2body[ispin](ipw) += 
                    potentials.jastrow.tc_2body(Gvect[ipw], ispin, jspin) * orbital(ipw);

                double uk = potentials.jastrow.uk(Gvect[ipw], ispin, jspin);
                for (int idim=0; idim<3; idim++)
                {

                    // [2b_h2] sum_q <*,q| \nabla_2 u_21 \nabla_2 |*,q>
                    V_2body[ispin](ipw) += 
                        uk * Gvect[ipw](idim) * phiq_gphiq[jspin_ref][idim](ipw) / spin_factor;

                    // [2b_h1] sum_q <*,q| \nabla_1 u_12 |*,q>
                    dnu[ispin][idim](ipw) +=
                        uk * Gvect[ipw](idim) * orbital(ipw);
                }
            }
        }
        plane_wave_basis.FFT_backward(V_2body[ispin], V_2body[ispin]); // -> V_2body(R)
        for (int idim=0; idim<3; idim++)
        {
            plane_wave_basis.FFT_backward(dnu[ispin][idim], dnu[ispin][idim]); // -> dnu(R)
        }
    }

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

            for (int jband=0; jband<H2phi[ispin][ik].size(); jband++)
            {
                if (!parallelization.is_assigned_irreducible_kpoints_all_bands()[ispin][ik][jband]) { continue; }
                for (int jspinor=0; jspinor<num_spinor; jspinor++)
                {
                    // phi -> orbital on the FFT-grid
                    plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, // isym = 0
                                                         false, // time-rersal not used for isym=0
                                                         phi[ispin][ik][jband][jspinor], orbital,
                                                         method.calc_mode());
                    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                    {
                        for (int idim=0; idim<3; idim++)
                        {
                            grad_orbital[idim](ipw) = (kvect(idim) + Gvect[ipw](idim)) * orbital(ipw);
                        }
                    }
                    for (int idim=0; idim<3; idim++)
                    {
                        plane_wave_basis.FFT_backward(grad_orbital[idim], grad_orbital[idim]);
                    }
                    plane_wave_basis.FFT_backward(orbital, orbital); // -> phi_j(R)

                    // [2a_h + 2b_h2]
                    V_orbital = V_2body[ispin].array() * orbital.array();
                    // [2b_h1]
                    for (int idim=0; idim<3; idim++)
                    {
                        V_orbital = V_orbital.array() -
                            grad_orbital[idim].array() * dnu[ispin][idim].array(); // I^2 = -1
                    }

                    plane_wave_basis.FFT_forward(V_orbital, V_orbital);
                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                    {
                        H2phi[ispin][ik][jband][jspinor](ipw_at_k)
                            += two_body_factor * V_orbital(Gindex_at_k(ipw_at_k));
                    }
                } // jspinor
            } // jband
        } // ik
    } // ispin
}
