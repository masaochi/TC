// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

// *** tc3a2 ***
// calculate 0.5*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |j,q2,q1>
void calc_hamiltonian::tc3a2(const Parallelization &parallelization, 
                             const Method &method,
                             const CrystalStructure &crystal_structure,
                             const Symmetry &symmetry, const Potentials &potentials,
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
    const std::vector<int> size_FFT_grid_vec = plane_wave_basis.size_FFT_grid_vec();

    // (0.5) * I^2 = -0.5
    double three_body_factor = -0.5/(kpoints.num_kpoints()*kpoints.num_kpoints()*
                                     crystal_structure.unit_cell_volume()*
                                     crystal_structure.unit_cell_volume());
    // spin indices of x2,x3 should be the same: filling should be divided with 2
    // Also, "u_para" and "u_anti" are summed over. : divided with 2
    if (num_independent_spins==1) { three_body_factor /= 4.0; } 

    std::vector<Eigen::Vector3d> Gvect(plane_wave_basis.size_FFT_grid());
    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
    {
        Gvect[ipw] = crystal_structure.reciprocal_vectors().transpose()
            * plane_wave_basis.get_Gvector(ipw).cast<double>();
    }

    Eigen::Vector3d q1vect, q2vect, q2q1vect;
    Eigen::VectorXcd orbital(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd orbital2(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd phij(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd phiq1(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd phiq2(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd chiq1;
    if (is_bitc) { chiq1.resize(plane_wave_basis.size_FFT_grid()); }
    Eigen::VectorXcd chiq2;
    if (is_bitc) { chiq2.resize(plane_wave_basis.size_FFT_grid()); }
    Eigen::VectorXcd phiphi(plane_wave_basis.size_FFT_grid());

    std::vector<std::vector<Eigen::VectorXd> > du(2);
    for (int ispin=0; ispin<2; ispin++)
    {
        du[ispin].resize(3);
        for (int idim=0; idim<3; idim++)
        {
            du[ispin][idim] = Eigen::VectorXd::Zero(plane_wave_basis.size_FFT_grid());
        }
    }
    std::vector<Eigen::VectorXcd> qsum_temp_3a2(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        qsum_temp_3a2[ispin] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> qsum_temp_3a2_temp(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        qsum_temp_3a2_temp[ispin] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> V_3a2(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        V_3a2[ispin] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }

    for (int ispinq=0; ispinq<num_independent_spins; ispinq++) // spin of |q1> and |q2>
    {
        for (int iq2=0; iq2<num_irreducible_kpoints_scf; iq2++)
        {
            q2vect = crystal_structure.reciprocal_vectors().transpose() *
                kpoints.kvectors_scf()[iq2][0];
            for (int ibandq2=0; ibandq2<bloch_states.num_occupied_bands()[ispinq][iq2]; ibandq2++)
            {
                Eigen::VectorXcd &chiq2_ref = is_bitc ? chiq2 : phiq2; // bra orbital
                
                plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, 0, ibandq2,
                                                     kpoints.is_time_reversal_used_at_k()[iq2][0],
                                                     bloch_states.phik_scf()[ispinq][iq2][ibandq2][0], phiq2,
                                                     "SCF");
                plane_wave_basis.FFT_backward(phiq2, phiq2);
                
                if (is_bitc)
                {
                    plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, 0, ibandq2,
                                                         kpoints.is_time_reversal_used_at_k()[iq2][0],
                                                         bloch_states.phik_left_scf()[ispinq][iq2][ibandq2][0], chiq2,
                                                         "SCF");
                    plane_wave_basis.FFT_backward(chiq2, chiq2);
                }

                for (int ispin=0; ispin<num_independent_spins; ispin++)
                {                
                    qsum_temp_3a2[ispin] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                }
                for (int iq1=0; iq1<num_irreducible_kpoints_scf; iq1++)
                {
                    for (int isymq1=0; isymq1<kpoints.kvectors_scf()[iq1].size(); isymq1++)
                    {
                        bool is_q1_assigned = false;
                        for (int ibandq1=0; ibandq1<bloch_states.num_occupied_bands()[ispinq][iq1]; ibandq1++)
                        {
                            if (parallelization.is_assigned_all_kpoints_occupied_bands()[ispinq][iq1][isymq1][ibandq1]) 
                            {
                                is_q1_assigned = true;
                                break;
                            }
                        }
                        if (!is_q1_assigned) { continue; } // if no ibandq1 is assinged to this node
                        
                        q1vect = crystal_structure.reciprocal_vectors().transpose() *
                            kpoints.kvectors_scf()[iq1][isymq1];
                        q2q1vect = q2vect - q1vect;
                        
                        for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                        {
                            double uk0 = potentials.jastrow.uk(q2q1vect + Gvect[ipw], ispinq, 0);
                            double uk1 = potentials.jastrow.uk(q2q1vect + Gvect[ipw], ispinq, 1);
                            for (int idim=0; idim<3; idim++)
                            {
                                du[0][idim](ipw) = (q2q1vect(idim) + Gvect[ipw](idim)) * uk0;
                                du[1][idim](ipw) = (q2q1vect(idim) + Gvect[ipw](idim)) * uk1;
                            }
                        }
                        for (int ibandq1=0; ibandq1<bloch_states.num_occupied_bands()[ispinq][iq1]; ibandq1++)
                        {
                            if (!parallelization.is_assigned_all_kpoints_occupied_bands()[ispinq][iq1][isymq1][ibandq1]) { continue; }
                            
                            Eigen::VectorXcd &chiq1_ref = is_bitc ? chiq1 : phiq1; // bra orbital
                            
                            plane_wave_basis.get_orbital_FFTgrid(ispinq, iq1, isymq1, ibandq1,
                                                                 kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                 bloch_states.phik_scf()[ispinq][iq1][ibandq1][0], phiq1,
                                                                 "SCF");
                            plane_wave_basis.FFT_backward(phiq1, phiq1);
                            
                            if (is_bitc)
                            {
                                plane_wave_basis.get_orbital_FFTgrid(ispinq, iq1, isymq1, ibandq1,
                                                                     kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                     bloch_states.phik_left_scf()[ispinq][iq1][ibandq1][0], chiq1,
                                                                     "SCF");
                                plane_wave_basis.FFT_backward(chiq1, chiq1);
                            }
                            
                            // <*,q1,*| |*,q2,*>
                            phiphi = chiq1_ref.conjugate().array() * phiq2.array();
                            plane_wave_basis.FFT_forward(phiphi, phiphi);

                            // (<*,*,q2| |*,*,q1>)^* = <*,*,phi_q1| |*,*,chi_q2>
                            if (is_bitc)
                            {
                                phij = phiq1.conjugate().array() * chiq2_ref.array();
                                plane_wave_basis.FFT_forward(phij, phij);
                            }

                            for (int ispin=0; ispin<2; ispin++)
                            {
                                int ispin_ref = num_independent_spins==1 ? 0 : ispin;
                                for (int idim=0; idim<3; idim++)
                                {
                                    // <*,q1,*| \nabla_1 u_12 |*,q2,*>
                                    orbital = du[ispin][idim].array() * phiphi.array();
                                    plane_wave_basis.FFT_backward(orbital, orbital); // idim-component

                                    if (is_bitc)
                                    {
                                        // <*,*,phi_q1| \nabla_1 u_13 |*,*,chi_q2>
                                        orbital2 = du[ispin][idim].array() * phij.array();
                                        plane_wave_basis.FFT_backward(orbital2, orbital2); // idim-component

                                        // <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |*,q2,q1>
                                        qsum_temp_3a2[ispin_ref] = qsum_temp_3a2[ispin_ref].array() - // -: conj(I) for orbital2
                                            bloch_states.filling()[ispinq][iq1][ibandq1] *
                                            orbital.array() * orbital2.conjugate().array();
                                    }
                                    else // orbital2 = orbital
                                    {
                                        // <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |*,q2,q1>
                                        qsum_temp_3a2[ispin_ref] = qsum_temp_3a2[ispin_ref].array() - // -: conj(I) for orbital2
                                            bloch_states.filling()[ispinq][iq1][ibandq1] *
                                            orbital.array().abs2();
                                    } // if (is_bitc)
                                } // idim
                            } // ispin
                        } // ibandq1
                    } // isymq1
                } // iq1
                for (int ispin=0; ispin<num_independent_spins; ispin++)
                {
                    MPI_Allreduce(qsum_temp_3a2[ispin].data(),
                                  qsum_temp_3a2_temp[ispin].data(),
                                  qsum_temp_3a2_temp[ispin].size(),
                                  MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
                    plane_wave_basis.FFT_forward(qsum_temp_3a2_temp[ispin], qsum_temp_3a2_temp[ispin]);
                }
            
                // divergence correction. Note: filling(iq2) will be multiplied when calculating V_3a2,
                // but the filling appearing below is filling(iq1) * delta_{q1,q2}.
                if (potentials.includes_div_correction())
                {
                    // \int d^3q/(2pi)^3 4pi*exp(-alpha q^2)/q^2 = 2/pi \int dq_0^infty exp(-alpha q^2)/q^2 = 1/sqrt(pi*alpha)
                    // Later divided by "num_kpoints * unit_cell_volume" (three_body_factor)
                    double div_corr =
                        1.0/std::sqrt(PI*potentials.alpha_Vaux()) * kpoints.num_kpoints() * crystal_structure.unit_cell_volume() 
                        - potentials.sum_of_Vaux_scf()[iq2];

                    // <chi|phi>
                    phiphi = chiq2_ref.conjugate().array() * phiq2.array();
                    plane_wave_basis.FFT_forward(phiphi, phiphi);

                    phiphi *= (bloch_states.filling()[ispinq][iq2][ibandq2] *
                               FourPI * 2); // *2: each Jastrow shows divergence in <...| \nabla u \nabla u |...>

                    for (int ispin=0; ispin<2; ispin++)
                    {
                        int ispin_ref = num_independent_spins==1 ? 0 : ispin;

                        qsum_temp_3a2_temp[ispin_ref](0) -= // const component (G=0)
                            FourPI * potentials.jastrow.A_long()[ispinq][ispin] *
                            (potentials.jastrow.A_long()[ispinq][ispin] * div_corr +
                             potentials.jastrow.uk({0.0, 0.0, 0.0}, ispinq, ispin)) *
                            bloch_states.filling()[ispinq][iq2][ibandq2];

                        for (int ipw=1; ipw<plane_wave_basis.size_FFT_grid(); ipw++) // no G=0 term! (ipw!=0)
                        {
                            double uk = potentials.jastrow.uk(Gvect[ipw], ispinq, ispin);

                            qsum_temp_3a2_temp[ispin_ref](ipw) -= phiphi(ipw) *
                                potentials.jastrow.A_long()[ispinq][ispin] * uk;
                        }
                    }
                } // if (includes_div_correction)

                // symmetry transformation of qsum_temp_3a4: (iq2, 0) -> (iq2, isymq2)
                // (k'(iq2,isymq2) = U^+ k(iq2,0))
                for (int isymq2=0; isymq2<kpoints.kvectors_scf()[iq2].size(); isymq2++)
                {
                    int ind_sym = kpoints.index_of_rotation_at_k()[iq2][isymq2];
                    Eigen::Matrix3i Udag = symmetry.rotation()[ind_sym].transpose(); // U^+
                    Eigen::Vector3d translation = symmetry.translation()[ind_sym];                    

                    Eigen::Vector3i Gvector, Gvector_trans;
                    for (int ispin=0; ispin<num_independent_spins; ispin++)
                    {
                        qsum_temp_3a2[ispin] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                    }
                    if (kpoints.is_time_reversal_used_at_k()[iq2][isymq2])
                    {
                        for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                        {
                            Gvector = plane_wave_basis.get_Gvector(ipw);
                            Gvector_trans = -Udag * Gvector; // U^+ G                            
                            for (int idim=0; idim<3; idim++)
                            {
                                while (Gvector_trans(idim)<0) { Gvector_trans(idim) += size_FFT_grid_vec[idim]; }
                                while (Gvector_trans(idim)>=size_FFT_grid_vec[idim]) { Gvector_trans(idim) -= size_FFT_grid_vec[idim]; }
                            }
                            Complex phase = std::cos(2*PI*translation.dot(Gvector_trans.cast<double>())) +
                                I*std::sin(2*PI*translation.dot(Gvector_trans.cast<double>())); // (U^+ G) * r0

                            int ipw_trans = Gvector_trans(0)
                                + Gvector_trans(1) * size_FFT_grid_vec[0]
                                + Gvector_trans(2) * size_FFT_grid_vec[0] * size_FFT_grid_vec[1];
                            for (int ispin=0; ispin<num_independent_spins; ispin++)
                            {
                                qsum_temp_3a2[ispin](ipw_trans) += phase *
                                    std::conj(qsum_temp_3a2_temp[ispin](ipw));
                            }
                        } // ipw
                    }
                    else // time-reversal symmetry is not used for this k-point
                    {
                        for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                        {
                            Gvector = plane_wave_basis.get_Gvector(ipw);
                            Gvector_trans = Udag * Gvector; // U^+ G                           
                            for (int idim=0; idim<3; idim++)
                            {
                                while (Gvector_trans(idim)<0) { Gvector_trans(idim) += size_FFT_grid_vec[idim]; }
                                while (Gvector_trans(idim)>=size_FFT_grid_vec[idim]) { Gvector_trans(idim) -= size_FFT_grid_vec[idim]; }
                            }
                            Complex phase = std::cos(2*PI*translation.dot(Gvector_trans.cast<double>())) +
                                I*std::sin(2*PI*translation.dot(Gvector_trans.cast<double>())); // (U^+ G) * r0

                            int ipw_trans = Gvector_trans(0)
                                + Gvector_trans(1) * size_FFT_grid_vec[0]
                                + Gvector_trans(2) * size_FFT_grid_vec[0] * size_FFT_grid_vec[1];
                            for (int ispin=0; ispin<num_independent_spins; ispin++)
                            {
                                qsum_temp_3a2[ispin](ipw_trans) += phase *
                                    qsum_temp_3a2_temp[ispin](ipw);
                            }
                        } // ipw
                    } // if (time_reversal)

                    for (int ispin=0; ispin<num_independent_spins; ispin++)
                    {
                        V_3a2[ispin] = V_3a2[ispin].array() +
                            bloch_states.filling()[ispinq][iq2][ibandq2] * qsum_temp_3a2[ispin].array();
                    }
                } // isymq2
            } // ibandq2
        } // iq2
    } // ispinq
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        plane_wave_basis.FFT_backward(V_3a2[ispin], V_3a2[ispin]);
    }
    
    for (int ispin=0; ispin<num_independent_spins; ispin++) // spin of |j>
    {
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            const int num_G_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.num_G_at_k_scf()[ik] : plane_wave_basis.num_G_at_k_band()[ik];
            const Eigen::VectorXi Gindex_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.Gindex_at_k_scf()[ispin][ik][0] :  plane_wave_basis.Gindex_at_k_band()[ispin][ik][0];
                        
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

                    orbital = V_3a2[ispin].array() * phij.array();
                    plane_wave_basis.FFT_forward(orbital, orbital);
                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                    {
                        H3phi[ispin][ik][jband][0](ipw_at_k) += 
                            three_body_factor *
                            orbital(Gindex_at_k(ipw_at_k));
                    }
                } // jspinor
            } // jband
        } // ik
    } // ispin
}
