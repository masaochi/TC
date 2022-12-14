// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

// [3a2] 0.5*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |j,q2,q1>
// [3a4] -0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q1,q2,j>
//   *equivalent then including: [3a5] -0.5*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q2,j,q1>
// [3b2] 0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |j,q2,q1>
// [3b5] -0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |q2,j,q1>
void calc_hamiltonian::tc3a2a4b2b5(const Parallelization &parallelization, 
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

    double three_body_factor_base = 1.0/(kpoints.num_kpoints()*kpoints.num_kpoints()*
                                         crystal_structure.unit_cell_volume()*
                                         crystal_structure.unit_cell_volume());

    double three_body_factor_a2 = -0.5*three_body_factor_base; // (0.5) * I^2 * -0.5
    if (num_independent_spins==1) { three_body_factor_a2 /= 4.0; } // spin indices of x2,x3 are the same (/2.0). u_para/anti are summed over (/2.0).

    double three_body_factor_a4 = three_body_factor_base; // (-0.5*2) * I^2 = +1.0
    if (num_independent_spins==1) { three_body_factor_a4 /= 4.0; } // spin indices of x1,x2,x3 are the same

    double three_body_factor_b2 = three_body_factor_base; // (0.5*2) * I^2 * (-1) (convolution \nabla_2 u_21 = -\nabla_1 u_12) = +1.0
    if (num_independent_spins==1) { three_body_factor_b2 /= 4.0; } // spin indices of x2,x3 are the same (/2.0). u_para/anti are summed over (/2.0).

    double three_body_factor_b5 = -three_body_factor_base; // (-0.5*2) * I^2 * (-1) (convolution \nabla_2 u_21 = -\nabla_1 u_12) = -1.0
    if (num_independent_spins==1) { three_body_factor_b5 /= 4.0; } // spin indices of x1,x2,x3 are the same

    // When you would like to debug each term:
//    three_body_factor_a2 = 0.0;
//    three_body_factor_a4 = 0.0;
//    three_body_factor_b2 = 0.0;
//    three_body_factor_b5 = 0.0;

    std::vector<Eigen::Vector3d> Gvect(plane_wave_basis.size_FFT_grid());
    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
    {
        Gvect[ipw] = crystal_structure.reciprocal_vectors().transpose()
            * plane_wave_basis.get_Gvector(ipw).cast<double>();
    }

    Eigen::Vector3d kvect, q1vect, q2vect, q2symvect;
    Eigen::Vector3d kqvect, qkvect, q2q1vect;
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
    std::vector<Eigen::VectorXcd> qsum_temp_3a4(3);
    for (int idim=0; idim<3; idim++)
    {
        qsum_temp_3a4[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> qsum_temp_temp(3); // also used as qsum_temp_temp(num_independent_spins)
    for (int idim=0; idim<3; idim++)
    {
        qsum_temp_temp[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> V_3a2(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        V_3a2[ispin] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> qsum_temp_3b5(3);
    for (int idim=0; idim<3; idim++)
    {
        qsum_temp_3b5[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> qsum_temp_3b5_temp(3);
    for (int idim=0; idim<3; idim++)
    {
        qsum_temp_3b5_temp[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> qsum_temp_3b2(3);
    for (int idim=0; idim<3; idim++)
    {
        qsum_temp_3b2[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }

    // divergence correction
    Eigen::Vector3d kVaux;

    // Start!
    for (int ispinq=0; ispinq<num_independent_spins; ispinq++) // spin of |q1> and |q2>
    {
        for (int idim=0; idim<3; idim++)
        {
            qsum_temp_3b2[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
        }
        for (int iq2=0; iq2<num_irreducible_kpoints_scf; iq2++)
        {
            q2vect = crystal_structure.reciprocal_vectors().transpose() *
                kpoints.kvectors_scf()[iq2][0];
            const int nbands2_old = bloch_states.filling_old()[ispinq][iq2].size();
            for (int ibandq2=-nbands2_old; ibandq2<bloch_states.num_occupied_bands()[ispinq][iq2]; ibandq2++)
            {
                double filq2 = ibandq2>=0 ? 
                    bloch_states.filling()[ispinq][iq2][ibandq2] : bloch_states.filling_old()[ispinq][iq2][-1-ibandq2];

                Eigen::VectorXcd &chiq2_ref = is_bitc ? chiq2 : phiq2; // bra orbital
                
                if (ibandq2>=0)
                {
                    plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, 0,
                                                         kpoints.is_time_reversal_used_at_k()[iq2][0],
                                                         bloch_states.phik_scf()[ispinq][iq2][ibandq2][0], phiq2,
                                                         "SCF");
                }
                else
                {
                    plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, 0,
                                                         kpoints.is_time_reversal_used_at_k()[iq2][0],
                                                         bloch_states.phik_scf_old()[ispinq][iq2][-1-ibandq2][0], phiq2,
                                                         "SCF");
                }
                plane_wave_basis.FFT_backward(phiq2, phiq2);
                
                if (is_bitc)
                {
                    if (ibandq2>=0)
                    {
                        plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, 0,
                                                             kpoints.is_time_reversal_used_at_k()[iq2][0],
                                                             bloch_states.phik_left_scf()[ispinq][iq2][ibandq2][0], chiq2,
                                                             "SCF");
                    }
                    else
                    {
                        plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, 0,
                                                             kpoints.is_time_reversal_used_at_k()[iq2][0],
                                                             bloch_states.phik_left_scf_old()[ispinq][iq2][-1-ibandq2][0], chiq2,
                                                             "SCF");
                    }
                    plane_wave_basis.FFT_backward(chiq2, chiq2);
                }

                for (int ispin=0; ispin<num_independent_spins; ispin++)
                {                
                    qsum_temp_3a2[ispin] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                }
                for (int idim=0; idim<3; idim++)
                {
                    qsum_temp_3a4[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                    qsum_temp_3b5[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                }
                for (int iq1=0; iq1<num_irreducible_kpoints_scf; iq1++)
                {
                    for (int isymq1=0; isymq1<kpoints.kvectors_scf()[iq1].size(); isymq1++)
                    {
                        bool is_q1_assigned = false;
                        const int nbands1_old = bloch_states.filling_old()[ispinq][iq1].size();
                        for (int ibandq1=-nbands1_old; ibandq1<bloch_states.num_occupied_bands()[ispinq][iq1]; ibandq1++)
                        {
                            if (ibandq1>=0)
                            {
                                if (parallelization.is_assigned_all_kpoints_occupied_bands()[ispinq][iq1][isymq1][ibandq1]) 
                                {
                                    is_q1_assigned = true;
                                    break;
                                }
                            }
                            else
                            {
                                if (parallelization.is_assigned_all_kpoints_occupied_bands_old()[ispinq][iq1][isymq1][-1-ibandq1]) 
                                {
                                    is_q1_assigned = true;
                                    break;
                                }
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
                        for (int ibandq1=-nbands1_old; ibandq1<bloch_states.num_occupied_bands()[ispinq][iq1]; ibandq1++)
                        {
                            if (ibandq1>=0)
                            {
                                if (!parallelization.is_assigned_all_kpoints_occupied_bands()[ispinq][iq1][isymq1][ibandq1]) { continue; }
                            }
                            else
                            {
                                if (!parallelization.is_assigned_all_kpoints_occupied_bands_old()[ispinq][iq1][isymq1][-1-ibandq1]) { continue; }
                            }

                            Eigen::VectorXcd &chiq1_ref = is_bitc ? chiq1 : phiq1; // bra orbital
                            
                            if (ibandq1>=0)
                            {
                                plane_wave_basis.get_orbital_FFTgrid(ispinq, iq1, isymq1,
                                                                     kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                     bloch_states.phik_scf()[ispinq][iq1][ibandq1][0], phiq1,
                                                                     "SCF");
                            }
                            else
                            {
                                plane_wave_basis.get_orbital_FFTgrid(ispinq, iq1, isymq1,
                                                                     kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                     bloch_states.phik_scf_old()[ispinq][iq1][-1-ibandq1][0], phiq1,
                                                                     "SCF");
                            }
                            plane_wave_basis.FFT_backward(phiq1, phiq1);
                            
                            if (is_bitc)
                            {
                                if (ibandq1>=0)
                                {
                                    plane_wave_basis.get_orbital_FFTgrid(ispinq, iq1, isymq1,
                                                                         kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                         bloch_states.phik_left_scf()[ispinq][iq1][ibandq1][0], chiq1,
                                                                         "SCF");
                                }
                                else
                                {
                                    plane_wave_basis.get_orbital_FFTgrid(ispinq, iq1, isymq1,
                                                                         kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                         bloch_states.phik_left_scf_old()[ispinq][iq1][-1-ibandq1][0], chiq1,
                                                                         "SCF");
                                }
                                plane_wave_basis.FFT_backward(chiq1, chiq1);
                            }
                            
                            // <*,q1,*| |*,q2,*>
                            phiphi = chiq1_ref.conjugate().array() * phiq2.array();
                            plane_wave_basis.FFT_forward(phiphi, phiphi);

                            if (is_bitc)
                            {
                                // (<*,*,q2| |*,*,q1>)^* = <*,*,phi_q1| |*,*,chi_q2>
                                phij = phiq1.conjugate().array() * chiq2_ref.array();
                                plane_wave_basis.FFT_forward(phij, phij); // = phiphi (for !is_bitc)
                            }

                            double filq1 = ibandq1>=0 ? 
                                bloch_states.filling()[ispinq][iq1][ibandq1] : bloch_states.filling_old()[ispinq][iq1][-1-ibandq1];

                            for (int ispin=0; ispin<2; ispin++)
                            {
                                int ispin_ref = num_independent_spins==1 ? 0 : ispin;
                                for (int idim=0; idim<3; idim++)
                                {
                                    // <*,q1,*| \nabla_1 u_12 |*,q2,*>
                                    orbital = du[ispin][idim].array() * phiphi.array();
                                    plane_wave_basis.FFT_backward(orbital, orbital); // idim-component

                                    /*** [3a2] ***/
                                    if (is_bitc)
                                    {
                                        // <*,*,phi_q1| \nabla_1 u_13 |*,*,chi_q2>
                                        orbital2 = du[ispin][idim].array() * phij.array();
                                        plane_wave_basis.FFT_backward(orbital2, orbital2); // idim-component

                                        // <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |*,q2,q1>
                                        qsum_temp_3a2[ispin_ref] = qsum_temp_3a2[ispin_ref].array() - // -: conj(I) for orbital2
                                            filq1* orbital.array() * orbital2.conjugate().array();
                                    }
                                    else
                                    {
                                        // <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |*,q2,q1>
                                        qsum_temp_3a2[ispin_ref] = qsum_temp_3a2[ispin_ref].array() - // -: conj(I) for orbital2
                                            filq1 * orbital.array().abs2();
                                    }
                                    /*** [3a2] ***/

                                    if (ispin==ispinq) // parallel-spin only
                                    {
                                        /*** [3a4(,b2,b5 for TC)] ***/
                                        // <*,q1,*| \nabla_1 u_12 |q1,q2,*>
                                        qsum_temp_3a4[idim] = qsum_temp_3a4[idim].array() +
                                            filq1 * orbital.array() * phiq1.array();
                                        /*** [3a4(,b2,b5 for TC)] ***/

                                        /*** [3b2,b5 for BITC] ***/
                                        if (is_bitc)
                                        {
                                            // <*,q1,q2| \nabla_2 u_23 |*,*,q1> (= conj(qsum_temp_3a4) for TC) 
                                            qsum_temp_3b5[idim] = qsum_temp_3b5[idim].array() - // -: conj(I)
                                                filq1 * orbital2.conjugate().array() * chiq1_ref.conjugate().array();
                                        }
                                        /*** [3b2,b5 for BITC] ***/
                                    }
                                } // idim
                            } // ispin
                        } // ibandq1
                    } // isymq1
                } // iq1

                /*** [3a2] ***/
                for (int ispin=0; ispin<num_independent_spins; ispin++)
                {
                    MPI_Allreduce(qsum_temp_3a2[ispin].data(),
                                  qsum_temp_temp[ispin].data(),
                                  qsum_temp_temp[ispin].size(),
                                  MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
                    plane_wave_basis.FFT_forward(qsum_temp_temp[ispin], qsum_temp_temp[ispin]);
                }
                if (potentials.includes_div_correction()) 
                {
                    // \int d^3q/(2pi)^3 4pi*exp(-alpha q^2)/q^2 = 2/pi \int dq_0^infty exp(-alpha q^2)/q^2 = 1/sqrt(pi*alpha)
                    // Later divided by "num_kpoints * unit_cell_volume" (three_body_factor)
                    double div_corr =
                        1.0/std::sqrt(PI*potentials.alpha_Vaux()) * kpoints.num_kpoints() * crystal_structure.unit_cell_volume() 
                        - potentials.sum_of_Vaux_scf()[iq2];

                    const int nbands1_old = bloch_states.filling_old()[ispinq][iq2].size(); // Note! q1 = q2
                    for (int ibandq1=-nbands1_old; ibandq1<bloch_states.num_occupied_bands()[ispinq][iq2]; ibandq1++)
                    {
                        Eigen::VectorXcd &chiq1_ref = is_bitc ? chiq1 : phiq1; // bra orbital
                            
                        if (ibandq1>=0)
                        {
                            plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, 0,
                                                                 kpoints.is_time_reversal_used_at_k()[iq2][0],
                                                                 bloch_states.phik_scf()[ispinq][iq2][ibandq1][0], phiq1,
                                                                 "SCF");
                        }
                        else
                        {
                            plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, 0,
                                                                 kpoints.is_time_reversal_used_at_k()[iq2][0],
                                                                 bloch_states.phik_scf_old()[ispinq][iq2][-1-ibandq1][0], phiq1,
                                                                 "SCF");
                        }
                        plane_wave_basis.FFT_backward(phiq1, phiq1);
                            
                        if (is_bitc)
                        {
                            if (ibandq1>=0)
                            {
                                plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, 0,
                                                                     kpoints.is_time_reversal_used_at_k()[iq2][0],
                                                                     bloch_states.phik_left_scf()[ispinq][iq2][ibandq1][0], chiq1,
                                                                     "SCF");
                            }
                            else
                            {
                                plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, 0,
                                                                     kpoints.is_time_reversal_used_at_k()[iq2][0],
                                                                     bloch_states.phik_left_scf_old()[ispinq][iq2][-1-ibandq1][0], chiq1,
                                                                     "SCF");
                            }
                            plane_wave_basis.FFT_backward(chiq1, chiq1);
                        }
                            
                        // <q1|q2>
                        phiphi = chiq1_ref.conjugate().array() * phiq2.array();
                        plane_wave_basis.FFT_forward(phiphi, phiphi);

                        // <q2|q1> (unnecessary for bitc)
                        phij = chiq2_ref.conjugate().array() * phiq1.array();
                        plane_wave_basis.FFT_forward(phij, phij);

                        double filq1 = ibandq1>=0 ? 
                            bloch_states.filling()[ispinq][iq2][ibandq1] : bloch_states.filling_old()[ispinq][iq2][-1-ibandq1];

                        for (int ispin=0; ispin<2; ispin++)
                        {
                            int ispin_ref = num_independent_spins==1 ? 0 : ispin;
                            
                            qsum_temp_temp[ispin_ref](0) -= // const component (G=0)
                                FourPI * potentials.jastrow.A_long()[ispinq][ispin] *
                                (potentials.jastrow.A_long()[ispinq][ispin] * div_corr +
                                 potentials.jastrow.uk({0.0, 0.0, 0.0}, ispinq, ispin)) *
                                filq1 * phiphi(0) * phij(0);

                            double coefftmp = FourPI * potentials.jastrow.A_long()[ispinq][ispin] * filq1;
                            for (int ipw=1; ipw<plane_wave_basis.size_FFT_grid(); ipw++) // no G=0 term! (ipw!=0)
                            {
                                double uk = potentials.jastrow.uk(Gvect[ipw], ispinq, ispin);

                                qsum_temp_temp[ispin_ref](ipw) -= 
                                    coefftmp * uk * (phiphi(ipw)*phij(0) + phiphi(0)*phij(ipw));
                            } // ipw
                        } // ispin
                    } // ibandq1
                } // if (includes_div_correction)
                /*** [3a2] ***/

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
                                    std::conj(qsum_temp_temp[ispin](ipw));
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
                                    qsum_temp_temp[ispin](ipw);
                            }
                        } // ipw
                    } // if (time_reversal)

                    for (int ispin=0; ispin<num_independent_spins; ispin++)
                    {
                        V_3a2[ispin] = V_3a2[ispin].array() +
                            filq2 * qsum_temp_3a2[ispin].array();
                    }
                } // isymq2
                /*** [3a2] ***/

                /*** [3a4,b2,b5] ***/
                for (int idim=0; idim<3; idim++)
                {
                    MPI_Allreduce(qsum_temp_3a4[idim].data(),
                                  qsum_temp_temp[idim].data(),
                                  qsum_temp_temp[idim].size(),
                                  MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
                    if (is_bitc)
                    {
                        MPI_Allreduce(qsum_temp_3b5[idim].data(),
                                      qsum_temp_3b5_temp[idim].data(),
                                      qsum_temp_3b5_temp[idim].size(),
                                      MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
                    }
                    else
                    {
                        qsum_temp_3b5_temp[idim] = -qsum_temp_temp[idim].conjugate(); // -: conj(I)
                    }
                    plane_wave_basis.FFT_forward(qsum_temp_temp[idim], qsum_temp_temp[idim]);
                    plane_wave_basis.FFT_forward(qsum_temp_3b5_temp[idim], qsum_temp_3b5_temp[idim]);
                }

                // symmetry transformation of qsum_temp_3a4: (iq2, 0) -> (iq2, isymq2)
                // (k'(iq2,isymq2) = U^+ k(iq2,0))
                for (int isymq2=0; isymq2<kpoints.kvectors_scf()[iq2].size(); isymq2++)
                {
                    int ind_sym = kpoints.index_of_rotation_at_k()[iq2][isymq2];
                    Eigen::Matrix3i Udag = symmetry.rotation()[ind_sym].transpose(); // U^+
                    Eigen::Vector3d translation = symmetry.translation()[ind_sym];                    

                    Eigen::Matrix3d Udag_cart = // U^+ in the cartesian coordinate
                        (crystal_structure.reciprocal_vectors().transpose() * Udag.cast<double>() 
                         * crystal_structure.lattice_vectors()) / (2*PI);

                    Eigen::Vector3i Gvector, Gvector_trans;
                    for (int idim=0; idim<3; idim++)
                    {
                        qsum_temp_3a4[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                        qsum_temp_3b5[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
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
                            Complex phase = std::cos(2*PI*(translation.dot(Gvector_trans.cast<double>()))) +
                                I*std::sin(2*PI*(translation.dot(Gvector_trans.cast<double>()))); // (U^+ G) * r0

                            int ipw_trans = Gvector_trans(0)
                                + Gvector_trans(1) * size_FFT_grid_vec[0]
                                + Gvector_trans(2) * size_FFT_grid_vec[0] * size_FFT_grid_vec[1];
                            for (int idim=0; idim<3; idim++)
                            {
                                qsum_temp_3a4[idim](ipw_trans) -= phase *
                                    (Udag_cart(idim,0)*std::conj(qsum_temp_temp[0](ipw)) +
                                     Udag_cart(idim,1)*std::conj(qsum_temp_temp[1](ipw)) +
                                     Udag_cart(idim,2)*std::conj(qsum_temp_temp[2](ipw)));

                                qsum_temp_3b5[idim](ipw_trans) -= phase *
                                    (Udag_cart(idim,0)*std::conj(qsum_temp_3b5_temp[0](ipw)) +
                                     Udag_cart(idim,1)*std::conj(qsum_temp_3b5_temp[1](ipw)) +
                                     Udag_cart(idim,2)*std::conj(qsum_temp_3b5_temp[2](ipw)));
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
                            Complex phase = std::cos(2*PI*(translation.dot(Gvector_trans.cast<double>()))) +
                                I*std::sin(2*PI*(translation.dot(Gvector_trans.cast<double>()))); // (U^+ G) * r0

                            int ipw_trans = Gvector_trans(0)
                                + Gvector_trans(1) * size_FFT_grid_vec[0]
                                + Gvector_trans(2) * size_FFT_grid_vec[0] * size_FFT_grid_vec[1];
                            for (int idim=0; idim<3; idim++)
                            {
                                qsum_temp_3a4[idim](ipw_trans) += phase *
                                    (Udag_cart(idim,0)*qsum_temp_temp[0](ipw) +
                                     Udag_cart(idim,1)*qsum_temp_temp[1](ipw) +
                                     Udag_cart(idim,2)*qsum_temp_temp[2](ipw));

                                qsum_temp_3b5[idim](ipw_trans) += phase *
                                    (Udag_cart(idim,0)*qsum_temp_3b5_temp[0](ipw) +
                                     Udag_cart(idim,1)*qsum_temp_3b5_temp[1](ipw) +
                                     Udag_cart(idim,2)*qsum_temp_3b5_temp[2](ipw));
                            }
                        } // ipw
                    } // if (time_reversal)
                    for (int idim=0; idim<3; idim++)
                    {
                        plane_wave_basis.FFT_backward(qsum_temp_3a4[idim], qsum_temp_3a4[idim]);
                        plane_wave_basis.FFT_backward(qsum_temp_3b5[idim], qsum_temp_3b5[idim]);
                    }
                    /*** [3a4,b2,b5] ***/

                    q2symvect = crystal_structure.reciprocal_vectors().transpose() *
                        kpoints.kvectors_scf()[iq2][isymq2];

                    // recalculate phiq2 & chiq2 (including isymq2)
                    if (ibandq2>=0)
                    {
                        plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, isymq2,
                                                             kpoints.is_time_reversal_used_at_k()[iq2][isymq2],
                                                             bloch_states.phik_scf()[ispinq][iq2][ibandq2][0], phiq2,
                                                             "SCF");
                    }
                    else
                    {
                        plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, isymq2,
                                                             kpoints.is_time_reversal_used_at_k()[iq2][isymq2],
                                                             bloch_states.phik_scf_old()[ispinq][iq2][-1-ibandq2][0], phiq2,
                                                             "SCF");
                    }
                    plane_wave_basis.FFT_backward(phiq2, phiq2);
                
                    if (is_bitc)
                    {
                        if (ibandq2>=0)
                        {
                            plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, isymq2,
                                                                 kpoints.is_time_reversal_used_at_k()[iq2][isymq2],
                                                                 bloch_states.phik_left_scf()[ispinq][iq2][ibandq2][0], chiq2,
                                                                 "SCF");
                        }
                        else
                        {
                            plane_wave_basis.get_orbital_FFTgrid(ispinq, iq2, isymq2,
                                                                 kpoints.is_time_reversal_used_at_k()[iq2][isymq2],
                                                                 bloch_states.phik_left_scf_old()[ispinq][iq2][-1-ibandq2][0], chiq2,
                                                                 "SCF");
                        }
                        plane_wave_basis.FFT_backward(chiq2, chiq2);
                    }

                    /*** [3b2] ***/
                    for (int idim=0; idim<3; idim++)
                    {
                        qsum_temp_3b2[idim] = qsum_temp_3b2[idim].array() +
                            filq2 * qsum_temp_3b5[idim].array() * phiq2.array();
                    }
                    /*** [3b2] ***/

                    for (int ik=0; ik<num_irreducible_kpoints; ik++)
                    {
                        const Eigen::Vector3d kvector_ref = method.calc_mode()=="SCF" ? 
                            kpoints.kvectors_scf()[ik][0] : kpoints.kvectors_band()[ik][0];
                        const int num_G_at_k = method.calc_mode()=="SCF" ?
                            plane_wave_basis.num_G_at_k_scf()[ik] : plane_wave_basis.num_G_at_k_band()[ik];
                        const Eigen::VectorXi Gindex_at_k = method.calc_mode()=="SCF" ?
                            plane_wave_basis.Gindex_at_k_scf()[ispinq][ik][0] :  plane_wave_basis.Gindex_at_k_band()[ispinq][ik][0];
                        
                        bool is_k_assigned = false;
                        for (int jband=0; jband<H3phi[ispinq][ik].size(); jband++) // Note: spin of |j> = spin of |q1,q2> (spinq)
                        {
                            if (parallelization.is_assigned_irreducible_kpoints_all_bands()[ispinq][ik][jband]) 
                            {
                                is_k_assigned = true;
                                break;
                            }
                        }
                        if (!is_k_assigned) { continue; } // if no jband is assinged to this node

                        kvect = crystal_structure.reciprocal_vectors().transpose() * kvector_ref;
                        kqvect = kvect - q2symvect;

                        for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                        {
                            double uk = potentials.jastrow.uk(kqvect + Gvect[ipw], ispinq, ispinq);
                            for (int idim=0; idim<3; idim++)
                            {
                                du[0][idim](ipw) = (kqvect(idim) + Gvect[ipw](idim)) * uk;
                            }
                        }
                        for (int jband=0; jband<H3phi[ispinq][ik].size(); jband++)
                        {
                            if (!parallelization.is_assigned_irreducible_kpoints_all_bands()[ispinq][ik][jband]) { continue; }
                            for (int jspinor=0; jspinor<num_spinor; jspinor++)
                            {
                                // phi -> phij on the FFT-grid
                                plane_wave_basis.get_orbital_FFTgrid(ispinq, ik, 0, // isym = 0
                                                                     false, // time-rersal not used for isym=0
                                                                     phi[ispinq][ik][jband][jspinor], phij,
                                                                     method.calc_mode());
                                plane_wave_basis.FFT_backward(phij, phij); // -> phij(R)

                                orbital2 = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid()); // dummy of H3phi here

                                /*** [3a4] ***/
                                // <*,*,q2| |*,*,j>
                                phiphi = chiq2_ref.conjugate().array() * phij.array();
                                plane_wave_basis.FFT_forward(phiphi, phiphi);
                                phiphi *= filq2 * three_body_factor_a4;
                                for (int idim=0; idim<3; idim++)
                                {
                                    // <*,*,q2| \nabla_1 u_13 |*,*,j>
                                    orbital = phiphi.array() * du[0][idim].array();
                                    plane_wave_basis.FFT_backward(orbital, orbital); // idim-component

                                    // <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q1,q2,j>
                                    orbital2 = orbital2.array() +
                                        qsum_temp_3a4[idim].array() * orbital.array();
                                }
                                /*** [3a4] ***/

                                /*** [3b5] ***/
                                phiphi = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                                for (int idim=0; idim<3; idim++)
                                {
                                    // <*,q1,q2| \nabla_2 u_23 |*,j,q1>
                                    orbital = qsum_temp_3b5[idim].array() * phij.array();
                                    plane_wave_basis.FFT_forward(orbital, orbital); // idim-component

                                    // <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |*,j,q1>
                                    phiphi = phiphi.array() + orbital.array() * du[0][idim].array();
                                }
                                plane_wave_basis.FFT_backward(phiphi, phiphi);

                                // <*,q1,q2| | \nabla_2 u_21 \nabla_2 u_23 |q2,j,q1>
                                orbital2 = orbital2.array() +
                                    three_body_factor_b5 * filq2 *
                                    phiphi.array() * phiq2.array();
                                /*** [3b5] ***/

                                /*** [3a4,b5] ***/
                                plane_wave_basis.FFT_forward(orbital2, orbital2);
                                for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                                {
                                    H3phi[ispinq][ik][jband][0](ipw_at_k) +=  // Note: spin of |j> = spin of |q1,q2> (spinq)
                                        orbital2(Gindex_at_k(ipw_at_k));
                                }
                                /*** [3a4,b5] ***/
                            } // jspinor
                        } // jband
                    } // ik
                } // isymq2
            } // ibandq2
        } // iq2

        /*** [3b2] ***/
        for (int idim=0; idim<3; idim++)
        {
            plane_wave_basis.FFT_forward(qsum_temp_3b2[idim], qsum_temp_3b2[idim]);
        }
        for (int ispin=0; ispin<num_independent_spins; ispin++)
        {
            if (num_independent_spins==1)
            {
                for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                {
                    double uk = potentials.jastrow.uk(Gvect[ipw], ispin, 0) +
                        potentials.jastrow.uk(Gvect[ipw], ispin, 1); // summation over ispinq
                    for (int idim=0; idim<3; idim++)
                    {
                        du[0][idim](ipw) = Gvect[ipw](idim) * uk;
                    }
                }
            }
            else // num_independent_spins==2
            {
                for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                {
                    double uk = potentials.jastrow.uk(Gvect[ipw], ispin, ispinq);
                    for (int idim=0; idim<3; idim++)
                    {
                        du[0][idim](ipw) = Gvect[ipw](idim) * uk;
                    }
                }
            } // if (num_independent_spins)

            // <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |*,q2,q1>
            orbital = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
            for (int idim=0; idim<3; idim++)
            {
                orbital = orbital.array() +
                    du[0][idim].array() * qsum_temp_3b2[idim].array();
            }
            plane_wave_basis.FFT_backward(orbital, orbital);

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
                                                             false, // time-rersal not used for isym=0
                                                             phi[ispin][ik][jband][jspinor], phij,
                                                             method.calc_mode());
                        plane_wave_basis.FFT_backward(phij, phij); // -> phij(R)

                        // <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |j,q2,q1>
                        phiphi = orbital.array() * phij.array();
                        plane_wave_basis.FFT_forward(phiphi, phiphi);
                        for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                        {
                            H3phi[ispin][ik][jband][0](ipw_at_k) += 
                                three_body_factor_b2 * phiphi(Gindex_at_k(ipw_at_k));
                        }
                    } // jspinor
                } // jband
            } // ik
        } // ispin
        /*** [3b2] ***/
    } // ispinq

    /*** [3a2] ***/
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
                                                         false, // time-rersal not used for isym=0
                                                         phi[ispin][ik][jband][jspinor], phij,
                                                         method.calc_mode());
                    plane_wave_basis.FFT_backward(phij, phij); // -> phij(R)

                    orbital = V_3a2[ispin].array() * phij.array();
                    plane_wave_basis.FFT_forward(orbital, orbital);
                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                    {
                        H3phi[ispin][ik][jband][0](ipw_at_k) += 
                            three_body_factor_a2 *
                            orbital(Gindex_at_k(ipw_at_k));
                    }
                } // jspinor
            } // jband
        } // ik
    } // ispin
    /*** [3a2] ***/

    /*** [3a4,b5] ***/
    if (potentials.includes_div_correction()) 
    {
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
            
                kVaux = method.calc_mode()=="BAND" ?
                    -potentials.sum_of_kVaux_band()[ik] :
                    -potentials.sum_of_kVaux_scf()[ik];
                //  Here, "-1" comes from: \int \nabla Vaux - \sum \nabla Vaux = 0 - sum_of_kVaux = "-1" * sum_of_kVaux.
                kVaux *= potentials.jastrow.A_long()[ispin][ispin]; // since we consider not 1/r but A/r-like divergence

                if (kVaux.squaredNorm() > 1e-8) // if = 0 then no correction is needed
                {
                    const int nbands2_old = 
                        method.calc_mode()=="BAND" ? 0 : bloch_states.filling_old()[ispin][ik].size();
                    for (int ibandk2=-nbands2_old; ibandk2<num_bands_tc[ispin]; ibandk2++)
                    {
                        double filling_ibandk2 = 
                            ibandk2>=0 ? (method.calc_mode()=="BAND" ?
                                          kpoints.return_band_filling(spin, bloch_states.eigenvalues_band()[ispin][ik][ibandk2].real(),
                                                                      bloch_states.fermi_energy(),
                                                                      bloch_states.num_electrons(), ibandk2) :
                                          bloch_states.filling()[ispin][ik][ibandk2]) :
                            bloch_states.filling_old()[ispin][ik][-1-ibandk2];

                        if (ibandk2>=0 && method.calc_mode()=="SCF")
                        {
                            // above filling_ibandk2 (=filling[ispin][ik][ibandk2]) is not correct for zero-weight k-points in SCF calc.
                            if (kpoints.kweight_scf()[ik] < 1e-8)
                            {
                                filling_ibandk2 =
                                    kpoints.return_band_filling(spin, bloch_states.eigenvalues_scf()[ispin][ik][ibandk2].real(),
                                                                bloch_states.fermi_energy(),
                                                                bloch_states.num_electrons(), ibandk2);
                            }
                        }

                        if (filling_ibandk2 < 1e-8) { continue; } // no filling                                   

                        const Eigen::VectorXcd &phik2_ref =
                            ibandk2>=0 ? ( method.calc_mode()=="BAND" ?
                                           bloch_states.phik_band()[ispin][ik][ibandk2][0] :
                                           bloch_states.phik_scf()[ispin][ik][ibandk2][0]) :
                            bloch_states.phik_scf_old()[ispin][ik][-1-ibandk2][0];

                        const Eigen::VectorXcd &chik2_ref = !is_bitc ? phik2_ref :
                            (ibandk2>=0 ? (method.calc_mode()=="BAND" ? bloch_states.phik_left_band()[ispin][ik][ibandk2][0] :
                                           bloch_states.phik_left_scf()[ispin][ik][ibandk2][0]) :
                             bloch_states.phik_left_scf_old()[ispin][ik][-1-ibandk2][0]);

                        Eigen::VectorXcd &chiq2_ref = is_bitc ? chiq2 : phiq2;

                        plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, 
                                                             false,
                                                             phik2_ref, phiq2,
                                                             method.calc_mode());
                        plane_wave_basis.FFT_backward(phiq2, phiq2);

                        if (is_bitc)
                        {
                            plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, 
                                                                 false,
                                                                 chik2_ref, chiq2,
                                                                 method.calc_mode());
                            plane_wave_basis.FFT_backward(chiq2, chiq2);
                        }

                        orbital = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());  /*** [3a4] ***/
                        orbital2 = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());  /*** [3b5] ***/
                        for (int iq1=0; iq1<num_irreducible_kpoints_scf; iq1++)
                        {
                            for (int isymq1=0; isymq1<kpoints.kvectors_scf()[iq1].size(); isymq1++)
                            {
                                bool is_q1_assigned = false;
                                const int nbands1_old = bloch_states.filling_old()[ispin][iq1].size();
                                for (int ibandq1=-nbands1_old; ibandq1<bloch_states.num_occupied_bands()[ispin][iq1]; ibandq1++)
                                {
                                    if (ibandq1>=0)
                                    {
                                        if (parallelization.is_assigned_all_kpoints_occupied_bands()[ispin][iq1][isymq1][ibandq1]) 
                                        {
                                            is_q1_assigned = true;
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (parallelization.is_assigned_all_kpoints_occupied_bands_old()[ispin][iq1][isymq1][-1-ibandq1]) 
                                        {
                                            is_q1_assigned = true;
                                            break;
                                        }
                                    }
                                }
                                if (!is_q1_assigned) { continue; } // if no ibandq1 is assinged to this node

                                q1vect = crystal_structure.reciprocal_vectors().transpose() *
                                    kpoints.kvectors_scf()[iq1][isymq1];
                                kqvect = kvect - q1vect;

                                for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                                {
                                    double uk = potentials.jastrow.uk(kqvect + Gvect[ipw], ispin, ispin);
                                    du[0][0](ipw) = kVaux.dot(kqvect + Gvect[ipw]) * uk; // kVaux * du
                                }
                                for (int ibandq1=-nbands1_old; ibandq1<bloch_states.num_occupied_bands()[ispin][iq1]; ibandq1++)
                                {
                                    if (ibandq1>=0)
                                    {
                                        if (!parallelization.is_assigned_all_kpoints_occupied_bands()[ispin][iq1][isymq1][ibandq1]) { continue; }
                                    }
                                    else
                                    {
                                        if (!parallelization.is_assigned_all_kpoints_occupied_bands_old()[ispin][iq1][isymq1][-1-ibandq1]) { continue; }
                                    }

                                    Eigen::VectorXcd &chiq1_ref = is_bitc ? chiq1 : phiq1; // bra orbital
                                    
                                    if (ibandq1>=0)
                                    {
                                        plane_wave_basis.get_orbital_FFTgrid(ispin, iq1, isymq1,
                                                                             kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                             bloch_states.phik_scf()[ispin][iq1][ibandq1][0], phiq1,
                                                                             "SCF");
                                    }
                                    else
                                    {
                                        plane_wave_basis.get_orbital_FFTgrid(ispin, iq1, isymq1,
                                                                             kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                             bloch_states.phik_scf_old()[ispin][iq1][-1-ibandq1][0], phiq1,
                                                                             "SCF");
                                    }
                                    plane_wave_basis.FFT_backward(phiq1, phiq1);
                                    
                                    if (is_bitc)
                                    {
                                        if (ibandq1>=0)
                                        {
                                            plane_wave_basis.get_orbital_FFTgrid(ispin, iq1, isymq1,
                                                                                 kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                                 bloch_states.phik_left_scf()[ispin][iq1][ibandq1][0], chiq1,
                                                                                 "SCF");
                                        }
                                        else
                                        {
                                            plane_wave_basis.get_orbital_FFTgrid(ispin, iq1, isymq1,
                                                                                 kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                                 bloch_states.phik_left_scf_old()[ispin][iq1][-1-ibandq1][0], chiq1,
                                                                                 "SCF");
                                        }
                                        plane_wave_basis.FFT_backward(chiq1, chiq1);
                                    }

                                    phiphi = phiq2.array() * chiq1_ref.conjugate().array();
                                    plane_wave_basis.FFT_forward(phiphi, phiphi);

                                    phij = phiphi.array() * du[0][0].array(); // <*,q1| \nabla u |*,ibandq2>
                                    plane_wave_basis.FFT_backward(phij, phij);

                                    double filq1 = ibandq1>=0 ? bloch_states.filling()[ispin][iq1][ibandq1] : bloch_states.filling_old()[ispin][iq1][-1-ibandq1];

                                    // *** [3a4] ***
                                    orbital = orbital.array() +
                                        filq1 * phij.array() * phiq1.array(); // <*,q1| \nabla u |q1,ibandq2>
                                    // *** [3a4] ***

                                    // *** [3b5] ***
                                    if (is_bitc) // else, phiphi and phij = those calculated above
                                    {
                                        phiphi = chiq2_ref.array() * phiq1.conjugate().array();
                                        plane_wave_basis.FFT_forward(phiphi, phiphi);

                                        phij = phiphi.array() * du[0][0].array(); // conj(<*,ibandq2| \nabla u |*,q1>)
                                        plane_wave_basis.FFT_backward(phij, phij);
                                    }
                                    orbital2 = orbital2.array() - // -: conj(I)
                                        filq1 * phij.array() * chiq1_ref.array(); // conj(<q1,ibandq2| \nabla u |*,q1>)
                                    // *** [3b5] ***

                                } // ibandq1
                            } // isymq1
                        } // iq1
                        plane_wave_basis.FFT_forward(orbital, orbital); // *** [3a4] ***
                        plane_wave_basis.FFT_forward(orbital2, orbital2); // *** [3b5] ***

                        for (int jband=0; jband<H3phi[ispin][ik].size(); jband++)
                        {
                            for (int jspinor=0; jspinor<num_spinor; jspinor++)
                            {
                                // phi -> phij on the FFT-grid
                                plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, // isym = 0
                                                                     false, // time-rersal not used for isym=0
                                                                     phi[ispin][ik][jband][jspinor], phij,
                                                                     method.calc_mode());
                                // plane_wave_basis.FFT_backward(phij, phij); // -> phij(R)
                                Complex chik2_j = chik2_ref.dot(phi[ispin][ik][jband][jspinor]); // <ibandq2|j>

                                // *** [3a4] ***
                                Complex coeff_div_1_3a4 = three_body_factor_a4 * filling_ibandk2 * chik2_j;
                                Complex coeff_div_2_3a4 = three_body_factor_a4 * filling_ibandk2 *
                                    filling_ibandk2 * chik2_j * kVaux.dot(kVaux) / // filling_ibandk1 * delta_{ibandk1, ibandk2}
                                    static_cast<double>(parallelization.num_mpi_processes()); // since this term is not parallelized

                                for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                                {
                                    H3phi[ispin][ik][jband][0](ipw_at_k) += 
                                        coeff_div_1_3a4 * orbital(Gindex_at_k(ipw_at_k));
                                }
                                H3phi[ispin][ik][jband][0] += coeff_div_2_3a4 * phik2_ref;
                                // *** [3a4] ***

                                // *** [3b5] ***
                                Complex coeff_div_3b5 = orbital2.dot(phij) + // <q1,ibandq2| \nabla u |j,q1>, in G-space
                                    filling_ibandk2 * chik2_j * kVaux.dot(kVaux) /
                                    static_cast<double>(parallelization.num_mpi_processes()); // since this term is not parallelized

                                coeff_div_3b5 *= three_body_factor_b5 * filling_ibandk2;
                                H3phi[ispin][ik][jband][0] += coeff_div_3b5 * phik2_ref;
                            } // jspinor
                        } // jband
                    } // ibandk2
                } // if (kVaux!=0)
            } // ik
        } // ispin
    } // if (includes_div_correction)
}
