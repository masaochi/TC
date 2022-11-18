// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

// *** tc3a4 ***
// calculate -0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q1,q2,j>
//   [equivalent then including: 3a5] -0.5*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q2,j,q1>
void calc_hamiltonian::tc3a4(const Parallelization &parallelization, 
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

    // (-0.5*2) * I^2 = +1.0
    double three_body_factor = 1.0/(kpoints.num_kpoints()*kpoints.num_kpoints()*
                                    crystal_structure.unit_cell_volume()*
                                    crystal_structure.unit_cell_volume());
    // spin indices of x1,x2,x3 should be the same: filling should be divided with 4
    if (num_independent_spins==1) { three_body_factor /= 4.0; } 

    std::vector<Eigen::Vector3d> Gvect(plane_wave_basis.size_FFT_grid());
    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
    {
        Gvect[ipw] = crystal_structure.reciprocal_vectors().transpose()
            * plane_wave_basis.get_Gvector(ipw).cast<double>();
    }

    Eigen::Vector3d kqvect, kvect, q1vect, q2vect, q2symvect, q2q1vect, kq1vect;
    Eigen::VectorXcd orbital(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXd orbitalr(plane_wave_basis.size_FFT_grid()); // r = real
    Eigen::VectorXcd phij(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd phiq1(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd phiq2(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd chiq1;
    if (is_bitc) { chiq1.resize(plane_wave_basis.size_FFT_grid()); }
    Eigen::VectorXcd phiphi(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd Hphij_sub(plane_wave_basis.size_FFT_grid());

    std::vector<Eigen::VectorXd> du(3);
    for (int idim=0; idim<3; idim++)
    {
        du[idim] = Eigen::VectorXd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> qsum_temp_3a4(3);
    for (int idim=0; idim<3; idim++)
    {
        qsum_temp_3a4[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> qsum_temp_3a4_temp(3);
    for (int idim=0; idim<3; idim++)
    {
        qsum_temp_3a4_temp[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }

    // divergence correction
    Eigen::Vector3d kVaux;
    Eigen::VectorXcd temp_div_3a4;

    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int iq2=0; iq2<num_irreducible_kpoints_scf; iq2++)
        {
            q2vect = crystal_structure.reciprocal_vectors().transpose() *
                kpoints.kvectors_scf()[iq2][0];
            for (int ibandq2=0; ibandq2<bloch_states.num_occupied_bands()[ispin][iq2]; ibandq2++)
            {
                plane_wave_basis.get_orbital_FFTgrid(ispin, iq2, 0, ibandq2,
                                                     kpoints.is_time_reversal_used_at_k()[iq2][0],
                                                     bloch_states.phik_scf()[ispin][iq2][ibandq2][0], phiq2,
                                                     "SCF");
                plane_wave_basis.FFT_backward(phiq2, phiq2);

                for (int idim=0; idim<3; idim++)
                {
                    qsum_temp_3a4[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                }
                for (int iq1=0; iq1<num_irreducible_kpoints_scf; iq1++)
                {
                    for (int isymq1=0; isymq1<kpoints.kvectors_scf()[iq1].size(); isymq1++)
                    {
                        bool is_q1_assigned = false;
                        for (int ibandq1=0; ibandq1<bloch_states.num_occupied_bands()[ispin][iq1]; ibandq1++)
                        {
                            if (parallelization.is_assigned_all_kpoints_occupied_bands()[ispin][iq1][isymq1][ibandq1]) 
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
                            double uk = potentials.jastrow.uk(q2q1vect + Gvect[ipw], ispin, ispin);
                            for (int idim=0; idim<3; idim++)
                            {
                                du[idim](ipw) = (q2q1vect(idim) + Gvect[ipw](idim)) * uk;
                            }
                        }
                        for (int ibandq1=0; ibandq1<bloch_states.num_occupied_bands()[ispin][iq1]; ibandq1++)
                        {
                            if (!parallelization.is_assigned_all_kpoints_occupied_bands()[ispin][iq1][isymq1][ibandq1]) { continue; }
                            
                            Eigen::VectorXcd &chiq1_ref = is_bitc ? chiq1 : phiq1; // bra orbital
                            
                            plane_wave_basis.get_orbital_FFTgrid(ispin, iq1, isymq1, ibandq1,
                                                                 kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                 bloch_states.phik_scf()[ispin][iq1][ibandq1][0], phiq1,
                                                                 "SCF");
                            plane_wave_basis.FFT_backward(phiq1, phiq1);
                            
                            if (is_bitc)
                            {
                                plane_wave_basis.get_orbital_FFTgrid(ispin, iq1, isymq1, ibandq1,
                                                                     kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                     bloch_states.phik_left_scf()[ispin][iq1][ibandq1][0], chiq1,
                                                                     "SCF");
                                plane_wave_basis.FFT_backward(chiq1, chiq1);
                            }
                            
                            // <*,q1,*| |*,q2,*>
                            phiphi = chiq1_ref.conjugate().array() * phiq2.array();
                            plane_wave_basis.FFT_forward(phiphi, phiphi);
                            
                            phiphi *=  bloch_states.filling()[ispin][iq1][ibandq1];
                            for (int idim=0; idim<3; idim++)
                            {
                                // <*,q1,*| \nabla_1 u_12 |*,q2,*>
                                orbital = du[idim].array() * phiphi.array();
                                plane_wave_basis.FFT_backward(orbital, orbital); // idim-component
                                
                                // <*,q1,*| \nabla_1 u_12 |q1,q2,*>
                                qsum_temp_3a4[idim] = qsum_temp_3a4[idim].array() +
                                    orbital.array() * phiq1.array();
                            }
                        } // ibandq1
                    } // isymq1
                } // iq1

                // Allreduce!
                for (int idim=0; idim<3; idim++)
                {
                    MPI_Allreduce(qsum_temp_3a4[idim].data(),
                                  qsum_temp_3a4_temp[idim].data(),
                                  qsum_temp_3a4_temp[idim].size(),
                                  MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
                    plane_wave_basis.FFT_forward(qsum_temp_3a4_temp[idim], qsum_temp_3a4_temp[idim]);
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
                                    (Udag_cart(idim,0)*std::conj(qsum_temp_3a4_temp[0](ipw)) +
                                     Udag_cart(idim,1)*std::conj(qsum_temp_3a4_temp[1](ipw)) +
                                     Udag_cart(idim,2)*std::conj(qsum_temp_3a4_temp[2](ipw)));
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
                                    (Udag_cart(idim,0)*qsum_temp_3a4_temp[0](ipw) +
                                     Udag_cart(idim,1)*qsum_temp_3a4_temp[1](ipw) +
                                     Udag_cart(idim,2)*qsum_temp_3a4_temp[2](ipw));
                            }
                        } // ipw
                    } // if (time_reversal)
                    for (int idim=0; idim<3; idim++)
                    {
                        plane_wave_basis.FFT_backward(qsum_temp_3a4[idim], qsum_temp_3a4[idim]);
                    }

                    q2symvect = crystal_structure.reciprocal_vectors().transpose() *
                        kpoints.kvectors_scf()[iq2][isymq2];

                    // calculate chiq2
                    if (is_bitc)
                    {
                        plane_wave_basis.get_orbital_FFTgrid(ispin, iq2, isymq2, ibandq2,
                                                             kpoints.is_time_reversal_used_at_k()[iq2][isymq2],
                                                             bloch_states.phik_left_scf()[ispin][iq2][ibandq2][0], phiq2,
                                                             "SCF");
                    }
                    else
                    {
                        plane_wave_basis.get_orbital_FFTgrid(ispin, iq2, isymq2, ibandq2,
                                                             kpoints.is_time_reversal_used_at_k()[iq2][isymq2],
                                                             bloch_states.phik_scf()[ispin][iq2][ibandq2][0], phiq2,
                                                             "SCF");
                    }
                    plane_wave_basis.FFT_backward(phiq2, phiq2);

                    for (int ik=0; ik<num_irreducible_kpoints; ik++)
                    {
                        const Eigen::Vector3d kvector_ref = method.calc_mode()=="SCF" ? 
                            kpoints.kvectors_scf()[ik][0] : kpoints.kvectors_band()[ik][0];
                        const int num_G_at_k = method.calc_mode()=="SCF" ?
                            plane_wave_basis.num_G_at_k_scf()[ik] : plane_wave_basis.num_G_at_k_band()[ik];
                        const Eigen::VectorXi Gindex_at_k = method.calc_mode()=="SCF" ?
                            plane_wave_basis.Gindex_at_k_scf()[ispin][ik][0] :  plane_wave_basis.Gindex_at_k_band()[ispin][ik][0];
                        
                        bool is_k_assigned = false;
                        for (int jband=0; jband<H3phi[ispin][ik].size(); jband++)
                        {
                            if (parallelization.is_assigned_irreducible_kpoints_all_bands()[ispin][ik][jband]) 
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
                            double uk = potentials.jastrow.uk(kqvect + Gvect[ipw], ispin, ispin);
                            for (int idim=0; idim<3; idim++)
                            {
                                du[idim](ipw) = (kqvect(idim) + Gvect[ipw](idim)) * uk;
                            }
                        }
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

                                // <*,*,q2| |*,*,j>
                                phiphi = phiq2.conjugate().array() * phij.array();
                                plane_wave_basis.FFT_forward(phiphi, phiphi);
                                phiphi *= bloch_states.filling()[ispin][iq2][ibandq2] * three_body_factor;

                                Hphij_sub = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                                for (int idim=0; idim<3; idim++)
                                {
                                    // <*,*,q2| \nabla_1 u_13 |*,*,j>
                                    orbital = phiphi.array() * du[idim].array();
                                    plane_wave_basis.FFT_backward(orbital, orbital); // idim-component

                                    // <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q1,q2,j>
                                    Hphij_sub = Hphij_sub.array() +
                                        qsum_temp_3a4[idim].array() * orbital.array();
                                }
                                plane_wave_basis.FFT_forward(Hphij_sub, Hphij_sub);
                                for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                                {
                                    H3phi[ispin][ik][jband][0](ipw_at_k) += 
                                        Hphij_sub(Gindex_at_k(ipw_at_k));
                                }
                            } // jspinor
                        } // jband
                    } // ik
                } // isymq2
            } // ibandq2
        } // iq2
    } // ispin

    // divergence correction for \nabla u
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
                    for (int ibandk2=0; ibandk2<num_bands_tc[ispin]; ibandk2++)
                    {
                        double filling_ibandk2 = method.calc_mode()=="BAND" ?
                            kpoints.return_band_filling(spin, bloch_states.eigenvalues_band()[ispin][ik][ibandk2].real(),
                                                        bloch_states.fermi_energy(),
                                                        bloch_states.num_electrons(), ibandk2) :
                            kpoints.return_band_filling(spin, bloch_states.eigenvalues_scf()[ispin][ik][ibandk2].real(),
                                                        bloch_states.fermi_energy(),
                                                        bloch_states.num_electrons(), ibandk2);
                        if (filling_ibandk2 < 1e-8) { continue; } // no filling                                   
                        // NOTE! filling_ibandk2 is non-zero even for zero-weight k-points (in fake scf) 

                        const Eigen::VectorXcd &phik2_ref = method.calc_mode()=="BAND" ?
                            bloch_states.phik_band()[ispin][ik][ibandk2][0] :
                            bloch_states.phik_scf()[ispin][ik][ibandk2][0];

                        const Eigen::VectorXcd &chik2_ref = !is_bitc ? phik2_ref :
                            (method.calc_mode()=="BAND" ? bloch_states.phik_left_band()[ispin][ik][ibandk2][0] :
                             bloch_states.phik_left_scf()[ispin][ik][ibandk2][0]);

                        plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, 
                                                             ibandk2, false,
                                                             phik2_ref, phiq2,
                                                             method.calc_mode());
                        plane_wave_basis.FFT_backward(phiq2, phiq2);

                        temp_div_3a4 = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                        for (int iq1=0; iq1<num_irreducible_kpoints_scf; iq1++)
                        {
                            for (int isymq1=0; isymq1<kpoints.kvectors_scf()[iq1].size(); isymq1++)
                            {
                                bool is_q1_assigned = false;
                                for (int ibandq1=0; ibandq1<bloch_states.num_occupied_bands()[ispin][iq1]; ibandq1++)
                                {
                                    if (parallelization.is_assigned_all_kpoints_occupied_bands()[ispin][iq1][isymq1][ibandq1]) 
                                    {
                                        is_q1_assigned = true;
                                        break;
                                    }
                                }
                                if (!is_q1_assigned) { continue; } // if no ibandq1 is assinged to this node

                                q1vect = crystal_structure.reciprocal_vectors().transpose() *
                                    kpoints.kvectors_scf()[iq1][isymq1];
                                kq1vect = kvect - q1vect;

                                for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                                {
                                    double uk = potentials.jastrow.uk(kq1vect + Gvect[ipw], ispin, ispin);
                                    orbitalr(ipw) = kVaux.dot(kq1vect + Gvect[ipw]) * uk; // kVaux * du
                                }
                                for (int ibandq1=0; ibandq1<bloch_states.num_occupied_bands()[ispin][iq1]; ibandq1++)
                                {
                                    if (!parallelization.is_assigned_all_kpoints_occupied_bands()[ispin][iq1][isymq1][ibandq1]) { continue; }

                                    Eigen::VectorXcd &chiq1_ref = is_bitc ? chiq1 : phiq1; // bra orbital
                                    
                                    plane_wave_basis.get_orbital_FFTgrid(ispin, iq1, isymq1, ibandq1,
                                                                         kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                         bloch_states.phik_scf()[ispin][iq1][ibandq1][0], phiq1,
                                                                         "SCF");
                                    plane_wave_basis.FFT_backward(phiq1, phiq1);
                                    
                                    if (is_bitc)
                                    {
                                        plane_wave_basis.get_orbital_FFTgrid(ispin, iq1, isymq1, ibandq1,
                                                                             kpoints.is_time_reversal_used_at_k()[iq1][isymq1],
                                                                             bloch_states.phik_left_scf()[ispin][iq1][ibandq1][0], chiq1,
                                                                             "SCF");
                                        plane_wave_basis.FFT_backward(chiq1, chiq1);
                                    }

                                    phiphi = phiq2.array() * chiq1_ref.conjugate().array();
                                    plane_wave_basis.FFT_forward(phiphi, phiphi);

                                    phij = phiphi.array() * orbitalr.array(); // <*,q1| \nabla u |*,ibandq2>
                                    plane_wave_basis.FFT_backward(phij, phij);

                                    temp_div_3a4 = temp_div_3a4.array() +
                                        bloch_states.filling()[ispin][iq1][ibandq1] *
                                        phij.array() * phiq1.array(); // <*,q1| \nabla u |q1,ibandq2>
                                } // ibandq1
                            } // isymq1
                        } // iq1
                        plane_wave_basis.FFT_forward(temp_div_3a4, temp_div_3a4);

                        for (int jband=0; jband<H3phi[ispin][ik].size(); jband++)
                        {
                            for (int jspinor=0; jspinor<num_spinor; jspinor++)
                            {
                                Complex chik2_j = chik2_ref.dot(phi[ispin][ik][jband][jspinor]); // <ibandq2|j>

                                for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                                {
                                    H3phi[ispin][ik][jband][0](ipw_at_k) += 
                                        (three_body_factor * filling_ibandk2 * 
                                         chik2_j) * temp_div_3a4(Gindex_at_k(ipw_at_k));
                                }

                                // Not a double-count of filling_ibandk2. filling_ibandk1 * delta_{ibandk1, ibandk2}
                                H3phi[ispin][ik][jband][0] += 
                                    (three_body_factor * filling_ibandk2 *
                                     filling_ibandk2 * chik2_j *
                                     kVaux.dot(kVaux)/static_cast<double>(parallelization.num_mpi_processes()))* // since this term is not parallelized
                                    phik2_ref;
                            } // jspinor
                        } // jband
                    } // ibandk2
                } // if (kVaux!=0)
            } // ik
        } // ispin
    } // if (includes_div_correction)
}
