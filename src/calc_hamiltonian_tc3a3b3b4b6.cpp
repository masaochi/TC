// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

// [3a3] 0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q1,j,q2>
//  *equivalent then including: [3a6] 0.5*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q2,q1,j> (x2<->x3 && q1<->q2)
// [3b3] 0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |q1,j,q2>
// [3b4] -0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |q1,q2,j>
// [3b6] 0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |q2,q1,j>
void calc_hamiltonian::tc3a3b3b4b6(const Parallelization &parallelization, 
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

    double three_body_factor_base = 1.0/(kpoints.num_kpoints()*kpoints.num_kpoints()*
                                         crystal_structure.unit_cell_volume()*
                                         crystal_structure.unit_cell_volume());
    
    double three_body_factor_a3 = -three_body_factor_base; // (0.5*2) * I^2 = -1.0
    if (num_independent_spins==1) { three_body_factor_a3 /= 2.0; }  // spin indices of x1,x2 are the same

    double three_body_factor_b3 = three_body_factor_base; // (0.5*2) * I^2 * (-1) (convolution \nabla_2 u_21 = -\nabla_1 u_12) = 1.0
    if (num_independent_spins==1) { three_body_factor_b3 /= 2.0; }  // spin indices of x1,x2 are the same

    double three_body_factor_b4 = -three_body_factor_base; // (-0.5*2) * I^2 * (-1) (convolution \nabla_2 u_21 = -\nabla_1 u_12) = -1.0
    if (num_independent_spins==1) { three_body_factor_b4 /= 4.0; }  // spin indices of x1,x2,x3 are the same

    double three_body_factor_b6 = three_body_factor_base; // (0.5*2) * I^2 * (-1) (convolution \nabla_2 u_21 = -\nabla_1 u_12) = 1.0
    if (num_independent_spins==1) { three_body_factor_b6 /= 2.0; }  // spin indices of x1,x2 are the same

    // When you would like to debug each term:
    //three_body_factor_a3 = 0.0;
    //three_body_factor_b3 = 0.0;
    //three_body_factor_b4 = 0.0;
    //three_body_factor_b6 = 0.0;

    std::vector<Eigen::Vector3d> Gvect(plane_wave_basis.size_FFT_grid());
    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
    {
        Gvect[ipw] = crystal_structure.reciprocal_vectors().transpose()
            * plane_wave_basis.get_Gvector(ipw).cast<double>();
    }

    Eigen::Vector3d kqvect, kvect, qvect;
    Eigen::VectorXcd orbital(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd orbital2(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd phij(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd phiq(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd chiq;
    if (is_bitc) { chiq.resize(plane_wave_basis.size_FFT_grid()); }
    Eigen::VectorXcd phiphi(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd Hphij_sub(plane_wave_basis.size_FFT_grid());

    std::vector<std::vector<Eigen::VectorXd> > du(2); // spin
    for (int ispin=0; ispin<2; ispin++)
    {
        du[ispin].resize(3);
        for (int idim=0; idim<3; idim++)
        {
            du[ispin][idim] = Eigen::VectorXd::Zero(plane_wave_basis.size_FFT_grid());
        }
    }
    std::vector<Eigen::VectorXcd> dnu(3);
    for (int idim=0; idim<3; idim++)
    {
        dnu[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> dnu_phi(3);
    for (int idim=0; idim<3; idim++)
    {
        dnu_phi[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }
    std::vector<Eigen::VectorXcd> qsum_temp_3a3b4(3);
    for (int idim=0; idim<3; idim++)
    {
        qsum_temp_3a3b4[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }

    // divergence correction
    Eigen::Vector3d kVaux;
    std::vector<Eigen::VectorXd> AuG(2); // *** [3b6] ***
    if (potentials.includes_div_correction()) 
    {
        for (int jspin=0; jspin<2; jspin++)
        {
            AuG[jspin] = Eigen::VectorXd::Zero(plane_wave_basis.size_FFT_grid());
            for (int ipw=1; ipw<plane_wave_basis.size_FFT_grid(); ipw++) // Note: ipw=0 is excluded
            {
                AuG[jspin](ipw) = potentials.jastrow.uk(Gvect[ipw],0,jspin);
            }
            AuG[jspin] = AuG[jspin] * potentials.jastrow.A_long()[0][jspin];
        }
    }
    Eigen::VectorXcd temp_div_3b6_G(plane_wave_basis.size_FFT_grid()); // in G-space
    Eigen::VectorXcd temp_div_3b6_R(plane_wave_basis.size_FFT_grid()); // in R-space

    for (int ispin=0; ispin<num_independent_spins; ispin++) // = spin index of x1
    {
        /*** [3a3,b3] ***/
        for (int ispin3=0; ispin3<2; ispin3++) // spin index of x3                                                      
        {
            orbital = bloch_states.density()[ispin3]; // density(R). Name "orbital" is meaningless here
            plane_wave_basis.FFT_forward(orbital, orbital); // -> density(G)

            // [3a3] <*,*,q2| \nabla_1 u_13 |*,*,q2> [3b3] <*,*,q2| \nabla_2 u_23 |*,*,q2>
            for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
            {
                double uk = potentials.jastrow.uk(Gvect[ipw], ispin, ispin3);
                for (int idim=0; idim<3; idim++)
                {
                    dnu[idim](ipw) +=
                        uk * Gvect[ipw](idim) * orbital(ipw);
                }
            }
        }
        for (int idim=0; idim<3; idim++)
        {
            plane_wave_basis.FFT_backward(dnu[idim], dnu[idim]); // -> dnu(R)
        }
        /*** [3a3,b3] ***/

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

                    /*** [3a3,b4] ***/
                    for (int idim=0; idim<3; idim++)
                    {
                        qsum_temp_3a3b4[idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                    }
                    /*** [3a3,b4] ***/

                    /*** [3b3] ***/
                    // <*,*,q2| \nabla_2 u_23 |*,j,q2>
                    for (int idim=0; idim<3; idim++) 
                    {
                        dnu_phi[idim] = dnu[idim].array() * phij.array();
                    }
                    Hphij_sub = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                    /*** [3b3] ***/

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
                                double uk0 = potentials.jastrow.uk(kqvect + Gvect[ipw], 0, ispin); // [3b6] 0 or 1 = spin index of x2
                                double uk1 = potentials.jastrow.uk(kqvect + Gvect[ipw], 1, ispin);
                                for (int idim=0; idim<3; idim++)
                                {
                                    du[0][idim](ipw) = (kqvect(idim) + Gvect[ipw](idim)) * uk0;
                                    du[1][idim](ipw) = (kqvect(idim) + Gvect[ipw](idim)) * uk1;
                                }
                            }

                            for (int ibandq=0; ibandq<bloch_states.num_occupied_bands()[ispin][iq]; ibandq++)
                            {
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

                                phiphi = chiq_ref.conjugate().array() * phij.array();
                                plane_wave_basis.FFT_forward(phiphi, phiphi);

                                /*** [3b3] ***/
                                orbital2 = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                                for (int idim=0; idim<3; idim++)
                                {
                                    // <*,q1,q2| \nabla_2 u_23 |*,j,q2>
                                    orbital = dnu_phi[idim].array() * chiq_ref.conjugate().array();
                                    plane_wave_basis.FFT_forward(orbital, orbital); // idim-component

                                    // <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |*,j,q2>
                                    orbital2 = orbital2.array() + du[ispin][idim].array() * orbital.array();
                                }
                                plane_wave_basis.FFT_backward(orbital2, orbital2);

                                // <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |q1,j,q2>
                                Hphij_sub = Hphij_sub.array() +
                                    (bloch_states.filling()[ispin][iq][ibandq] * three_body_factor_b3) *
                                    orbital2.array() * phiq.array();
                                /*** [3b3] ***/

                                /*** [3a3,b4,b6] ***/
                                phiphi *= bloch_states.filling()[ispin][iq][ibandq];
                                orbital2 = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                                for (int ispin2=0; ispin2<2; ispin2++) // spin index of x2
                                {
                                    for (int idim=0; idim<3; idim++)
                                    {
                                        // [3a3] <*,q1,*| \nabla_1 u_12 |*,j,*>, [3b4,b6] <*,*,q2| \nabla_2 u_23 |*,*,j>
                                        orbital = phiphi.array() * du[ispin2][idim].array();
                                        plane_wave_basis.FFT_backward(orbital, orbital); // idim-component

                                        /*** [3a3,b4] ***/
                                        if (ispin2==ispin)
                                        {
                                            // [3a3] <*,q1,*| \nabla_1 u_12 |q1,j,*>, [3b4] <*,*,q2| \nabla_2 u_23 |*,q2,j>
                                            qsum_temp_3a3b4[idim] = qsum_temp_3a3b4[idim].array() +
                                                orbital.array() * phiq.array();
                                        }
                                        /*** [3a3,b4] ***/

                                        /*** [3b6] ***/
                                        // <*,q1,q2| \nabla_2 u_23 |*,q1,j>
                                        orbital = orbital.array() * bloch_states.density()[ispin2].array();
                                        plane_wave_basis.FFT_forward(orbital, orbital);

                                        // <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |*,q1,j>
                                        orbital2 = orbital2.array() +
                                            orbital.array() * du[ispin2][idim].array();
                                        /*** [3b6] ***/
                                    } // idim
                                } // ispin2
                                /*** [3a3,b4,b6] ***/

                                /*** [3b6] ***/
                                plane_wave_basis.FFT_backward(orbital2, orbital2);

                                // <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |q2,q1,j>
                                Hphij_sub = Hphij_sub.array() +
                                    three_body_factor_b6 * orbital2.array() * phiq.array();
                                /*** [3b6] ***/
                            } // ibandq
                        } // isymq
                    } // iq

                    /*** [3b6] ***/
                    if (potentials.includes_div_correction()) 
                    {
                        bool if_match = false; // whether ik is included in SCF k-mesh
                        for (int iq=0; iq<num_irreducible_kpoints_scf; iq++)
                        {
                            for (int isymq=0; isymq<kpoints.kvectors_scf()[iq].size(); isymq++)
                            {
                                if (kpoints.kweight_scf()[iq] < 1e-8) { continue; } // zero k-weight (i.e. band k-weight)
                                Eigen::Vector3d kq = kvector_ref - kpoints.kvectors_scf()[iq][isymq];
                                for (int idim=0; idim<3; idim++)
                                {
                                    while (kq(idim) < -1e-8) { kq(idim) += 1.0; }
                                    while (kq(idim) > 1.0-1e-8) { kq(idim) -= 1.0; }
                                }
                                if (kq.squaredNorm()<1e-8) { if_match = true; }
                            }
                        }

                        // \int d^3q/(2pi)^3 4pi*exp(-alpha q^2)/q^2 = 2/pi \int dq_0^infty exp(-alpha q^2)/q^2 = 1/sqrt(pi*alpha)
                        // Later divided by "num_kpoints * unit_cell_volume" (three_body_factor)
                        double div_corr = method.calc_mode()=="BAND" ?
                            1.0/std::sqrt(PI*potentials.alpha_Vaux()) * kpoints.num_kpoints() * crystal_structure.unit_cell_volume() 
                            - potentials.sum_of_Vaux_band()[ik] :
                            1.0/std::sqrt(PI*potentials.alpha_Vaux()) * kpoints.num_kpoints() * crystal_structure.unit_cell_volume() 
                            - potentials.sum_of_Vaux_scf()[ik];

                        // set kVaux (jastrow-A will multiplied later)
                        kVaux = method.calc_mode()=="BAND" ?
                            -potentials.sum_of_kVaux_band()[ik] :
                            -potentials.sum_of_kVaux_scf()[ik];
                        //  Here, "-1" comes from: \int \nabla Vaux - \sum \nabla Vaux = 0 - sum_of_kVaux = "-1" * sum_of_kVaux.

                        // set "temp_div_3b6"
                        Complex temp_div_3b6_1 = 0.0;
                        temp_div_3b6_G = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                        for (int jspin=0; jspin<2; jspin++)
                        {
                            orbital = bloch_states.density()[jspin];
                            plane_wave_basis.FFT_forward(orbital, orbital); // density(G)

                            double Adiv_corr = div_corr * potentials.jastrow.A_long()[ispin][jspin];
                            if (if_match) { Adiv_corr += potentials.jastrow.uk({0.0, 0.0, 0.0}, ispin, jspin); }
                            temp_div_3b6_1 += FourPI * Adiv_corr * potentials.jastrow.A_long()[ispin][jspin]
                                * orbital(0); // n(G=0)
                            int ispin_jspin = ispin==jspin ? 0 : 1; // parallel(0) or anti-parallel(1)
                            if (if_match)
                            {
                                temp_div_3b6_G = temp_div_3b6_G.array() +
                                    FourPI * orbital.array() * AuG[ispin_jspin].array(); // A*n(G)*u(G) where u(0)=0 here
                            }
                            else
                            {
                                for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
                                {
                                    temp_div_3b6_G(ipw) += kVaux.dot(Gvect[ipw]) * orbital(ipw) * AuG[ispin_jspin](ipw);
                                }
                            } // if_match
                        } // jspin
                        plane_wave_basis.FFT_backward(temp_div_3b6_G, temp_div_3b6_R);

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

                            Eigen::VectorXcd &chiq_ref = is_bitc ? chiq : phiq; // bra orbital

                            const Eigen::VectorXcd &phik_ref = method.calc_mode()=="BAND" ?
                                bloch_states.phik_band()[ispin][ik][ibandk][0] :
                                bloch_states.phik_scf()[ispin][ik][ibandk][0];
                            plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, // isym = 0
                                                                 ibandk, false, // time-rersal not used for isym=0
                                                                 phik_ref, phiq, method.calc_mode());
                            plane_wave_basis.FFT_backward(phiq, phiq);
                            if (is_bitc)
                            {
                                const Eigen::VectorXcd &chik_ref = method.calc_mode()=="BAND" ?
                                    bloch_states.phik_left_band()[ispin][ik][ibandk][0] :
                                    bloch_states.phik_left_scf()[ispin][ik][ibandk][0];
                                
                                plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, ibandk, false,
                                                                     chik_ref, chiq, method.calc_mode());
                                plane_wave_basis.FFT_backward(chiq, chiq);
                            }

                            // chiq * conj(phij). Note! conj is applied to phij (not chiq) here.
                            phiphi = chiq_ref.array() * phij.conjugate().array(); 
                            plane_wave_basis.FFT_forward(phiphi, phiphi);
                            Complex chiq_phij_integral = std::conj(phiphi(0)); // G=0 component is <q|j>

                            // contribution 1
                            H3phi[ispin][ik][jband][0] = H3phi[ispin][ik][jband][0] +
                                (filling_ibandk * three_body_factor_b6 * temp_div_3b6_1 * 
                                 chiq_phij_integral) * phik_ref;
                            int sign = if_match ? 1: -1;
                            // contribution 3 or 5 (sign: originating from G -> -G)
                            H3phi[ispin][ik][jband][0] = H3phi[ispin][ik][jband][0] +
                                (sign * filling_ibandk * three_body_factor_b6 * phiphi.dot(temp_div_3b6_G)) *
                                phik_ref;

                            // contribution 2 or 4 (in r-space thus Hphij_sub)
                            Hphij_sub = Hphij_sub.array() +
                                (filling_ibandk * chiq_phij_integral * three_body_factor_b6) *
                                temp_div_3b6_R.array() * phiq.array();                            
                        } // ibandk
                    } // if (includes_div_correction)
                    /*** [3b6] ***/

                    /*** [3b3,b6] ***/
                    plane_wave_basis.FFT_forward(Hphij_sub, Hphij_sub);
                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                    {
                        H3phi[ispin][ik][jband][0](ipw_at_k) += 
                            Hphij_sub(Gindex_at_k(ipw_at_k));
                    }
                    /*** [3b3,b6] ***/

                    /*** [3a3,b3,b4] ***/
                    if (potentials.includes_div_correction()) 
                    {
                        kVaux = method.calc_mode()=="BAND" ?
                            -potentials.sum_of_kVaux_band()[ik] :
                            -potentials.sum_of_kVaux_scf()[ik] ;
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

                                const Eigen::VectorXcd &chik_ref = !is_bitc ? phik_ref :
                                    (method.calc_mode()=="BAND" ? bloch_states.phik_left_band()[ispin][ik][ibandk][0] :
                                     bloch_states.phik_left_scf()[ispin][ik][ibandk][0]);

                                plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, // isym = 0
                                                                     ibandk, false, // time-rersal not used for isym=0
                                                                     phik_ref, phiq, method.calc_mode());
                                plane_wave_basis.FFT_backward(phiq, phiq);

                                if (is_bitc)
                                {
                                    plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, ibandk, false,
                                                                         chik_ref, chiq, method.calc_mode());
                                    plane_wave_basis.FFT_backward(chiq, chiq);
                                }
                                Eigen::VectorXcd &chiq_ref = is_bitc ? chiq : phiq; // bra orbital

                                /*** [3b3] ***/
                                Eigen::Vector3cd coeff_vec_3b3 = {0.0, 0.0, 0.0};
                                for (int idim=0; idim<3; idim++)
                                {
                                    coeff_vec_3b3[idim] = chiq_ref.dot(dnu_phi[idim]); // conj(chiq)*dnu_phi
                                }
                                coeff_vec_3b3 *= filling_ibandk * three_body_factor_b3 /
                                    plane_wave_basis.size_FFT_grid(); // integral in r-space: normalized const.

                                Complex temp_div_3b3 = 0.0;
                                for (int idim=0; idim<3; idim++)
                                {
                                    temp_div_3b3 += coeff_vec_3b3(idim) * kVaux(idim);
                                }
                                H3phi[ispin][ik][jband][0] += temp_div_3b3 * phik_ref;
                                /*** [3b3] ***/

                                /*** [3a3,b4] ***/
                                Eigen::Vector3cd coeff_vec_3a3b4 =
                                    (chik_ref.dot(phi[ispin][ik][jband][jspinor]) * filling_ibandk) * kVaux;
                                for (int idim=0; idim<3; idim++)
                                {
                                    qsum_temp_3a3b4[idim] = qsum_temp_3a3b4[idim].array() +
                                        coeff_vec_3a3b4(idim) * phiq.array();
                                }
                                /*** [3a3,b4] ***/
                            } // ibandk
                        } // if (kVaux != 0)
                    } // if (includes_div_correction)
                    /*** [3a3,b3] ***/

                    /*** [3a3] ***/
                    // <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q1,j,q2>
                    orbital = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid()); // "orbital" is a meaningless name here
                    for (int idim=0; idim<3; idim++)
                    {
                        orbital = orbital.array() + qsum_temp_3a3b4[idim].array() * dnu[idim].array();
                    }
                    plane_wave_basis.FFT_forward(orbital, orbital);
                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                    {
                        H3phi[ispin][ik][jband][0](ipw_at_k) += 
                            three_body_factor_a3 * orbital(Gindex_at_k(ipw_at_k));
                    }
                    /*** [3a3] ***/

                    /*** [3b4] ***/
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
                                    du[0][idim](ipw) = (kqvect(idim) + Gvect[ipw](idim)) * uk;
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

                                phiphi = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
                                for (int idim=0; idim<3; idim++)
                                {
                                    // < *,q1,q2| \nabla_2 u_23 |*,q2,j>
                                    orbital = chiq_ref.conjugate().array() * qsum_temp_3a3b4[idim].array();
                                    plane_wave_basis.FFT_forward(orbital, orbital); // idim-component

                                    // <*,q1,q2| \nablba_2 u_21 \nabla_2 u_23 |*,q2,j>
                                    phiphi = phiphi.array() + orbital.array() * du[0][idim].array();
                                }
                                plane_wave_basis.FFT_backward(phiphi, phiphi);

                                // <*,q1,q2| | \nabla_2 u_21 \nabla_2 u_23 |q1,q2,j>
                                Hphij_sub = Hphij_sub.array() +
                                    bloch_states.filling()[ispin][iq][ibandq] * phiphi.array() * phiq.array();
                            } // ibandq
                        } // isymq
                    } // iq

                    plane_wave_basis.FFT_forward(Hphij_sub, Hphij_sub);
                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                    {
                        H3phi[ispin][ik][jband][0](ipw_at_k) += 
                            three_body_factor_b4 * Hphij_sub(Gindex_at_k(ipw_at_k));
                    }

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

                                Eigen::Vector3cd coeff_vec_3b4 = {0.0, 0.0, 0.0};
                                for (int idim=0; idim<3; idim++)
                                {
                                    coeff_vec_3b4[idim] = phiq.dot(qsum_temp_3a3b4[idim]); // conj(phiq)*qsum_temp_3a3b4
                                }
                                coeff_vec_3b4 *= three_body_factor_b4 * filling_ibandk
                                    /plane_wave_basis.size_FFT_grid(); // integral in r-space: normalized const.

                                Complex temp_div_3b4 = 0.0;
                                for (int idim=0; idim<3; idim++)
                                {
                                    temp_div_3b4 += coeff_vec_3b4(idim) * kVaux(idim);
                                }
                                H3phi[ispin][ik][jband][0] += temp_div_3b4 * phik_ref;
                            } // ibandk
                        } // if (kVaux!=0)
                    } // if (includes_div_correction)
                    /*** [3b4] ***/

                } // jspinor
            } // jband
        } // ik
    } // ispin
}
