// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

void calc_hamiltonian::pseudo(const Parallelization &parallelization, 
                              const Method &method,
                              const CrystalStructure &crystal_structure,
                              const Potentials &potentials,
                              const Spin &spin, const Kpoints &kpoints,
                              PlaneWaveBasis &plane_wave_basis, 
                              const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
                              std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H1phi,
                              std::ostream *ost)
{
    const int num_independent_spins = spin.num_independent_spins();
    const int num_irreducible_kpoints = method.calc_mode()=="SCF" ? 
        kpoints.num_irreducible_kpoints_scf() : kpoints.num_irreducible_kpoints_band();
    const int num_spinor = (spin.is_spinor()==false) ? 1 : 2;

    Eigen::Vector3d kGvect;
    Eigen::VectorXcd orbital(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd V_orbital(plane_wave_basis.size_FFT_grid());

    Eigen::VectorXcd phase_atom;
    std::vector<std::vector<Eigen::VectorXcd> > ylm(4); // s, p, d, f
    for (int il=0; il<ylm.size(); il++) { ylm[il].resize(il+1); } // we calculate only m = 0, 1,..., l (negative m can be obtained by positive m)
    std::vector<std::vector<double> > ylm_prefactor(ylm.size());
    for (int il=0; il<ylm.size(); il++)
    {
        ylm_prefactor[il].resize(ylm[il].size());
        for (int im=0; im<ylm[il].size(); im++)
        {
            int ifac1 = 1; // (il-im)!
            for (int i=1; i<=(il-im); i++) { ifac1 *= i; }
            int ifac2 = 1; // (il+im)!
            for (int i=1; i<=(il+im); i++) { ifac2 *= i; }

            ylm_prefactor[il][im] = std::sqrt(static_cast<double>((2*il+1)*ifac1)/(FourPI*ifac2));
        }
    }
    std::vector<std::vector<std::vector<Eigen::VectorXcd> > > projector(crystal_structure.num_atomic_species());

    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            const Eigen::Vector3d kvector = method.calc_mode()=="SCF" ? 
                kpoints.kvectors_scf()[ik][0] : kpoints.kvectors_band()[ik][0];
            const int num_G_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.num_G_at_k_scf()[ik] : plane_wave_basis.num_G_at_k_band()[ik];
            const Eigen::VectorXi Gindex_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.Gindex_at_k_scf()[ispin][ik][0] :  plane_wave_basis.Gindex_at_k_band()[ispin][ik][0];

            bool projector_calculated = false;
            for (int jband=0; jband<H1phi[ispin][ik].size(); jband++)
            {
                if (!parallelization.is_assigned_irreducible_kpoints_all_bands()[ispin][ik][jband]) { continue; }
                if (!projector_calculated) // projector (non-local pseudopot.) is independent of the band index
                {
                    // calculate Y_lm(e_{k+G})
                    for (int il=0; il<ylm.size(); il++) 
                    {
                        for (int im=0; im<ylm[il].size(); im++)
                        {
                            ylm[il][im].resize(num_G_at_k);
                        }
                    }
                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++) 
                    {
                        kGvect = crystal_structure.reciprocal_vectors().transpose()
                            *(kvector + plane_wave_basis.get_Gvector(Gindex_at_k(ipw_at_k)).cast<double>());
                        double kG = kGvect.norm();
                        
                        if (kG < 1e-8) // we can set ylm=0 for l>=1 since j_l=0 (l>=1) in projector[][]
                                       // also see a small note for potentials.projector_nonlocal(). (potentials.cpp)
                        {
                            ylm[0][0](ipw_at_k) = 1.0/std::sqrt(FourPI);
                            for (int il=1; il<ylm.size(); il++)
                            {
                                for (int im=0; im<ylm[il].size(); im++)
                                {
                                    ylm[il][im](ipw_at_k) = 0.0;
                                }
                            }
                        }
                        else
                        {
                            double cos_theta = kGvect(2) / kG; // z/r
                            double rsin_theta = std::sqrt(kGvect(0)*kGvect(0) + kGvect(1)*kGvect(1)); // r sin(theta)
                            Complex phase = 1;
                            if (rsin_theta > 1e-8) { phase = (kGvect(0) + I*kGvect(1)) / rsin_theta; } // No phase required for rsin_theta=0. (P_lm(+-1)=0 for m!=0)
                            for (int il=0; il<ylm.size(); il++) 
                            {
                                for (int im=0; im<ylm[il].size(); im++)
                                {
                                    ylm[il][im](ipw_at_k) = ylm_prefactor[il][im] 
                                        * boost::math::legendre_p(il, im, cos_theta);
                                    for (int i=0; i<im; i++) { ylm[il][im](ipw_at_k) *= phase; }
                                }
                            }
                        } // if (kG<1e-8)
                    } // ipw_at_k
                    // calculate <beta|k+G> (except the exp[i(k+G)R_atom] factor)
                    for (int iatomic_species=0; iatomic_species<crystal_structure.num_atomic_species(); iatomic_species++)
                    {
                        int num_projectors = potentials.pseudo_lbeta()[iatomic_species].size();
                        projector[iatomic_species].resize(num_projectors);
                        for (int iproj=0; iproj<num_projectors; iproj++)
                        {
                            projector[iatomic_species][iproj] 
                                = potentials.projector_nonlocal(crystal_structure, plane_wave_basis,
                                                                Gindex_at_k,
                                                                kvector, ylm,
                                                                ispin, ik, iatomic_species, iproj);
                        }
                    }
                    projector_calculated = true;
                } // if (!projector_calculated)

                for (int jspinor=0; jspinor<num_spinor; jspinor++)
                {
                    // (1/2) pseudo potential (non-local part)
                    for (int iatom=0; iatom<crystal_structure.num_atoms(); iatom++) // atom sum
                    {
                        int iatomic_species = crystal_structure.index_of_atoms()[iatom];
                        int num_projectors = potentials.pseudo_lbeta()[iatomic_species].size();

                        phase_atom.resize(num_G_at_k);
                        for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++) 
                        {
                            kGvect = crystal_structure.reciprocal_vectors().transpose()
                                *(kvector + plane_wave_basis.get_Gvector(Gindex_at_k(ipw_at_k)).cast<double>());
                            double phase = kGvect.transpose() * crystal_structure.atomic_position_cartesian()[iatom];

                            phase_atom(ipw_at_k) = std::cos(phase) + I*std::sin(phase);
                        } // set phase_atom

                        for (int jproj=0; jproj<num_projectors; jproj++)
                        {
                            int lbeta = potentials.pseudo_lbeta()[iatomic_species][jproj];
                            for (int im=0; im<2*lbeta+1; im++)
                            {
                                Complex coeff = 
                                    (phase_atom.array() * projector[iatomic_species][jproj][im].array() * 
                                     phi[ispin][ik][jband][jspinor].array()).sum();
                                for (int iproj=0; iproj<num_projectors; iproj++)
                                {
                                    if (std::abs(potentials.pseudo_dij()[iatomic_species](iproj, jproj))>1e-8)
                                    {
                                        if (potentials.pseudo_lbeta()[iatomic_species][iproj] != lbeta) { error_messages::stop("dij is not diagonal w.r.t. l in pseudopot."); }

                                        H1phi[ispin][ik][jband][jspinor] = H1phi[ispin][ik][jband][jspinor].array() + 
                                            (potentials.pseudo_dij()[iatomic_species](iproj, jproj) * coeff) *
                                            (phase_atom.array() * projector[iatomic_species][iproj][im].array()).conjugate();
                                    } // dij > 1e-8
                                } // iproj
                            } // im
                        } // jproj
                    } // iatom

                    // (2/2) pseudo potential (local part)
                    plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, // isym = 0
                                                         false, // time-rersal not used for isym=0
                                                         phi[ispin][ik][jband][jspinor], orbital,
                                                         method.calc_mode());
                    plane_wave_basis.FFT_backward(orbital, orbital); // G-space to R-space
                    V_orbital = potentials.Vpp_local().array() * orbital.array();
                    plane_wave_basis.FFT_forward(V_orbital, V_orbital); // R-space to G-space
                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                    {
                        H1phi[ispin][ik][jband][jspinor](ipw_at_k) += V_orbital(Gindex_at_k(ipw_at_k));
                    }
                } // jspinor
            } // jband
        } // ik
    } // ispin
}
