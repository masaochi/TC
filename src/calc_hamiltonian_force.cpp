// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

// calculate the Hellmann-Feynman force
std::vector<Eigen::Vector3cd> calc_hamiltonian::force(const Parallelization &parallelization, 
                                                      const Method &method,
                                                      const CrystalStructure &crystal_structure,
                                                      const Potentials &potentials,
                                                      const Spin &spin, const Kpoints &kpoints,
                                                      PlaneWaveBasis &plane_wave_basis,
                                                      const BlochStates &bloch_states,
                                                      const TotalEnergy &total_energy,
                                                      std::ostream *ost)
{
    const bool am_i_mpi_rank0 = parallelization.am_i_mpi_rank0();

    const bool is_bitc = method.calc_method()=="BITC" ? true : false;
    const int num_atoms = crystal_structure.num_atoms();
    const int num_independent_spins = spin.num_independent_spins();
    const int num_irreducible_kpoints_scf = kpoints.num_irreducible_kpoints_scf();
    const int num_spinor = (spin.is_spinor()==false) ? 1 : 2;
    const int size_FFT_grid = plane_wave_basis.size_FFT_grid();

    // initialize return_vector    
    std::vector<Eigen::Vector3cd> return_vector(num_atoms);
    for (int iatom=0; iatom<num_atoms; iatom++) { return_vector[iatom] = { 0.0, 0.0, 0.0 }; }

    // temporary variables
    Eigen::Vector3d kvect, kGvect;
    Eigen::VectorXcd phi, chi, phase_atom;
    std::vector<Eigen::VectorXcd> phase_atom_3d(3);
    std::vector<Eigen::VectorXcd> projector;
    std::vector<std::vector<std::vector<Eigen::VectorXcd> > > projector_fftgrid(crystal_structure.num_atomic_species());
    std::vector<std::vector<Complex> > coeff; // coeff[iproj][im]
    std::vector<std::vector<Eigen::Vector3cd> > coeff_3d; // coeff_3d[iproj][im](0:2)
    std::vector<std::vector<Complex> > coeff_left; // used only for BITC
    std::vector<std::vector<Eigen::Vector3cd> > coeff_3d_left;

    std::vector<std::vector<Complex> > &coeff_ref = is_bitc ? coeff_left : coeff;
    std::vector<std::vector<Eigen::Vector3cd> > &coeff_3d_ref = is_bitc ? coeff_3d_left : coeff_3d;

    // start calculation
    if (am_i_mpi_rank0)
    {
        for (int iatom=0; iatom<num_atoms; iatom++)
        {
            // (1) Ewald
            return_vector[iatom] = return_vector[iatom] 
                + total_energy.differentiate_ewald_energy(crystal_structure)[iatom];

            // (2) local pseudopot.
            std::vector<std::vector<Eigen::VectorXcd> > derivative_pseudo_local; // [iatom][3](ipw)
            potentials.differentiate_local_potential(crystal_structure, plane_wave_basis, derivative_pseudo_local);
            for (int idim=0; idim<3; idim++)
            {
                return_vector[iatom](idim) +=
                    ((bloch_states.density()[0] + bloch_states.density()[1]).array()
                     * derivative_pseudo_local[iatom][idim].array()).sum()
                    / static_cast<double>(size_FFT_grid * kpoints.num_kpoints()); // density should be divided with Nk
            }
        } // iatom
    } // am_i_mpi_rank0

    // (3) non-local pseudopot.
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

    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
        {
            for (int isym=0; isym<kpoints.kvectors_scf()[ik].size(); isym++)
            {
                kvect = kpoints.kvectors_scf()[ik][isym];
                const int num_G_at_k = plane_wave_basis.num_G_at_k_scf()[ik];
                const Eigen::VectorXi Gindex_at_k = plane_wave_basis.Gindex_at_k_scf()[ispin][ik][isym];

                bool projector_calculated = false;
                for (int iband=0; iband<bloch_states.num_occupied_bands()[ispin][ik]; iband++)
                {
                    if (!parallelization.is_assigned_all_kpoints_occupied_bands()[ispin][ik][isym][iband]) { continue; }
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
                                * (kvect + plane_wave_basis.get_Gvector(Gindex_at_k(ipw_at_k)).cast<double>());
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
                            } // if (qG<1e-8)
                        } // ipw_at_k
                        
                        // calculate <beta|k+G> (except the exp[i(k+G)R_atom] factor)
                        for (int iatomic_species=0; iatomic_species<crystal_structure.num_atomic_species(); iatomic_species++)
                        {
                            int num_projectors = potentials.pseudo_lbeta()[iatomic_species].size();

                            projector_fftgrid[iatomic_species].resize(num_projectors);
                            for (int iproj=0; iproj<num_projectors; iproj++)
                            {
                                projector = potentials.projector_nonlocal(crystal_structure, plane_wave_basis,
                                                                          Gindex_at_k,
                                                                          kvect, ylm,                    
                                                                          iatomic_species, iproj);
                                
                                int lbeta = potentials.pseudo_lbeta()[iatomic_species][iproj];

                                projector_fftgrid[iatomic_species][iproj].resize(2*lbeta+1);
                                for (int im=0; im<2*lbeta+1; im++)
                                {
                                    projector_fftgrid[iatomic_species][iproj][im] = Eigen::VectorXcd::Zero(size_FFT_grid);
                                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                                    {
                                        projector_fftgrid[iatomic_species][iproj][im](Gindex_at_k(ipw_at_k))
                                            = projector[im](ipw_at_k);
                                    }
                                } // im
                            } // iproj
                        } // iatomic_species
                        projector_calculated = true;
                    } // if (!projector_calculated)
 
                    // calculate orbitals
                    Eigen::VectorXcd &chi_ref = is_bitc ? chi : phi;
                    plane_wave_basis.get_orbital_FFTgrid(ispin, ik, isym,
                                                         kpoints.is_time_reversal_used_at_k()[ik][isym],
                                                         bloch_states.phik_scf()[ispin][ik][iband][0], phi,
                                                         "SCF");
                    if (is_bitc)
                    {
                        plane_wave_basis.get_orbital_FFTgrid(ispin, ik, isym,
                                                             kpoints.is_time_reversal_used_at_k()[ik][isym],
                                                             bloch_states.phik_left_scf()[ispin][ik][iband][0], chi,
                                                             "SCF");
                    }
                    
                    for (int jspinor=0; jspinor<num_spinor; jspinor++)
                    {
                        for (int iatom=0; iatom<num_atoms; iatom++)
                        {
                            int iatomic_species = crystal_structure.index_of_atoms()[iatom];
                            int num_projectors = potentials.pseudo_lbeta()[iatomic_species].size();
                        
                            // set phase_atom on the FFT grid
                            phase_atom = Eigen::VectorXcd::Zero(size_FFT_grid);
                            for (int idim=0; idim<3; idim++) { phase_atom_3d[idim] = Eigen::VectorXcd::Zero(size_FFT_grid); }
                            for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++) 
                            {
                                kGvect = crystal_structure.reciprocal_vectors().transpose()
                                    * (kvect + plane_wave_basis.get_Gvector(Gindex_at_k(ipw_at_k)).cast<double>());
                                double phase = kGvect.transpose() * crystal_structure.atomic_position_cartesian()[iatom];
                                
                                phase_atom(Gindex_at_k(ipw_at_k)) = std::cos(phase) + I*std::sin(phase);
                                for (int idim=0; idim<3; idim++)
                                {
                                    phase_atom_3d[idim](Gindex_at_k(ipw_at_k)) = 
                                        I * kGvect(idim) * phase_atom(Gindex_at_k(ipw_at_k));
                                }
                            }
                            
                            // initialize coefficients
                            coeff.resize(num_projectors); // coeff[iproj][im]
                            coeff_3d.resize(num_projectors); // coeff_3d[iproj][im](0:2)
                            if (is_bitc)
                            {
                                coeff_left.resize(num_projectors);
                                coeff_3d_left.resize(num_projectors);
                            }
                            for (int iproj=0; iproj<num_projectors; iproj++)
                            {
                                int lbeta = potentials.pseudo_lbeta()[iatomic_species][iproj];
                                
                                coeff[iproj].resize(2*lbeta+1);
                                coeff_3d[iproj].resize(2*lbeta+1);
                                if (is_bitc)
                                {
                                    coeff_left[iproj].resize(2*lbeta+1);
                                    coeff_3d_left[iproj].resize(2*lbeta+1);
                                }
                            }
                            
                            // set coeff, coeff_3d
                            for (int iproj=0; iproj<num_projectors; iproj++)
                            {
                                int lbeta = potentials.pseudo_lbeta()[iatomic_species][iproj];
                                for (int im=0; im<2*lbeta+1; im++)
                                {
                                    coeff[iproj][im] = 
                                        (phase_atom.array() * projector_fftgrid[iatomic_species][iproj][im].array() * 
                                         phi.array()).sum();
                                
                                    for (int idim=0; idim<3; idim++)
                                    {
                                        coeff_3d[iproj][im](idim) =
                                            (phase_atom_3d[idim].array() * projector_fftgrid[iatomic_species][iproj][im].array() * 
                                             phi.array()).sum();
                                    }
                                    if (is_bitc)
                                    {
                                        coeff_left[iproj][im] = 
                                            (phase_atom.array() * projector_fftgrid[iatomic_species][iproj][im].array() * 
                                             chi_ref.array()).sum();

                                        for (int idim=0; idim<3; idim++)
                                        {                                        
                                            coeff_3d_left[iproj][im](idim) =
                                                (phase_atom_3d[idim].array() * projector_fftgrid[iatomic_species][iproj][im].array() * 
                                                 chi_ref.array()).sum();
                                        }
                                    } // is_bitc
                                } // im
                            } // iproj                            
                            
                            // calculate return_vector
                            for (int jproj=0; jproj<num_projectors; jproj++)
                            {
                                int lbeta = potentials.pseudo_lbeta()[iatomic_species][jproj];
                                for (int im=0; im<2*lbeta+1; im++)
                                {
                                    for (int iproj=0; iproj<num_projectors; iproj++)
                                    {
                                        if (std::abs(potentials.pseudo_dij()[iatomic_species](iproj, jproj))>1e-8)
                                        {
                                            if (potentials.pseudo_lbeta()[iatomic_species][iproj] != lbeta) { error_messages::stop("dij is not diagonal w.r.t. l in pseudopot."); }
                                            
                                            for (int idim=0; idim<3; idim++)
                                            {
                                                return_vector[iatom](idim) +=
                                                    bloch_states.filling()[ispin][ik][iband] / kpoints.num_kpoints() *
                                                    potentials.pseudo_dij()[iatomic_species](iproj, jproj) *
                                                    (std::conj(coeff_ref[iproj][im]) * coeff_3d[jproj][im](idim) +
                                                     std::conj(coeff_3d_ref[iproj][im](idim)) * coeff[jproj][im]);
                                            } // idim
                                        } // dij > 1e-8
                                    } // iproj
                                } // im
                            } // jproj
                        } // iatom
                    } // jspinor
                } // iband
            } // isym
        } // ik
    } // ispin

    // reduce (not allreduce!)
    auto return_vector_local = return_vector;
    for (int iatom=0; iatom<num_atoms; iatom++)
    {
        MPI_Reduce(return_vector_local[iatom].data(),
                   return_vector[iatom].data(),
                   return_vector[iatom].size(),
                   MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    return return_vector;
}
