// [class Diagonalization]
// scf_loop, convergence check, control paramters, etc.

#include "include/header.hpp"

// non-collinear spin not supported
void Diagonalization::block_davidson(Parallelization &parallelization, 
                                     MyClock &my_clock, const Method &method,
                                     const CrystalStructure &crystal_structure,
                                     const Symmetry &symmetry, Potentials &potentials,
                                     const Spin &spin, const Kpoints &kpoints,
                                     PlaneWaveBasis &plane_wave_basis, 
                                     BlochStates &bloch_states, TotalEnergy &total_energy,
                                     std::ostream *ost) const
{
    const bool am_i_mpi_rank0 = parallelization.am_i_mpi_rank0();
    const bool is_bitc = method.calc_method()=="BITC" ? true : false;
    const bool uses_3body = (method.calc_method()=="TC" || method.calc_method()=="BITC") ? true : false;

    assert(is_bitc || !biortho_david_);

    const int num_independent_spins = spin.num_independent_spins();
    const int num_irreducible_kpoints = method.calc_mode()=="SCF" ? 
        kpoints.num_irreducible_kpoints_scf() : kpoints.num_irreducible_kpoints_band();
    const int num_spinor = (spin.is_spinor()==false) ? 1 : 2;
    if (num_spinor==2) { error_messages::stop("Non-collinear calculation is not supported..."); }
    const std::vector<int> num_bands_tc = method.calc_mode()=="SCF" ?
        bloch_states.num_bands_scf() : bloch_states.num_bands_band();

    const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_ref
        = method.calc_mode()=="SCF" ? bloch_states.phik_scf() : bloch_states.phik_band();
    const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_left_ref
        = method.calc_mode()=="SCF" ? bloch_states.phik_left_scf() : bloch_states.phik_left_band();
    const std::vector<std::vector<std::vector<Complex> > > &eigenvalues_ref
        = method.calc_mode()=="SCF" ? bloch_states.eigenvalues_scf() : bloch_states.eigenvalues_band();

    // subspaace Hamiltonian: H_subspace[spin][irreducible kpoints](band, band)
    std::vector<std::vector<Eigen::MatrixXcd> > H_subspace, H1_subspace, H2_subspace, H3_subspace; // 1,2,3-body terms and total
    if (am_i_mpi_rank0)
    {
        H_subspace.resize(num_independent_spins);
        H1_subspace.resize(num_independent_spins);
        H2_subspace.resize(num_independent_spins);
        if (uses_3body) { H3_subspace.resize(num_independent_spins); }
        for (int ispin=0; ispin<num_independent_spins; ispin++)
        {
            H_subspace[ispin].resize(num_irreducible_kpoints);
            H1_subspace[ispin].resize(num_irreducible_kpoints);
            H2_subspace[ispin].resize(num_irreducible_kpoints);
            if (uses_3body) { H3_subspace[ispin].resize(num_irreducible_kpoints); }
        }
    } // am_i_mpi_rank0

    // V[spin][irreducible kpoints][subspace(band)](num_G_at_k)
    std::vector<std::vector<std::vector<Eigen::VectorXcd> > > V(num_independent_spins); // trial vectors in Krylov subspace
    std::vector<std::vector<std::vector<Eigen::VectorXcd> > > H1V, H2V, H3V; // Hamiltonian (1,2,3-body terms)*V
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        const int nbands = num_bands_tc[ispin];

        V[ispin].resize(num_irreducible_kpoints);
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            V[ispin][ik].resize(nbands); // dimension of V will be later increased
            for (int iband=0; iband<nbands; iband++) // initialize V by phik
            {
                V[ispin][ik][iband] = phik_ref[ispin][ik][iband][0];
                if (is_bitc && !biortho_david_) // orthonormalization required
                {
                    if (!Gram_Schmidt(V[ispin][ik], iband, V[ispin][ik][iband]))
                    {
                        error_messages::stop("Linear independence breaks in block-Davidson algorithm.");
                    }
                } // is_bitc
            } // iband
        } // ik
    } // ispin
    if (am_i_mpi_rank0)
    {
        H1V.resize(num_independent_spins);
        H2V.resize(num_independent_spins);
        if (uses_3body) { H3V.resize(num_independent_spins); }
        for (int ispin=0; ispin<num_independent_spins; ispin++)
        {
            H1V[ispin].resize(num_irreducible_kpoints);
            H2V[ispin].resize(num_irreducible_kpoints);
            if (uses_3body) { H3V[ispin].resize(num_irreducible_kpoints);}
        }
    } // am_i_mpi_rank0

    // left orbitals
    std::vector<std::vector<std::vector<Eigen::VectorXcd> > > V_left; // trial vectors in Krylov subspace
    if (biortho_david_)
    {
        V_left.resize(num_independent_spins);
        for (int ispin=0; ispin<num_independent_spins; ispin++)
        {
            const int nbands = num_bands_tc[ispin];
            
            V_left[ispin].resize(num_irreducible_kpoints);
            for (int ik=0; ik<num_irreducible_kpoints; ik++)
            {
                V_left[ispin][ik].resize(nbands); // dimension of V will be later increased
                for (int iband=0; iband<nbands; iband++) // initialize V by phik
                {
                    V_left[ispin][ik][iband] = phik_left_ref[ispin][ik][iband][0];
                } // iband
            } // ik
        } // ispin
    } // biortho_david_

    // phi[spin][irreducible kpoint][band][ispinor][num_G_at_k] Note: G-space vector
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > phi(num_independent_spins);
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > H1phi(num_independent_spins); // one-body terms
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > H2phi(num_independent_spins); // two-body terms
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > H3phi(num_independent_spins); // three-body terms
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        phi[ispin].resize(num_irreducible_kpoints);
        H1phi[ispin].resize(num_irreducible_kpoints);
        H2phi[ispin].resize(num_irreducible_kpoints);
        H3phi[ispin].resize(num_irreducible_kpoints);
    }

    // left orbitals
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > phi_left;
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > H1phi_left;
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > H2phi_left;
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > H3phi_left;
    if (biortho_david_)
    {
        phi_left.resize(num_independent_spins);
        H1phi_left.resize(num_independent_spins);
        H2phi_left.resize(num_independent_spins);
        H3phi_left.resize(num_independent_spins);

        for (int ispin=0; ispin<num_independent_spins; ispin++)
        {
            phi_left[ispin].resize(num_irreducible_kpoints);
            H1phi_left[ispin].resize(num_irreducible_kpoints);
            H2phi_left[ispin].resize(num_irreducible_kpoints);
            H3phi_left[ispin].resize(num_irreducible_kpoints);
        }
    }

    for (int iter_inner=0; iter_inner<num_refresh_david_; iter_inner++) // update trial vectors while keeping the Fock operator unchanged (refresh)
    {
        // resize H*V, phi, H*phi
        for (int ispin=0; ispin<num_independent_spins; ispin++)
        {
            const int nbands = num_bands_tc[ispin];
            for (int ik=0; ik<num_irreducible_kpoints; ik++)
            {
                const int num_G_at_k = method.calc_mode()=="SCF" ?
                    plane_wave_basis.num_G_at_k_scf()[ik] : plane_wave_basis.num_G_at_k_band()[ik];

                // resize phi, H*phi
                phi[ispin][ik].resize(nbands);
                H1phi[ispin][ik].resize(nbands);
                H2phi[ispin][ik].resize(nbands);
                H3phi[ispin][ik].resize(nbands);
                for (int iband=0; iband<nbands; iband++)
                {
                    phi[ispin][ik][iband].resize(num_spinor);
                    H1phi[ispin][ik][iband].resize(num_spinor);
                    H2phi[ispin][ik][iband].resize(num_spinor);
                    H3phi[ispin][ik][iband].resize(num_spinor);
                    for (int ispinor=0; ispinor<num_spinor; ispinor++)
                    {
                        H1phi[ispin][ik][iband][ispinor].resize(num_G_at_k);
                        H2phi[ispin][ik][iband][ispinor].resize(num_G_at_k);
                        if (uses_3body) { H3phi[ispin][ik][iband][ispinor].resize(num_G_at_k); }
                    }
                }

                if (biortho_david_)
                {
                    phi_left[ispin][ik].resize(nbands);
                    H1phi_left[ispin][ik].resize(nbands);
                    H2phi_left[ispin][ik].resize(nbands);
                    H3phi_left[ispin][ik].resize(nbands);
                    for (int iband=0; iband<nbands; iband++)
                    {
                        phi_left[ispin][ik][iband].resize(num_spinor);
                        H1phi_left[ispin][ik][iband].resize(num_spinor);
                        H2phi_left[ispin][ik][iband].resize(num_spinor);
                        H3phi_left[ispin][ik][iband].resize(num_spinor);
                        for (int ispinor=0; ispinor<num_spinor; ispinor++)
                        {
                            H1phi_left[ispin][ik][iband][ispinor].resize(num_G_at_k);
                            H2phi_left[ispin][ik][iband][ispinor].resize(num_G_at_k);
                            if (uses_3body) { H3phi_left[ispin][ik][iband][ispinor].resize(num_G_at_k); }
                        }
                    }
                } // biortho_david_

                if (am_i_mpi_rank0)
                {
                    // resize H*V
                    H1V[ispin][ik].resize(nbands); // dimension of H*V will be later increased
                    H2V[ispin][ik].resize(nbands);
                    if (uses_3body) { H3V[ispin][ik].resize(nbands); } 
                } // am_i_mpi_rank0
            } // ik
        } // ispin

        for (int iblock=0; iblock<(max_num_blocks_david_-1); iblock++) // increase subspace dimension
        {
            // set phi by V
            for (int ispin=0; ispin<num_independent_spins; ispin++)
            {
                for (int ik=0; ik<num_irreducible_kpoints; ik++)
                {
                    const int nbands = phi[ispin][ik].size(); // handled in this iblock-loop
                    const int nbands_fixed = V[ispin][ik].size() - nbands; // handled in previous iblock-loops
                    for (int iband=0; iband<nbands; iband++)
                    {
                        phi[ispin][ik][iband][0] = V[ispin][ik][iband + nbands_fixed];
                        if (biortho_david_) { phi_left[ispin][ik][iband][0] = V_left[ispin][ik][iband + nbands_fixed]; }
                    }
                }
            }

            // set parallelization
            parallelization.set_orbital_parallelization(bloch_states, kpoints, method.calc_mode(), phi);

            // Hphi = Hamiltonian * phi
            calc_hamiltonian::all(parallelization, my_clock, method,
                                  crystal_structure, symmetry,
                                  potentials, spin, kpoints,
                                  plane_wave_basis, bloch_states,
                                  phi, H1phi, H2phi, H3phi, ost);

            if (biortho_david_)
            {
                // Hphi_left = H^{dag} * phi_left
                
                // Take Hermitian-conjugate of Hamiltonian:       
                potentials.jastrow.take_Hermitian_conj();
                bloch_states.switch_left_right_orbitals(method.calc_mode());

                calc_hamiltonian::all(parallelization, my_clock, method,
                                      crystal_structure, symmetry,
                                      potentials, spin, kpoints,
                                      plane_wave_basis, bloch_states,
                                      phi_left, H1phi_left, H2phi_left, H3phi_left, ost);

                // Undo Hermitian-conjugate:
                potentials.jastrow.take_Hermitian_conj();
                bloch_states.switch_left_right_orbitals(method.calc_mode());
            } // biortho_david_

            // For setting new trial vectors
            std::vector<std::vector<int> > nbands_old(num_independent_spins);
            std::vector<std::vector<int> > nbands_new(num_independent_spins);
            for (int ispin=0; ispin<num_independent_spins; ispin++)
            {
                nbands_old[ispin].resize(num_irreducible_kpoints);
                nbands_new[ispin].resize(num_irreducible_kpoints);
                for (int ik=0; ik<num_irreducible_kpoints; ik++)
                {
                    nbands_old[ispin][ik] = V[ispin][ik].size();
                }
            }

            if (am_i_mpi_rank0)
            {
                for (int ispin=0; ispin<num_independent_spins; ispin++)
                {
                    for (int ik=0; ik<num_irreducible_kpoints; ik++)
                    {
                        const int nbands = phi[ispin][ik].size();
                        const int nbands_fixed = V[ispin][ik].size() - nbands; // handled in previous iblock-loops

                        const Eigen::Vector3d kvector = method.calc_mode()=="SCF" ? 
                            kpoints.kvectors_scf()[ik][0] : kpoints.kvectors_band()[ik][0];
                        const Eigen::VectorXi Gindex_at_k = method.calc_mode()=="SCF" ?
                            plane_wave_basis.Gindex_at_k_scf()[ispin][ik][0] :  plane_wave_basis.Gindex_at_k_band()[ispin][ik][0];

                        Eigen::Vector3d kGvect; // k+G
                        Eigen::VectorXd kGvect2(Gindex_at_k.size()); // |k+G|^2
                        for (int ipw_at_k=0; ipw_at_k<Gindex_at_k.size(); ipw_at_k++)
                        {
                            kGvect = crystal_structure.reciprocal_vectors().transpose()
                                *(kvector + plane_wave_basis.get_Gvector(Gindex_at_k(ipw_at_k)).cast<double>());
                            kGvect2(ipw_at_k) = kGvect.squaredNorm();
                        }

                        for (int iband=0; iband<nbands; iband++)
                        {
                            int ndim_total_temporary = V[ispin][ik].size();
                        
                            H1V[ispin][ik][iband + nbands_fixed] = H1phi[ispin][ik][iband][0];
                            H2V[ispin][ik][iband + nbands_fixed] = H2phi[ispin][ik][iband][0];
                            if (uses_3body) { H3V[ispin][ik][iband + nbands_fixed] = H3phi[ispin][ik][iband][0]; }

                            // make a new trial vector (V) by preconditioned phi
                            make_a_new_trial_vector(kGvect2,
                                                    phi[ispin][ik][iband][0], H1phi[ispin][ik][iband][0],
                                                    H2phi[ispin][ik][iband][0], H3phi[ispin][ik][iband][0],
                                                    eigenvalues_ref[ispin][ik][iband], uses_3body);
                            if (biortho_david_)
                            {
                                make_a_new_trial_vector(kGvect2,
                                                        phi_left[ispin][ik][iband][0], H1phi_left[ispin][ik][iband][0],
                                                        H2phi_left[ispin][ik][iband][0], H3phi_left[ispin][ik][iband][0],
                                                        std::conj(eigenvalues_ref[ispin][ik][iband]), uses_3body);
                            }

                            bool is_new_vector_linearly_independent = false;
                            if (!biortho_david_)
                            {
                                is_new_vector_linearly_independent 
                                    = Gram_Schmidt(V[ispin][ik], ndim_total_temporary, phi[ispin][ik][iband][0]);
                            }
                            else
                            {
                                is_new_vector_linearly_independent 
                                    = Gram_Schmidt_biortho(V[ispin][ik], V_left[ispin][ik], ndim_total_temporary, 
                                                           phi[ispin][ik][iband][0], phi_left[ispin][ik][iband][0]);
                            }

                            if (is_new_vector_linearly_independent)
                            {
                                V[ispin][ik].push_back(phi[ispin][ik][iband][0]);
                                if (biortho_david_) { V_left[ispin][ik].push_back(phi_left[ispin][ik][iband][0]); }
                            }
                        } // iband
                        nbands_new[ispin][ik] = V[ispin][ik].size();
                    } // ik
                } // ispin
                my_clock.print_time_from_save(ost, "preparation for the next loop (preconditioning & Gram-Schmidt)");
            } // am_i_mpi_rank0

            // Bcast new trial vectors
            for (int ispin=0; ispin<num_independent_spins; ispin++)
            {
                for (int ik=0; ik<num_irreducible_kpoints; ik++)
                {
                    const int num_G_at_k = method.calc_mode()=="SCF" ?
                        plane_wave_basis.num_G_at_k_scf()[ik] : plane_wave_basis.num_G_at_k_band()[ik];

                    // Bcast nbands_new
                    MPI_Bcast(&nbands_new[ispin][ik], 1, MPI_INT, 0, MPI_COMM_WORLD);

                    // Bcast V
                    if (!am_i_mpi_rank0) { V[ispin][ik].resize(nbands_new[ispin][ik]); }
                    for (int iband=nbands_old[ispin][ik]; iband<nbands_new[ispin][ik]; iband++)
                    {
                        if (!am_i_mpi_rank0) { V[ispin][ik][iband].resize(num_G_at_k); }
                        MPI_Bcast(V[ispin][ik][iband].data(), V[ispin][ik][iband].size(),
                                  MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
                    }
                    if (biortho_david_)
                    {
                        if (!am_i_mpi_rank0) { V_left[ispin][ik].resize(nbands_new[ispin][ik]); }
                        for (int iband=nbands_old[ispin][ik]; iband<nbands_new[ispin][ik]; iband++)
                        {
                            if (!am_i_mpi_rank0) { V_left[ispin][ik][iband].resize(num_G_at_k); }
                            MPI_Bcast(V_left[ispin][ik][iband].data(), V_left[ispin][ik][iband].size(),
                                      MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
                        }
                    } // biortho_david_

                    // resize phi, H*phi
                    phi[ispin][ik].resize(nbands_new[ispin][ik] - nbands_old[ispin][ik]);
                    H1phi[ispin][ik].resize(nbands_new[ispin][ik] - nbands_old[ispin][ik]);
                    H2phi[ispin][ik].resize(nbands_new[ispin][ik] - nbands_old[ispin][ik]);
                    H3phi[ispin][ik].resize(nbands_new[ispin][ik] - nbands_old[ispin][ik]); 
                    for (int iband=0; iband<nbands_new[ispin][ik] - nbands_old[ispin][ik]; iband++)
                    {
                        phi[ispin][ik][iband].resize(num_spinor);
                        H1phi[ispin][ik][iband].resize(num_spinor); 
                        H2phi[ispin][ik][iband].resize(num_spinor);
                        H3phi[ispin][ik][iband].resize(num_spinor);
                        for (int ispinor=0; ispinor<num_spinor; ispinor++)
                        {
                            H1phi[ispin][ik][iband][ispinor].resize(num_G_at_k);
                            H2phi[ispin][ik][iband][ispinor].resize(num_G_at_k);
                            if (uses_3body) { H3phi[ispin][ik][iband][ispinor].resize(num_G_at_k); }
                        }
                    } // iband
                    if (biortho_david_)
                    {
                        phi_left[ispin][ik].resize(nbands_new[ispin][ik] - nbands_old[ispin][ik]);
                        H1phi_left[ispin][ik].resize(nbands_new[ispin][ik] - nbands_old[ispin][ik]);
                        H2phi_left[ispin][ik].resize(nbands_new[ispin][ik] - nbands_old[ispin][ik]);
                        H3phi_left[ispin][ik].resize(nbands_new[ispin][ik] - nbands_old[ispin][ik]); 
                        for (int iband=0; iband<nbands_new[ispin][ik] - nbands_old[ispin][ik]; iband++)
                        {
                            phi_left[ispin][ik][iband].resize(num_spinor);
                            H1phi_left[ispin][ik][iband].resize(num_spinor); 
                            H2phi_left[ispin][ik][iband].resize(num_spinor);
                            H3phi_left[ispin][ik][iband].resize(num_spinor);
                            for (int ispinor=0; ispinor<num_spinor; ispinor++)
                            {
                                H1phi_left[ispin][ik][iband][ispinor].resize(num_G_at_k);
                                H2phi_left[ispin][ik][iband][ispinor].resize(num_G_at_k);
                                if (uses_3body) { H3phi_left[ispin][ik][iband][ispinor].resize(num_G_at_k); }
                            }
                        } // iband
                    } // biortho_david_

                    // resize H*V
                    if (am_i_mpi_rank0)
                    {
                        H1V[ispin][ik].resize(nbands_new[ispin][ik]);
                        H2V[ispin][ik].resize(nbands_new[ispin][ik]);
                        if (uses_3body) { H3V[ispin][ik].resize(nbands_new[ispin][ik]); }
                    } // am_i_mpi_rank0
                } // ik
            } // ispin
            if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "preparation for the next loop (Bcast)"); }
        } // iblock

        // To construct a subspace Hamiltonian, we need to evaluate remaining matrix elements of HV
        for (int ispin=0; ispin<num_independent_spins; ispin++)
        {
            for (int ik=0; ik<num_irreducible_kpoints; ik++)
            {
                const int nbands = phi[ispin][ik].size();
                const int nbands_fixed = V[ispin][ik].size() - nbands;
                for (int iband=0; iband<nbands; iband++)
                {
                    phi[ispin][ik][iband][0] = V[ispin][ik][iband + nbands_fixed];
                }
            }
        }

        // set parallelization
        parallelization.set_orbital_parallelization(bloch_states, kpoints, method.calc_mode(), phi);

        calc_hamiltonian::all(parallelization, my_clock, method,
                              crystal_structure, symmetry,
                              potentials, spin, kpoints,
                              plane_wave_basis, bloch_states,
                              phi, H1phi, H2phi, H3phi, ost);
        if (am_i_mpi_rank0)
        {
            *ost << "  diagonalize subspace Hamiltonian" << std::endl;
            for (int ispin=0; ispin<num_independent_spins; ispin++)
            {
                for (int ik=0; ik<num_irreducible_kpoints; ik++)
                {
                    const int nbands = phi[ispin][ik].size();
                    const int nbands_fixed = V[ispin][ik].size() - nbands;
                    for (int iband=0; iband<nbands; iband++)
                    {
                        H1V[ispin][ik][iband + nbands_fixed] = H1phi[ispin][ik][iband][0];
                        H2V[ispin][ik][iband + nbands_fixed] = H2phi[ispin][ik][iband][0];
                        if (uses_3body) { H3V[ispin][ik][iband + nbands_fixed] = H3phi[ispin][ik][iband][0]; }
                    }

                    const int ndim_total = V[ispin][ik].size();

                    // resize H_subspace
                    H_subspace[ispin][ik].resize(ndim_total, ndim_total);
                    H1_subspace[ispin][ik].resize(ndim_total, ndim_total);
                    H2_subspace[ispin][ik].resize(ndim_total, ndim_total);
                    if (uses_3body) { H3_subspace[ispin][ik].resize(ndim_total, ndim_total); }

                    for (int idim=0; idim<ndim_total; idim++)
                    {
                        Eigen::VectorXcd &V_ref = biortho_david_ ? 
                            V_left[ispin][ik][idim] : V[ispin][ik][idim];
                        for (int jdim=0; jdim<ndim_total; jdim++)
                        {
                            H1_subspace[ispin][ik](idim, jdim) = V_ref.dot(H1V[ispin][ik][jdim]);
                            H2_subspace[ispin][ik](idim, jdim) = V_ref.dot(H2V[ispin][ik][jdim]);
                            if (uses_3body)
                            {
                                H3_subspace[ispin][ik](idim, jdim) = V_ref.dot(H3V[ispin][ik][jdim]);
                            }
                        }
                    }
                    H_subspace[ispin][ik] = H1_subspace[ispin][ik] + H2_subspace[ispin][ik];
                    if (uses_3body) { H_subspace[ispin][ik] += H3_subspace[ispin][ik]; }

                    // diagonalize H_subspace
                    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> diag;
                    diag.compute(H_subspace[ispin][ik]); // diagonalize
                    if (diag.info()!=Eigen::Success) { error_messages::stop("diagonalization of H_subspace failed"); }

                    // A = X*D*X^{-1} where X = diag.eigenvectors()
                    Eigen::VectorXcd eigenvalues = diag.eigenvalues(); // e.g. eigenvalues(0) = 0th eigenvalue
                    Eigen::MatrixXcd eigenvectors = diag.eigenvectors(); // e.g. eigenvectors(*,0) = 0th (normalized) eigenvector
                    if (!biortho_david_)
                    {
                        sort_eigen(eigenvalues, eigenvectors); // Re(e0) <= Re(e1) <= ...
                    }
                    else
                    {
                        sort_eigen_overlap(eigenvalues, eigenvectors); // large overlap of wfc between the present and preivous loops
                    }

                    // Note: X^{-1}*A = D*X^{-1}. eigenvectors_left(*,0) = 0th left-eigenvector 
                    Eigen::MatrixXcd eigenvectors_left;
                    if (is_bitc) { eigenvectors_left = eigenvectors.inverse().adjoint(); }

                    // reset trial vectors V_tmp = V*X: right eigenvectors of H.
                    //    V and V_left^{-1}^{dag} are the basis functions of right- and left-orbitals of H.
                    //    i.e., A = V_left^{dag}*H*V. X and X^{-1}^{dag} are coefficients of eigenvectors of H in terms of such basis functions.
                    //    Then, <V_left*X^{-1}^{dag} | H | V*X> = X^{-1} * A * X = D,
                    //    where <V_left*X^{-1}^{dag} | V*X> = X^{-1} * (V_left^{dag} * V) * X = I
                    //    Therefore, V*X and V_left*X^{-1}^{dag} are good estimates of the right- and left-eigenvectors of H.
                    // [Note: V_left = V for uses_biortho_david==false]
                    std::vector<Eigen::VectorXcd> V_tmp(num_bands_tc[ispin]);
                    for (int iband=0; iband<num_bands_tc[ispin]; iband++)
                    {
                        V_tmp[iband] = Eigen::VectorXcd::Zero(V[ispin][ik][0].size());
                        for (int ibasis=0; ibasis<ndim_total; ibasis++)
                        {
                            for (int ipw_at_k=0; ipw_at_k<V[ispin][ik][0].size(); ipw_at_k++)
                            {
                                V_tmp[iband](ipw_at_k) += eigenvectors(ibasis, iband) * V[ispin][ik][ibasis](ipw_at_k);
                            }
                        }
                    }
                    std::vector<Eigen::VectorXcd> V_tmp_left;
                    if (biortho_david_)
                    {
                        V_tmp_left.resize(num_bands_tc[ispin]);
                        for (int iband=0; iband<num_bands_tc[ispin]; iband++)
                        {
                            V_tmp_left[iband] = Eigen::VectorXcd::Zero(V_left[ispin][ik][0].size());
                            for (int ibasis=0; ibasis<ndim_total; ibasis++)
                            {
                                for (int ipw_at_k=0; ipw_at_k<V_left[ispin][ik][0].size(); ipw_at_k++)
                                {
                                    V_tmp_left[iband](ipw_at_k) += eigenvectors_left(ibasis, iband) * V_left[ispin][ik][ibasis](ipw_at_k);
                                }
                            }
                        }
                    } // biortho_david_

                    if (!biortho_david_)
                    {
                        if (!(is_bitc && iter_inner==(num_refresh_david_-1))) // Not the "BITC final loop" (where no orthonormalization is needed)
                        {
                            // We orhonormalize V_tmp rather than "eigenvectors" for numerical stability
                            for (int iband=0; iband<num_bands_tc[ispin]; iband++)
                            {
                                if (!Gram_Schmidt(V_tmp, iband, V_tmp[iband]))
                                {
                                    error_messages::stop("Linear independence breaks in block-Davidson algorithm.");
                                }
                            }

                            // update bloch_states.phik for the TC final loop
                            if (iter_inner==num_refresh_david_-1) { bloch_states.reset_phik(ispin, ik, V_tmp, "right", mixes_density_matrix_, method.calc_mode()); }
                        }
                        else // BITC final loop: update bloch_states.phik & phik_left
                        {
                            // Note: V^{-1}*V = id. <=> eigenvectors_left(*,i).conj() * eigenvectors(*,j) = delta_{i,j}
                            // i.e. left and right eigenvectors are already bi-orthonormalized. 
                            // Also right eigevectors are normalized by Eigen implementation of compute().

                            // update bloch_states.phik for the BITC final loop
                            bloch_states.reset_phik(ispin, ik, V_tmp, "right", mixes_density_matrix_, method.calc_mode());

                            for (int iband=0; iband<num_bands_tc[ispin]; iband++)
                            {
                                V_tmp[iband] = Eigen::VectorXcd::Zero(V[ispin][ik][0].size());
                                for (int ibasis=0; ibasis<ndim_total; ibasis++)
                                {
                                    for (int ipw_at_k=0; ipw_at_k<V[ispin][ik][0].size(); ipw_at_k++)
                                    {
                                        V_tmp[iband](ipw_at_k) += eigenvectors_left(ibasis, iband) * V[ispin][ik][ibasis](ipw_at_k);
                                    }
                                }
                            }
                            bloch_states.reset_phik(ispin, ik, V_tmp, "left", mixes_density_matrix_, method.calc_mode());
                        }

                        if (iter_inner==(num_refresh_david_-1)) // final loop
                        {
                            // update bloch_states.eigenvalues
                            bloch_states.reset_eigenvalues(ispin, ik, eigenvalues, method.calc_mode());
                            
                            if (method.calc_mode()=="SCF")
                            {
                                // update bloch_states.energy_*body (required for total-energy calculation)
                                Eigen::MatrixXcd energy_1body, energy_2body, energy_3body;
                                if (!is_bitc)
                                {
                                    energy_1body = eigenvectors.adjoint() * H1_subspace[ispin][ik] * eigenvectors;
                                    energy_2body = eigenvectors.adjoint() * H2_subspace[ispin][ik] * eigenvectors;
                                    if (uses_3body) { energy_3body = eigenvectors.adjoint() * H3_subspace[ispin][ik] * eigenvectors; }
                                }
                                else // BITC
                                {
                                    energy_1body = eigenvectors_left.adjoint() * H1_subspace[ispin][ik] * eigenvectors;
                                    energy_2body = eigenvectors_left.adjoint() * H2_subspace[ispin][ik] * eigenvectors;
                                    energy_3body = eigenvectors_left.adjoint() * H3_subspace[ispin][ik] * eigenvectors;
                                }
                                total_energy.reset_energies(ispin, ik, energy_1body, energy_2body, energy_3body, uses_3body);
                            } // SCF
                        } 
                        else // not the final loop: update V by V_tmp
                        {
                            for (int iband=0; iband<num_bands_tc[ispin]; iband++)
                            {
                                V[ispin][ik][iband] = V_tmp[iband];
                            }
                        } // final loop or not
                    }
                    else // biortho_david_
                    {
                        // Note: X^{-1}*X = id. <=> eigenvectors_left(*,i).conj() * eigenvectors(*,j) = delta_{i,j}
                        // i.e. left and right eigenvectors are already bi-orthogonal.
                        // Also right eigevectors are normalized by Eigen implementation of compute().

                        // However, note that V_tmp is NOT normalized even though eigenvectors() is normalized
                        // because V is no longer orthnormalized when uses_biortho_david==true

                        // Therefore, we need to multiply a constant!
                        for (int iband=0; iband<num_bands_tc[ispin]; iband++)
                        {
                            double norm = V_tmp[iband].norm();
                            if (norm < 1e-10) { error_messages::stop("Linear independence breaks in block-Davidson algorithm."); }
                            V_tmp[iband] /= norm;
                            V_tmp_left[iband] *= norm; // to keep <iband | iband> = 1
                        }

                        if (iter_inner==(num_refresh_david_-1)) // final loop
                        {
                            // update bloch_states.phik for the BITC final loop
                            bloch_states.reset_phik(ispin, ik, V_tmp, "right", mixes_density_matrix_, method.calc_mode());
                            bloch_states.reset_phik(ispin, ik, V_tmp_left, "left", mixes_density_matrix_, method.calc_mode());

                            // update bloch_states.eigenvalues
                            bloch_states.reset_eigenvalues(ispin, ik, eigenvalues, method.calc_mode());
                            
                            if (method.calc_mode()=="SCF")
                            {
                                // update bloch_states.energy_*body (required for total-energy calculation)
                                Eigen::MatrixXcd energy_1body, energy_2body, energy_3body;
                                energy_1body = eigenvectors_left.adjoint() * H1_subspace[ispin][ik] * eigenvectors;
                                energy_2body = eigenvectors_left.adjoint() * H2_subspace[ispin][ik] * eigenvectors;
                                energy_3body = eigenvectors_left.adjoint() * H3_subspace[ispin][ik] * eigenvectors;
                                total_energy.reset_energies(ispin, ik, energy_1body, energy_2body, energy_3body, uses_3body);
                            } // SCF
                        } 
                        else // not the final loop: update V by V_tmp
                        {
                            for (int iband=0; iband<num_bands_tc[ispin]; iband++)
                            {
                                V[ispin][ik][iband] = V_tmp[iband];
                                V_left[ispin][ik][iband] = V_tmp_left[iband];
                            }
                        } // final loop or not
                    } // biortho_david
                } // ik
            } // ispin
            my_clock.print_time_from_save(ost, "diagonalizing subspace Hamiltonian");
        } // am_i_mpi_rank0

        if (iter_inner!=(num_refresh_david_-1)) { // not the final loop
            // Bcast new trial vectors (for the next "iter_inner" loop)
            for (int ispin=0; ispin<num_independent_spins; ispin++)
            {
                for (int ik=0; ik<num_irreducible_kpoints; ik++)
                {
                    for (int iband=0; iband<num_bands_tc[ispin]; iband++)
                    {
                        MPI_Bcast(V[ispin][ik][iband].data(), V[ispin][ik][iband].size(),
                                  MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
                        if (biortho_david_)
                        {
                            MPI_Bcast(V_left[ispin][ik][iband].data(), V_left[ispin][ik][iband].size(),
                                      MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
                        }
                    }

                    // deallocate
                    V[ispin][ik].resize(num_bands_tc[ispin]);
                    if (biortho_david_) { V_left[ispin][ik].resize(num_bands_tc[ispin]); }
                }
            }
        }
        else // final loop: update bloch_states
        {
            bloch_states.bcast_phik(is_bitc, mixes_density_matrix_, method.calc_mode(), am_i_mpi_rank0);
            bloch_states.bcast_eigenvalues(method.calc_mode());
            if (method.calc_mode()=="SCF") { total_energy.bcast_energies(uses_3body); }
        }
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "diagonalizing subspace Hamiltonian (Bcast)"); }
    } // iter_inner
}
