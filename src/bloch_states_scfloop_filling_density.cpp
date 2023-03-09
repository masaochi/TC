// [class BlochStates]
// wave functions, eigenvalues, filling, num. of bands, etc.

#include "include/header.hpp"

// public

// set fermi_energy_ & filling_ & filling_old_ & num_occupied_bands_ 
void BlochStates::set_filling(const Spin &spin, const Kpoints &kpoints, 
                              const bool sets_filling_old, const double &mixing_beta,
                              const bool am_i_mpi_rank0, std::ostream *ost)
{
    const int num_independent_spins = filling_.size();
    const int num_irreducible_kpoints_scf = filling_[0].size();
    const double spin_factor = (num_independent_spins==1 && !spin.is_spinor()) ? 2.0 : 1.0; // 2 for no-spin

    if (am_i_mpi_rank0) {  *ost << "  bloch_states.set_filling()" << std::endl; }

    // set filling_old_
    if (sets_filling_old)
    {
        for (int ispin=0; ispin<num_independent_spins; ispin++)
        {
            for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
            {
                filling_old_[ispin][ik].resize(num_occupied_bands_[ispin][ik]);
                for (int iband=0; iband<num_occupied_bands_[ispin][ik]; iband++)
                {
                    filling_old_[ispin][ik][iband] = filling_[ispin][ik][iband];
                }
            }
        }
    }

    if (am_i_mpi_rank0)
    {
        if (kpoints.smearing_mode()=="fixed")
        {
            // find the valence-band top assuming the insulating band structure

            const int num_occupied_bands = spin.is_spinor() ?
                std::round(num_electrons_) : std::round(num_electrons_/2.0);

            fermi_energy_ = -1e5;
            for (int ispin=0; ispin<num_independent_spins; ispin++)
            {
                if (num_occupied_bands > num_bands_scf_[ispin]) { error_messages::stop("num_bands_tc too small in set_fermi_energy()"); }
                for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
                {
                    if (kpoints.kweight_scf()[ik] < 1e-8) // zero k-weight (e.g. band k-point)
                    {
                        for (int iband=0; iband<num_bands_scf_[ispin]; iband++)
                        {
                            filling_[ispin][ik][iband] = 0.0; 
                        }
                    }
                    else // non-zero k-weight
                    {
                        for (int iband=0; iband<num_occupied_bands; iband++)
                        {
                            filling_[ispin][ik][iband] = spin_factor;
                            
                            double orbital_energy = eigenvalues_scf_[ispin][ik][iband].real();
                            if (fermi_energy_ < orbital_energy) { fermi_energy_ = orbital_energy; }
                        }
                    } // zero or non-zero k-weight
                }
            }

            // to check consistency
            double conduction_band_minimum = 1e5;
            for (int ispin=0; ispin<num_independent_spins; ispin++)
            {
                for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
                {
                    if (kpoints.kweight_scf()[ik] > 1e-8) // non-zero k-weight (e.g. band k-point)
                    {
                        for (int iband=num_occupied_bands; iband<num_bands_scf_[ispin]; iband++)
                        {
                            double orbital_energy = eigenvalues_scf_[ispin][ik][iband].real();
                            if (conduction_band_minimum > orbital_energy) { conduction_band_minimum = orbital_energy; } 
                        }
                    } // non-zero k-weight
                }
            }
            if (conduction_band_minimum < fermi_energy_) 
            {
                *ost << "   CAUTION!! the conduction band minimum = " << conduction_band_minimum << " Hartree";
                *ost << " is lower than the valence band maximum = " << fermi_energy_ << " Hartree" << std::endl;
                *ost << "   while the fixed-occupation is applied here." << std::endl;
            }
        }
        else if (kpoints.smearing_mode()=="gaussian")
        {
            double upper_bound = 1e5;
            double lower_bound = -1e5;
            for (int iter=0; iter<60; iter++)
            {
                fermi_energy_ = (upper_bound + lower_bound)/2.0;

                double tmp_num_electrons = 0.0;
                for (int ispin=0; ispin<num_independent_spins; ispin++)
                {
                    for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
                    {
                        if (kpoints.kweight_scf()[ik] < 1e-8) // zero k-weight (e.g. band k-point)
                        {
                            for (int iband=0; iband<num_bands_scf_[ispin]; iband++)
                            {
                                filling_[ispin][ik][iband] = 0.0; 
                            }
                        }
                        else // non-zero k-weight
                        {
                            for (int iband=0; iband<num_bands_scf_[ispin]; iband++)
                            {
                                filling_[ispin][ik][iband] 
                                    = kpoints.return_band_filling(spin,
                                                                  eigenvalues_scf_[ispin][ik][iband].real(), fermi_energy_,
                                                                  num_electrons_, iband);

                                // "*kweight_scf": num. of symmetrically-equivalent k-points * spin-weight (spin_factor)
                                // "/spin_factor": prevent double-counting of spin-weight
                                tmp_num_electrons += 
                                    filling_[ispin][ik][iband]*(kpoints.kweight_scf()[ik]/spin_factor);
                            }
                        } // zero or non-zero k-weight
                    }
                }
                if (tmp_num_electrons > num_electrons_)
                {
                    upper_bound = fermi_energy_;
                }
                else
                {
                    lower_bound = fermi_energy_;
                }
            } // iter
        } // smearing_mode

        // print fermi_energy & filling
        *ost << "   Fermi energy: " << fermi_energy_ * Ht_in_eV << " eV" << std::endl;
        *ost << std::endl;
        for (int ispin=0; ispin<num_independent_spins; ispin++)
        {
            *ost << "   Spin " << ispin << std::endl;
            for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
            {
                *ost << "   k-vector " << kpoints.kvectors_scf()[ik][0].transpose() << std::endl;
                *ost << "   (band index, energy (eV), filling)" << std::endl;
                for (int iband=0; iband<num_bands_scf_[ispin]; iband++)
                {
                    *ost << "   " << iband << " " << eigenvalues_scf_[ispin][ik][iband].real() * Ht_in_eV;
                    *ost << " " << filling_[ispin][ik][iband] << std::endl;
                }
                *ost << std::endl;
            }
        }

        // check num. of electrons
        double tmp_num_elect = 0.0;
        for (int ispin=0; ispin<num_independent_spins; ispin++)
        {
            for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
            {
                for (int iband=0; iband<num_bands_scf_[ispin]; iband++)
                {
                 tmp_num_elect +=
                     filling_[ispin][ik][iband]*(kpoints.kweight_scf()[ik]/spin_factor);
                }
            }
        }

        if (std::abs(tmp_num_elect - num_electrons_) > 1e-8) { error_messages::stop("set Fermi energy failed (try: increase nbands_tc and/or change smearing_mode)."); }

        if (std::abs(fermi_energy_ - 1e5) < 1e-8)
        {
            *ost << "   Note: all bands are occupied and thus fermi_energy becomes a bit peculiar value." << std::endl;
            *ost << "   Please make sure that num_bands is sufficiently large... (we recommend increasing num_bands_tc)" << std::endl;
        }

        // set num_occupied_bands_
        for (int ispin=0; ispin<num_independent_spins; ispin++)
        {
            for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
            {
                num_occupied_bands_[ispin][ik] = 0;
                for (int iband=0; iband<num_bands_scf_[ispin]; iband++)
                {
                    // For gaussian smearing, bands with "energy > fermi_energy + (around) 4*smearing_wdith"
                    // are regarded as "unoccupied".
                    if (filling_[ispin][ik][iband] > 1e-8) { num_occupied_bands_[ispin][ik]++; }
                }
            }
        }
    } // am_i_mpi_rank0

    // bcast    
    MPI_Bcast(&fermi_energy_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        MPI_Bcast(&num_occupied_bands_[ispin][0], num_irreducible_kpoints_scf,
                  MPI_INT, 0, MPI_COMM_WORLD);
        for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
        {
            MPI_Bcast(&filling_[ispin][ik][0], filling_[ispin][ik].size(),
                      MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
}

// called only when diagonalization.mixes_density_matrix_==true && not the first loop.
void BlochStates::mix_density_matrix(const double &mixing_beta)
{
    assert(filling_.size() == filling_old_.size());
    assert(filling_[0].size() == filling_old_[0].size());

    for (int ispin=0; ispin<filling_.size(); ispin++)
    {
        for (int ik=0; ik<filling_[ispin].size(); ik++)
        {
            for (int iband=0; iband<filling_[ispin][ik].size(); iband++)
            {
                filling_[ispin][ik][iband] *= mixing_beta;
            }
            for (int iband=0; iband<filling_old_[ispin][ik].size(); iband++)
            {
                filling_old_[ispin][ik][iband] *= (1.0 - mixing_beta);
            }
        }
    }
}

void BlochStates::recover_filling_for_density_matrix(const double &mixing_beta)
{
    assert(filling_.size() == filling_old_.size());
    assert(filling_[0].size() == filling_old_[0].size());

    for (int ispin=0; ispin<filling_.size(); ispin++)
    {
        for (int ik=0; ik<filling_[ispin].size(); ik++)
        {
            for (int iband=0; iband<filling_[ispin][ik].size(); iband++)
            {
                filling_[ispin][ik][iband] /= mixing_beta;
            }
            for (int iband=0; iband<filling_old_[ispin][ik].size(); iband++)
            {
                filling_old_[ispin][ik][iband] = 0.0;
            }
        }
    }
}

// non-collinear calculation is not supported
void BlochStates::set_density(const CrystalStructure &crystal_structure,
                              const Kpoints &kpoints, PlaneWaveBasis &plane_wave_basis,
                              const bool &mixes_density_matrix, const double &mixing_beta,
                              const bool is_first_iter, const bool is_bitc,
                              const bool am_i_mpi_rank0, std::ostream *ost)
{
    const int num_independent_spins = filling_.size();
    const int num_irreducible_kpoints_scf = filling_[0].size();
    const int npw = plane_wave_basis.size_FFT_grid();

    if (is_first_iter) // allocate
    {
        density_.resize(2);
        for (int ispin=0; ispin<2; ispin++) { density_[ispin].resize(npw); }
    }

    // calculate density
    if (am_i_mpi_rank0)
    {
        *ost << "  bloch_states.set_density()" << std::endl;

        density_difference_ = 0.0;
        Eigen::VectorXcd tmp_density, orbital, orbital_left;
        Eigen::VectorXcd tmp_density_sub(npw);
        for (int ispin=0; ispin<num_independent_spins; ispin++) 
        {
            tmp_density  = Eigen::VectorXcd::Zero(npw);
            for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
            {
                for (int isym=0; isym<kpoints.kvectors_scf()[ik].size(); isym++)
                {
                    const int nbands_old = filling_old_[ispin][ik].size();
                    for (int iband=-nbands_old; iband<num_occupied_bands_[ispin][ik]; iband++) // negative = orbitals in the previous SCF loop
                    {
                        if (iband>=0) 
                        {
                            plane_wave_basis.get_orbital_FFTgrid(ispin, ik, isym,
                                                                 kpoints.is_time_reversal_used_at_k()[ik][isym],
                                                                 phik_scf_[ispin][ik][iband][0], orbital, "SCF");
                        } 
                        else
                        {
                            plane_wave_basis.get_orbital_FFTgrid(ispin, ik, isym,
                                                                 kpoints.is_time_reversal_used_at_k()[ik][isym],
                                                                 phik_scf_old_[ispin][ik][-1-iband][0], orbital, "SCF");
                        }
                        plane_wave_basis.FFT_backward(orbital, orbital);
                        if (!is_bitc)
                        {
                            tmp_density_sub = orbital.array().abs2();
                        }
                        else // BITC
                        {
                            if (iband>=0)
                            {
                                plane_wave_basis.get_orbital_FFTgrid(ispin, ik, isym,
                                                                     kpoints.is_time_reversal_used_at_k()[ik][isym],
                                                                     phik_left_scf_[ispin][ik][iband][0], orbital_left, "SCF");
                            }
                            else
                            {
                                plane_wave_basis.get_orbital_FFTgrid(ispin, ik, isym,
                                                                     kpoints.is_time_reversal_used_at_k()[ik][isym],
                                                                     phik_left_scf_old_[ispin][ik][-1-iband][0], orbital_left, "SCF");
                            }
                            plane_wave_basis.FFT_backward(orbital_left, orbital_left);
                            tmp_density_sub = orbital_left.array().conjugate() * orbital.array();
                        } // is_bitc

                        if (std::abs(tmp_density_sub.sum()/static_cast<double>(npw) - 1.0) > 1e-4) 
                        {
                            std::cout << ispin << " " << ik << " " << isym << " " << iband << std::endl;
                            std::cout << tmp_density_sub.sum()/static_cast<double>(npw) << ": Not equal to 1" << std::endl;
                            error_messages::stop("normalization breaks."); 
                        }

                        if (iband>=0)
                        {
                            tmp_density += filling_[ispin][ik][iband]*tmp_density_sub;
                        }
                        else
                        {
                            tmp_density += filling_old_[ispin][ik][-1-iband]*tmp_density_sub;
                        }
                    } // iband
                } // isym
            } // ik
            if (num_independent_spins==1) { tmp_density /= 2.0; } // prevent double-counting of the spin factor ("2" is already included in filling_)
            if (is_first_iter) // first loop
            {
                density_[ispin] = tmp_density;
            }
            else // not the first loop
            {
                if (!mixes_density_matrix) { tmp_density = (1.0-mixing_beta)*density_[ispin] + mixing_beta*tmp_density; } // density mixing (new density)
                density_difference_ += (tmp_density - density_[ispin]).array().abs().sum(); // L1 norm: \int |rho - rho_old|
                density_[ispin] = tmp_density;
            }
        } // ispin
        if (num_independent_spins==1) 
        {
            density_[1] = density_[0];
            density_difference_ *= 2;
        }
        density_difference_ /= static_cast<double>(npw*kpoints.num_kpoints());

        // num. of electrons for each spin (per unit cell)
        Complex num_electrons_up = density_[0].sum()/static_cast<double>(npw*kpoints.num_kpoints());
        Complex num_electrons_dn = density_[1].sum()/static_cast<double>(npw*kpoints.num_kpoints());
        *ost << "   num. of electrons (spin up) = " << num_electrons_up;
        *ost << " , num. of electrons (spin down) = " << num_electrons_dn << std::endl;
        *ost << "   total magnetization = " << num_electrons_up - num_electrons_dn << std::endl;

        integrate_local_density(crystal_structure, kpoints, plane_wave_basis, ost); // for showing local num. of electrons and local magnetic moment

        // density can be inconsistent when the normalization of the orbitals breaks (by some numerical problems in diagonalization or Gramd-Schmidt orthonormalization)
        if (std::abs(num_electrons_up + num_electrons_dn - num_electrons_) > 1e-4 ) { error_messages::stop("inconsistent num. of electrons in set_density()."); }

        if (!is_first_iter) { *ost << "   density difference (L1 norm) from the previous loop = " << density_difference_ << std::endl; }
    } // am_i_mpi_rank0

    for (int ispin=0; ispin<2; ispin++)
    {
        MPI_Bcast(density_[ispin].data(), density_[ispin].size(),
                  MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(&density_difference_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (am_i_mpi_rank0) { *ost << std::endl; }
}

// private

void BlochStates::integrate_local_density(const CrystalStructure &crystal_structure,
                                          const Kpoints &kpoints, const PlaneWaveBasis &plane_wave_basis,
                                          std::ostream *ost)
{

    const int num_independent_spins = filling_.size();
    const int npw = plane_wave_basis.size_FFT_grid();
    const int npw1 = plane_wave_basis.size_FFT_grid_vec()[0];
    const int npw2 = plane_wave_basis.size_FFT_grid_vec()[1];
    const int npw3 = plane_wave_basis.size_FFT_grid_vec()[2];

    const int nradius = 40; // Radius-mesh size for integrating the local quantities
    const double dradius = 0.1; // dr for the radius mesh (bohr)

    const int num_atoms = crystal_structure.num_atoms();
    // (0,0:2) = a1, (1,0:2) = a2, (2,0:2) = a3 in cartesian coordinate (in bohr)
    const Eigen::Matrix3d &lattice_vectors = crystal_structure.lattice_vectors();
    const std::vector<Eigen::Vector3d> &atomic_position_cartesian
        = crystal_structure.atomic_position_cartesian(); // cartesian coordinate of each atom (in bohr)

    // density_around_atom[ispin][iatom][iradius]
    std::vector<std::vector<std::vector<Complex> > > density_around_atom(2);
    for (int ispin=0; ispin<2; ispin++) 
    { 
        density_around_atom[ispin].resize(num_atoms); 
        for (int iatom=0; iatom<num_atoms; iatom++)
        {
            density_around_atom[ispin][iatom].resize(nradius);
        }
    }

    Eigen::Vector3d Rvector, Rvector_diff, Rvector_diff_trans;
    for (int ipw=0; ipw<npw; ipw++)
    {
        int i1 = ipw%npw1;
        int i2 = (ipw/npw1)%npw2;
        int i3 = ipw/(npw1*npw2);

        for (int idim=0; idim<3; idim++)
        {
            Rvector(idim) = static_cast<double>(i1)/static_cast<double>(npw1)*lattice_vectors(0,idim)
                + static_cast<double>(i2)/static_cast<double>(npw2)*lattice_vectors(1,idim)
                + static_cast<double>(i3)/static_cast<double>(npw3)*lattice_vectors(2,idim);
        }
        for (int iatom=0; iatom<num_atoms; iatom++)
        {
            Rvector_diff = Rvector - atomic_position_cartesian[iatom];
            for (int ii1=-3; ii1<4; ii1++)
            {
                for (int ii2=-3; ii2<4; ii2++)
                {
                    for (int ii3=-3; ii3<4; ii3++)
                    {
                        for (int idim=0; idim<3; idim++)
                        {
                            Rvector_diff_trans(idim) = Rvector_diff(idim) 
                                + ii1*lattice_vectors(0,idim)
                                + ii2*lattice_vectors(1,idim)
                                + ii3*lattice_vectors(2,idim);
                        }
                        double norm_diff = Rvector_diff_trans.norm();

                        for (int iradius=0; iradius<nradius; iradius++)
                        {
                            double radius = (iradius+1)*dradius;
                            if (norm_diff < radius) 
                            {
                                for (int ispin=0; ispin<2; ispin++)
                                {
                                    density_around_atom[ispin][iatom][iradius]
                                        += density_[ispin](ipw);
                                }
                            } // if
                        } // iradius
                    } // iix
                } // iiy
            } // iiz
        } // iatom
    } // ipw

    for (int ispin=0; ispin<2; ispin++)
    {
        for (int iatom=0; iatom<num_atoms; iatom++)
        {
            for (int iradius=0; iradius<nradius; iradius++)
            {
                density_around_atom[ispin][iatom][iradius] /= static_cast<double>(npw*kpoints.num_kpoints());
            }
        }
    }                

    *ost << "   num. of electrons inside a sphere with a radius r around each atom" << std::endl;
    for (int iatom=0; iatom<num_atoms; iatom++)
    {
        *ost << "    Atom " << iatom << ": " << std::endl;
        *ost << "    r (Bohr),  up,  down,  up - down,  up + down" << std::endl;
        for (int iradius=0; iradius<nradius; iradius++)
        {
            double radius = (iradius+1)*dradius;
            double up_density = density_around_atom[0][iatom][iradius].real();
            double dn_density = density_around_atom[1][iatom][iradius].real();
            *ost << "    " << radius << "  " << up_density << "  " << dn_density 
                 << "  " << up_density - dn_density << "  " << up_density + dn_density << std::endl;
        } // iradius
    } // iatom
}
