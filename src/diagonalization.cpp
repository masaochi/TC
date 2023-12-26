// [class Diagonalization]
// scf_loop, convergence check, control paramters, etc.

#include "include/header.hpp"

void Diagonalization::set(const bool restarts,
                          const double &energy_tolerance, const double &charge_tolerance,
                          const double &force_tolerance,
                          const int &max_num_iterations, const int &max_num_ionic_steps,
                          const bool &mixes_density_matrix, const double &mixing_beta,
                          const int &num_refresh_david, const int &max_num_blocks_david,
                          const bool &biortho_david, const bool &dumps_pwfn, const bool &reads_crystal_structure,
                          const std::string &calc_mode, const std::string &calc_method, const bool is_heg)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");

    restarts_ = restarts;

    if (energy_tolerance > -1e-8) 
    {
        energy_tolerance_ = energy_tolerance;
    }
    else
    {
        error_messages::inappropriate_argument("energy_tolerance", energy_tolerance,
                                               "should not be negative");
    }

    if (charge_tolerance > -1e-8)
    {
        charge_tolerance_ = charge_tolerance;
    }
    else
    {
        error_messages::inappropriate_argument("charge_tolerance", charge_tolerance,
                                               "should not be negative");
    }

    if (force_tolerance > -1e-8)
    {
        force_tolerance_ = force_tolerance;
    }
    else
    {
        error_messages::inappropriate_argument("force_tolerance", force_tolerance,
                                               "should not be negative");
    }

    if (max_num_iterations >= 0) 
    {
        max_num_iterations_ = max_num_iterations;
    }
    else // use a default value (see include/diagonalization.hpp)
    {
        if (calc_mode=="SCF")
        {
            max_num_iterations_ = 30;
        }
        else if (calc_mode=="BAND")
        {
            max_num_iterations_ = 15;
        }
//        error_messages::inappropriate_argument("max_num_iterations", max_num_iterations,
//                                               "should not be negative");
    }

    max_num_ionic_steps_ = max_num_ionic_steps;
    if (max_num_ionic_steps < 0) 
    { 
        error_messages::inappropriate_argument("max_num_ionic_steps", max_num_ionic_steps, "should not be negative");
    }
    if (calc_mode=="BAND" && max_num_ionic_steps > 0)
    {
        error_messages::inappropriate_argument("max_num_ionic_steps", max_num_ionic_steps, "should be zero for calc_mode == BAND (use calc_mode == SCF for structural opt.)");
    }
    if (is_heg && max_num_ionic_steps > 0)
    {
        error_messages::inappropriate_argument("max_num_ionic_steps", max_num_ionic_steps, "should be zero for is_heg == true since ionic potentials are neglected for HEG calc.");
    }
    if (calc_method=="TC" && max_num_ionic_steps > 0)
    {
        error_messages::inappropriate_argument("max_num_ionic_steps", max_num_ionic_steps, "should be zero for calc_method == TC because Hellmann-Feynman theorem does not hold. Use calc_method == BITC instead.");
    }

    mixes_density_matrix_ = mixes_density_matrix;

    if (mixing_beta > 1e-8) 
    {
        mixing_beta_ = mixing_beta;
    }
    else
    {
        error_messages::inappropriate_argument("mixing_beta", mixing_beta,
                                               "should be positive");
    }

    if (num_refresh_david >= 1)
    {
        num_refresh_david_ = num_refresh_david;
    }
    else
    {
        error_messages::inappropriate_argument("max_refresh_david", num_refresh_david,
                                               "should >=1");
    }

    if (max_num_blocks_david >= 2) 
    {
        max_num_blocks_david_ = max_num_blocks_david; 
    } 
    else
    {
        error_messages::inappropriate_argument("max_num_blocks_david", max_num_blocks_david,
                                               "should >=2");
    }

    if (biortho_david && calc_method!="BITC")
    {
        error_messages::inappropriate_argument("biortho_david", biortho_david,
                                               "should be false for calc_method!=BITC");
    }
    biortho_david_ = biortho_david;

    dumps_pwfn_ = dumps_pwfn;

    if (reads_crystal_structure && max_num_ionic_steps==0)
    {
        error_messages::inappropriate_argument("reads_crystal_structure", reads_crystal_structure,
                                               "should be false when max_num_ionic_steps is zero");
    }
    else
    {
        reads_crystal_structure_ = reads_crystal_structure;
    }
}

void Diagonalization::bcast(const bool am_i_mpi_rank0)
{
    MPI_Bcast(&restarts_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&energy_tolerance_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&charge_tolerance_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&force_tolerance_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_num_iterations_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_num_ionic_steps_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mixes_density_matrix_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mixing_beta_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_refresh_david_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_num_blocks_david_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&biortho_david_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dumps_pwfn_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&reads_crystal_structure_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
}

void Diagonalization::scf(Parallelization &parallelization,
                          const FileNames &file_names, MyClock &my_clock,
                          const Method &method, CrystalStructure &crystal_structure,
                          const Symmetry &symmetry, Potentials &potentials,
                          const Spin &spin, const Kpoints &kpoints,
                          PlaneWaveBasis &plane_wave_basis, 
                          BlochStates &bloch_states, TotalEnergy &total_energy,
                          std::ostream *ost) const
{
    const bool am_i_mpi_rank0 = parallelization.am_i_mpi_rank0();
    const bool is_bitc = method.calc_method()=="BITC" ? true : false;
    const bool uses_3body = (method.calc_method()=="TC" || method.calc_method()=="BITC") ? true : false;
    const bool is_heg = potentials.is_heg();
    assert(!(is_heg && max_num_ionic_steps_!=0)); // cannot perform structural opt. for homogeneous-electron-gas calc.
    assert(!(method.calc_method()=="TC" && max_num_ionic_steps_!=0)); // cannot perform structural opt. for TC calc. (Hellmann-Feynman theorem breaks)

    bool is_force_converged = false;
    for (int istep_ionic=0; istep_ionic<max_num_ionic_steps_+1; istep_ionic++) // ionic step
    {
        if (am_i_mpi_rank0) { *ost << " Start SCF calculation!" << std::endl; }

        bool is_energy_converged, is_energy_converged_prev, is_charge_converged, is_charge_converged_prev, is_scf_converged;
        is_energy_converged = is_energy_converged_prev = is_charge_converged = is_charge_converged_prev = is_scf_converged = false;
        for (int iter=0; iter<max_num_iterations_; iter++)
        {
            bool is_first_iter = iter==0 ? true : false;
            
            // Preparation for SCF calc.
            bloch_states.set_filling(spin, kpoints, (!is_first_iter && mixes_density_matrix_),
                                     am_i_mpi_rank0, ost); // set filling_old when (!is_fisrt_iter && mixes_...)==true
            if (!is_first_iter && mixes_density_matrix_) { bloch_states.mix_density_matrix(mixing_beta_); }
            bloch_states.set_density(crystal_structure, kpoints, plane_wave_basis, mixes_density_matrix_, 
                                     mixing_beta_, is_first_iter, is_bitc, am_i_mpi_rank0, ost);
            if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "preparation for the next SCF loop. (set filling and density)"); }            

            // convergence check and force calc.
            if (!is_first_iter) 
            {
                total_energy.calc_total_energy(spin, kpoints,
                                               bloch_states.filling(), bloch_states.eigenvalues_scf(),
                                               bloch_states.fermi_energy(), uses_3body, am_i_mpi_rank0, ost);
                is_energy_converged =
                    (std::abs(total_energy.total_energy_difference()) < energy_tolerance_);

                is_charge_converged = (bloch_states.density_difference() < charge_tolerance_);
                
                is_scf_converged = is_energy_converged && is_energy_converged_prev &&
                    is_charge_converged && is_charge_converged_prev;

                // compute force
                if (!is_heg && (method.calc_method() == "HF" || method.calc_method() == "BITC")) // H-F theorem does not hold for TC
                {
                    total_energy.set_force(calc_hamiltonian::force(parallelization, method,
                                                                   crystal_structure,
                                                                   potentials, spin, kpoints,
                                                                   plane_wave_basis,
                                                                   bloch_states, total_energy, ost));
                    if (am_i_mpi_rank0) { total_energy.print_force(is_scf_converged, ost); } // use is_scf_converged
 
                    is_force_converged = (total_energy.maximum_force_in_eV_ang() < force_tolerance_);
                } // if (calc_method)

                if (is_scf_converged)
                {
                    if (am_i_mpi_rank0) { *ost << "  convergence is achieved in SCF!" << std::endl; }
                    if (am_i_mpi_rank0) { total_energy.print_total_energy_final_loop(kpoints, ost); }

                    break; // go out from an "iter" loop
                }
            } // !is_first_iter
            is_energy_converged_prev = is_energy_converged;
            is_charge_converged_prev = is_charge_converged;
            
            if (am_i_mpi_rank0) { *ost << " Iteration " << iter+1 << std::endl; }
            
            // set Hamoiltonian and diagonalize it
            block_davidson(parallelization, my_clock, method,
                           crystal_structure, symmetry,
                           potentials, spin, kpoints,
                           plane_wave_basis, 
                           bloch_states, total_energy, ost);
            
            if (!is_first_iter && mixes_density_matrix_) { bloch_states.recover_filling_for_density_matrix(mixing_beta_); } // recover filling
            
            if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "diagonalization"); }
            if (am_i_mpi_rank0) { io_tc_files::dump_eigen(file_names, method, "SCF", bloch_states, ost); } // output for each iteration
            if (am_i_mpi_rank0 && !is_heg && dumps_pwfn_) // output for CASINO
            {
                io_qmc_files::dump_pwfn(file_names, crystal_structure, kpoints, plane_wave_basis, bloch_states, total_energy, ost);
            }
            if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "dumping wave functions etc."); }
            if (am_i_mpi_rank0) { *ost << std::endl; }
        } // iter

        if (!is_scf_converged) // When the convergence is NOT achieved
        {
            bloch_states.set_filling(spin, kpoints, false, am_i_mpi_rank0, ost); // no need to mix the density matrix
            bloch_states.set_density(crystal_structure, kpoints, plane_wave_basis, mixes_density_matrix_, 
                                     mixing_beta_, false, is_bitc, am_i_mpi_rank0, ost); // needed for calculating the force
            total_energy.calc_total_energy(spin, kpoints,
                                           bloch_states.filling(), bloch_states.eigenvalues_scf(),
                                           bloch_states.fermi_energy(), uses_3body, am_i_mpi_rank0, ost);
            
            if (!is_heg && (method.calc_method() == "HF" || method.calc_method() == "BITC")) // H-F theorem does not hold for TC
            {
                total_energy.set_force(calc_hamiltonian::force(parallelization, method,
                                                               crystal_structure,
                                                               potentials, spin, kpoints,
                                                               plane_wave_basis,
                                                               bloch_states, total_energy, ost));
                if (am_i_mpi_rank0) { total_energy.print_force(is_scf_converged, ost); } // use is_scf_converged
                
                is_force_converged = (total_energy.maximum_force_in_eV_ang() < force_tolerance_);
            }
            
            if (am_i_mpi_rank0) { *ost << " CAUTION!! convergence is not achieved in SCF..." << std::endl; }
        }

        // move atoms
        if (max_num_ionic_steps_!=0) // perform structural opt.
        {
            if (is_force_converged)
            {
                if (am_i_mpi_rank0) { *ost << "  Force is sufficiently small!" << std::endl; }
                return;   
            }
            else
            {
                if (am_i_mpi_rank0) { *ost << "  Force is NOT sufficiently small!" << std::endl; }
                if (istep_ionic<max_num_ionic_steps_) // not the final loop
                {
                    structural_optimization(method, crystal_structure, symmetry,
                                            total_energy, istep_ionic, am_i_mpi_rank0, ost);
                }
                total_energy.reset_ewald_energy(crystal_structure, is_heg, am_i_mpi_rank0);
                potentials.reset_Vpp_local(crystal_structure,
                                           plane_wave_basis,
                                           am_i_mpi_rank0);
                if (am_i_mpi_rank0) { io_tc_files::dump_crystal_structure(file_names, crystal_structure, ost); }
            } // if (is_force_converged)
        } // if (max_num_ionic_steps_!=0)
    } // istep_ionic
}

void Diagonalization::band(Parallelization &parallelization,
                           const FileNames &file_names, MyClock &my_clock,
                           const Method &method, const CrystalStructure &crystal_structure,
                           const Symmetry &symmetry, Potentials &potentials,
                           const Spin &spin, const Kpoints &kpoints,
                           PlaneWaveBasis &plane_wave_basis, 
                           BlochStates &bloch_states, TotalEnergy &total_energy,
                           std::ostream *ost) const
{
    const bool am_i_mpi_rank0 = parallelization.am_i_mpi_rank0();
    const bool is_bitc = method.calc_method()=="BITC" ? true : false;
    const bool uses_3body = (method.calc_method()=="TC" || method.calc_method()=="BITC") ? true : false;

    if (am_i_mpi_rank0) { *ost << " Start BAND calculation!" << std::endl; }

    bloch_states.set_filling(spin, kpoints, false, am_i_mpi_rank0, ost); // no need to mix the density matrix
    bloch_states.set_density(crystal_structure, kpoints, plane_wave_basis, mixes_density_matrix_,
                             mixing_beta_, true, is_bitc, am_i_mpi_rank0, ost);
    if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "preparation for BAND calc. (set filling and density)"); }

    // "iter" loop is required even for BAND calculation
    // because divergence correction depends on "phik_band".
    // (and also because of the accuracy of diagonalization while this reason is not relevant 
    //  when the loop numbers in block_davidson() is increased.)
    std::vector<std::vector<std::vector<Complex> > > eigenvalues_band_prev;
    bool is_energy_converged = false;
    for (int iter=0; iter<max_num_iterations_; iter++) 
    {
        if (am_i_mpi_rank0) { bloch_states.print_band_energies(kpoints, ost); }
        if (iter!=0) 
        {
            int icount = 0;
            double mean_absolute_error = 0.0;
            for (int ispin=0; ispin<spin.num_independent_spins(); ispin++)
            {
                for (int ik=0; ik<kpoints.num_irreducible_kpoints_band(); ik++)
                {
                    for (int iband=0; iband<bloch_states.num_bands_band()[ispin]; iband++)
                    {
                        mean_absolute_error += std::abs(bloch_states.eigenvalues_band()[ispin][ik][iband]
                                                        - eigenvalues_band_prev[ispin][ik][iband]);
                        icount++;
                    } // iband
                } // ik
            } // ispin
            mean_absolute_error /= static_cast<double>(icount);
            if (am_i_mpi_rank0) { *ost << "   mean absolute energy difference from the previous loop = " << mean_absolute_error << " Ht." << std::endl; }
            is_energy_converged = (mean_absolute_error < energy_tolerance_);
        } // iter!=0
        if (is_energy_converged)
        {
            if (am_i_mpi_rank0) { *ost << "  convergence is achieved!" << std::endl; }
            return;
        }
        eigenvalues_band_prev = bloch_states.eigenvalues_band(); // save old eigenvalues

        if (am_i_mpi_rank0) { *ost << " Iteration " << iter+1 << std::endl; }

        // set Hamoiltonian and diagonalize it
        block_davidson(parallelization, my_clock, method,
                       crystal_structure, symmetry,
                       potentials, spin, kpoints,
                       plane_wave_basis, bloch_states, total_energy, ost);

        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "diagonalization"); }
        if (am_i_mpi_rank0) { io_tc_files::dump_eigen(file_names, method, "BAND", bloch_states, ost); } // output for each iteration
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "dumping wave functions etc."); }
        if (am_i_mpi_rank0) { io_tc_files::dump_bandplot(file_names, crystal_structure,
                                                         kpoints, bloch_states, ost); } // dump tc_bandplot.dat
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "dumping tc_bandplot.dat"); }
        if (am_i_mpi_rank0) { *ost << std::endl; }
    } // iter

    // When the convergence is NOT achieved
    if (am_i_mpi_rank0) { bloch_states.print_band_energies(kpoints, ost); }
    if (am_i_mpi_rank0) { *ost << " CAUTION!! convergence is not achieved..." << std::endl; }
}
