// [class Diagonalization]
// scf_loop, convergence check, control paramters, etc.

#include "include/header.hpp"

void Diagonalization::set(const bool restarts,
                          const double &energy_tolerance, const double &charge_tolerance,
                          const int &max_num_iterations, const double &mixing_beta,
                          const int &num_refresh_david, const int &max_num_blocks_david,
                          const std::string &calc_mode, const bool is_heg)
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

    if (mixing_beta > -1e-8) 
    {
        mixing_beta_ = mixing_beta;
    }
    else
    {
        error_messages::inappropriate_argument("mixing_beta", mixing_beta,
                                               "should not be negative");
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
}

void Diagonalization::bcast(const bool am_i_mpi_rank0)
{
    MPI_Bcast(&restarts_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&energy_tolerance_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&charge_tolerance_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_num_iterations_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mixing_beta_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_refresh_david_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_num_blocks_david_, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void Diagonalization::scf(Parallelization &parallelization,
                          const FileNames &file_names, MyClock &my_clock,
                          const Method &method, const CrystalStructure &crystal_structure,
                          const Symmetry &symmetry, const Potentials &potentials,
                          const Spin &spin, const Kpoints &kpoints,
                          PlaneWaveBasis &plane_wave_basis, 
                          BlochStates &bloch_states, TotalEnergy &total_energy,
                          std::ostream *ost)
{
    const bool am_i_mpi_rank0 = parallelization.am_i_mpi_rank0();
    const bool is_bitc = method.calc_method()=="BITC" ? true : false;
    const bool uses_3body = (method.calc_method()=="TC" || method.calc_method()=="BITC") ? true : false;

    if (am_i_mpi_rank0) { *ost << " Start SCF calculation!" << std::endl; }

    bool is_energy_converged, is_energy_converged_prev, is_charge_converged, is_charge_converged_prev;
    is_energy_converged = is_energy_converged_prev = is_charge_converged = is_charge_converged_prev = false;
    for (int iter=0; iter<max_num_iterations_; iter++)
    {
        bool is_first_iter = iter==0 ? true : false;

        // Preparation for SCF calc.
        bloch_states.set_filling(spin, kpoints, am_i_mpi_rank0, ost);
        if (!is_first_iter) 
        {
            total_energy.calc_total_energy(spin, kpoints,
                                           bloch_states.filling(), bloch_states.eigenvalues_scf(),
                                           bloch_states.fermi_energy(), uses_3body, am_i_mpi_rank0, ost);
            is_energy_converged =
                (std::abs(total_energy.total_energy_difference()) < energy_tolerance_);
        }
        bloch_states.set_density(crystal_structure, kpoints, plane_wave_basis, mixing_beta_, 
                                 is_first_iter, is_bitc, am_i_mpi_rank0, ost);
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "preparation for SCF calc. (set filling and density)"); }

        // Convergence check
        if (!is_first_iter) { is_charge_converged = (bloch_states.density_difference() < charge_tolerance_); }

        if (is_energy_converged && is_energy_converged_prev &&
            is_charge_converged && is_charge_converged_prev)
        {
            if (am_i_mpi_rank0) { *ost << "  convergence is achieved!" << std::endl; }
            return;
        }
        is_energy_converged_prev = is_energy_converged;
        is_charge_converged_prev = is_charge_converged;
        
        if (am_i_mpi_rank0) { *ost << " Iteration " << iter+1 << std::endl; }

        // set Hamoiltonian and diagonalize it
        block_davidson(parallelization, my_clock, method,
                       crystal_structure, symmetry,
                       potentials, spin, kpoints,
                       plane_wave_basis, 
                       bloch_states, total_energy, ost);

        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "diagonalization"); }
        if (am_i_mpi_rank0) { io_tc_files::dump_eigen(file_names, method, "SCF", bloch_states, ost); } // output for each iteration
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "dumping wave functions etc."); }
        if (am_i_mpi_rank0) { *ost << std::endl; }
    } // iter

    // When the convergence is NOT achieved
    bloch_states.set_filling(spin, kpoints, am_i_mpi_rank0, ost);
    total_energy.calc_total_energy(spin, kpoints,
                                   bloch_states.filling(), bloch_states.eigenvalues_scf(),
                                   bloch_states.fermi_energy(), uses_3body, am_i_mpi_rank0, ost);
    
    if (am_i_mpi_rank0) { *ost << " CAUTION!! convergence is not achieved..." << std::endl; }
}

void Diagonalization::band(Parallelization &parallelization,
                           const FileNames &file_names, MyClock &my_clock,
                           const Method &method, const CrystalStructure &crystal_structure,
                           const Symmetry &symmetry, const Potentials &potentials,
                           const Spin &spin, const Kpoints &kpoints,
                           PlaneWaveBasis &plane_wave_basis, 
                           BlochStates &bloch_states, TotalEnergy &total_energy,
                           std::ostream *ost)
{
    const bool am_i_mpi_rank0 = parallelization.am_i_mpi_rank0();
    const bool is_bitc = method.calc_method()=="BITC" ? true : false;
    const bool uses_3body = (method.calc_method()=="TC" || method.calc_method()=="BITC") ? true : false;

    if (am_i_mpi_rank0) { *ost << " Start BAND calculation!" << std::endl; }

    bloch_states.set_filling(spin, kpoints, am_i_mpi_rank0, ost);
    bloch_states.set_density(crystal_structure, kpoints, plane_wave_basis, mixing_beta_, 
                             true, is_bitc, am_i_mpi_rank0, ost);
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
