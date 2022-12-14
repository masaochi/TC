#include "include/header.hpp"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    Parallelization parallelization;
    const bool am_i_mpi_rank0 = parallelization.am_i_mpi_rank0();

    FileNames file_names;
    std::ostream* ost = am_i_mpi_rank0 ?
        new std::ofstream(file_names.tc_output()) : &std::cout;
    if (!(ost)) { error_messages::cannot_open(file_names.tc_output()); }
    *ost << std::setprecision(10);
    if (am_i_mpi_rank0) { std::cout << std::setprecision(10); }

    if (am_i_mpi_rank0)
    {
        *ost << "TC++ ver.1.1" << std::endl;
        parallelization.print(ost);
        *ost << std::endl;
    }

    MyClock my_clock;
    Method method;
    CrystalStructure crystal_structure;
    Symmetry symmetry;
    Potentials potentials;
    Spin spin;
    Kpoints kpoints;
    PlaneWaveBasis plane_wave_basis;
    BlochStates bloch_states;
    TotalEnergy total_energy;
    Diagonalization diagonalization;

    io_tc_files::read_input_in(file_names, diagonalization,
                               method, potentials,
                               kpoints, bloch_states,
                               am_i_mpi_rank0, ost);
    read_qe::read_qe(file_names, method, crystal_structure, symmetry,
                     potentials, spin, kpoints, plane_wave_basis,
                     bloch_states, total_energy, am_i_mpi_rank0, ost);
    if (!potentials.is_heg()) // upf files are not read for the homogeneous-electron-gas mode
    {
        read_upf::read_upf(file_names, crystal_structure, potentials,
                           plane_wave_basis, am_i_mpi_rank0, ost);
    }
    if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "initialization"); }

    if (spin.is_spinor()) { error_messages::stop("Non-collinear calculation is not supported..."); }

    if (method.calc_mode()=="SCF" && am_i_mpi_rank0) // dump SCF information for the subsequent BAND calc.
    {
        io_tc_files::dump_scfinfo(file_names, kpoints, plane_wave_basis, 
                                  bloch_states, ost);
        my_clock.print_time_from_save(ost, "dumping tc_scfinfo.dat");
    }
    else if (method.calc_mode()=="BAND") // read SCF information in BAND calc.
    {
        io_tc_files::set_scfinfo(file_names, method, spin, kpoints, 
                                 plane_wave_basis, bloch_states, am_i_mpi_rank0, ost);
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "reading tc_scfinfo.dat"); }
    }

    if ((method.calc_mode()=="SCF" && diagonalization.restarts()) // restart SCF calc.
        || method.calc_mode()=="BAND") // or BAND calc.
    {
        io_tc_files::read_eigen(file_names, method, "SCF", bloch_states, am_i_mpi_rank0, ost);
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "reading SCF eigenvalues and eigenvectors calculated by TC++"); }
    }

    if (method.calc_mode()=="BAND" && diagonalization.restarts()) // restart BAND calc.
    {
        io_tc_files::read_eigen(file_names, method, "BAND", bloch_states, am_i_mpi_rank0, ost);
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "reading BAND eigenvalues and eigenvectors calculated by TC++"); }
    }

    potentials.set_sum_of_Vaux(parallelization, method,
                               crystal_structure, kpoints);
    if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "setting potentials"); }
    if (am_i_mpi_rank0) { *ost << std::endl; }

    if (method.calc_mode()=="SCF")
    {
        diagonalization.scf(parallelization,
                            file_names, my_clock, method,
                            crystal_structure, symmetry,
                            potentials, spin, kpoints,
                            plane_wave_basis, 
                            bloch_states, total_energy, ost);
    }
    else if (method.calc_mode()=="BAND")
    {
        diagonalization.band(parallelization,
                             file_names, my_clock, method,
                             crystal_structure, symmetry,
                             potentials, spin, kpoints,
                             plane_wave_basis, 
                             bloch_states, total_energy, ost);
    }

    if (am_i_mpi_rank0) { my_clock.print_total_time(ost); }
    MPI_Finalize();
    return 0;
}
