// [namespace io_tc_files]
// Read input.in, read/write tc_wfc.dat (wave functions) and tc_energy.dat (orbital energies), etc.

// This file contains functions dumping "tc_wfc.dat" etc..
// (For functions regarding "input.in", see io_tc_files_input_in.cpp)

#include "include/header.hpp"

// read eigenvalues & eigenvectors for restarting SCF calculation
void io_tc_files::read_eigen(const FileNames &file_names,
                             const Method &method,
                             const std::string &calc_mode,
                             BlochStates &bloch_states,
                             const bool am_i_mpi_rank0,
                             std::ostream *ost)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    const bool is_bitc = method.calc_method()=="BITC" ? true : false;

    if (am_i_mpi_rank0) { read_wfc(file_names, method, calc_mode, bloch_states, ost); }
    bloch_states.bcast_phik(is_bitc, calc_mode);
    if (am_i_mpi_rank0) { read_energy(file_names, method, calc_mode, bloch_states, ost); }
    bloch_states.bcast_eigenvalues(calc_mode);
}

// non-collinear calculation not supported
void io_tc_files::read_wfc(const FileNames &file_names,
                           const Method &method,
                           const std::string &calc_mode,
                           BlochStates &bloch_states,
                           std::ostream *ost)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    const std::string &file_name_ref = calc_mode=="SCF" ?
        file_names.tc_wfc_scf() : file_names.tc_wfc_band();
    const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_ref = calc_mode=="SCF" ?
        bloch_states.phik_scf() : bloch_states.phik_band();

    if (calc_mode=="SCF")
    {
        *ost << " Read SCF wave functions (" << file_name_ref << ")" << std::endl;
    }
    else if (calc_mode=="BAND")
    {
        *ost << " Read BAND wave functions (" << file_name_ref << ")" << std::endl;
    }
    std::ifstream ifs(file_name_ref, std::ios::in | std::ios::binary); // binary
    if (ifs.fail()) { error_messages::cannot_open(file_name_ref); }

    // Read a size of array
    int num_independent_spins, num_irreducible_kpoints_scf, num_spinor;
    ifs.read(reinterpret_cast<char*>(&num_independent_spins), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&num_irreducible_kpoints_scf), sizeof(int));

    std::vector<int> num_bands_scf(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        ifs.read(reinterpret_cast<char*>(&num_bands_scf[ispin]), sizeof(int));
    }
    ifs.read(reinterpret_cast<char*>(&num_spinor), sizeof(int));

    std::vector<int> num_G_at_k(num_irreducible_kpoints_scf);
    for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
    {
        ifs.read(reinterpret_cast<char*>(&num_G_at_k[ik]), sizeof(int)); 
    }

    // check consistency
    if (num_independent_spins != phik_ref.size()) { error_messages::stop("num_independent_spins is inconsistent (read_wfc)"); }
    if (num_irreducible_kpoints_scf != phik_ref[0].size()) { error_messages::stop("num_irreducible_kpoints_scf is inconsistent (read_wfc)"); }
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        if (num_bands_scf[ispin] != phik_ref[ispin][0].size()) { error_messages::stop("num_bands_tc is inconsistent (read_wfc)"); }
    }
    if (num_spinor != phik_ref[0][0][0].size()) { error_messages::stop("num_spinor is inconsistent (read_wfc)"); }
    for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
    {
        if (num_G_at_k[ik] != phik_ref[0][ik][0][0].size()) { error_messages::stop("num_G_at_k is inconsistent (read_wfc)"); }
    }

    // non-collinear calculation not supported
    assert(num_spinor==1);

    // read phik
    Eigen::VectorXcd wfc_temp;
    for (int ispin=0; ispin<phik_ref.size(); ispin++)
    {
        for (int ik=0; ik<phik_ref[ispin].size(); ik++)
        {
            wfc_temp.resize(num_G_at_k[ik]);
            for (int iband=0; iband<phik_ref[ispin][ik].size(); iband++)
            {
                for (int ispinor=0; ispinor<phik_ref[ispin][ik][iband].size(); ispinor++)
                {
                    for (int ipw_at_k=0; ipw_at_k<phik_ref[ispin][ik][iband][ispinor].size(); ipw_at_k++)
                    {
                        double re, im;
                        ifs.read(reinterpret_cast<char*>(&re), sizeof(double));
                        ifs.read(reinterpret_cast<char*>(&im), sizeof(double));
                        wfc_temp(ipw_at_k) = re + I*im;
                    } // ipw_at_k
                } // ispinor
                bloch_states.set_phik_from_tc_wfc(ispin, ik, iband, wfc_temp, "right", calc_mode);
            } // iband
        } // ik
    } // ispin

    if (method.calc_method()=="BITC")
    {
        // read phik_left
        for (int ispin=0; ispin<phik_ref.size(); ispin++)
        {
            for (int ik=0; ik<phik_ref[ispin].size(); ik++)
            {
                wfc_temp.resize(num_G_at_k[ik]);
                for (int iband=0; iband<phik_ref[ispin][ik].size(); iband++)
                {
                    for (int ispinor=0; ispinor<phik_ref[ispin][ik][iband].size(); ispinor++)
                    {
                        for (int ipw_at_k=0; ipw_at_k<phik_ref[ispin][ik][iband][ispinor].size(); ipw_at_k++)
                        {
                            double re, im;
                            ifs.read(reinterpret_cast<char*>(&re), sizeof(double));
                            ifs.read(reinterpret_cast<char*>(&im), sizeof(double));
                            wfc_temp(ipw_at_k) = re + I*im;
                        } // ipw_at_k
                    } // ispinor
                    bloch_states.set_phik_from_tc_wfc(ispin, ik, iband, wfc_temp, "left", calc_mode);
                } // iband
            } // ik
        } // ispin
        if (ifs.eof()) { error_messages::stop("BITC should be restarted after BITC. (Restarts of HF->BITC and TC->BITC are not supported.)"); }
    }
    else // HF or TC
    {
        double re;
        ifs.read(reinterpret_cast<char*>(&re), sizeof(double));
        if (!ifs.eof()) { error_messages::stop("HF,TC should not be restarted after BITC."); }
    } // if (BITC)

    ifs.close();
}

void io_tc_files::read_energy(const FileNames &file_names,
                              const Method &method,
                              const std::string &calc_mode,
                              BlochStates &bloch_states,
                              std::ostream *ost)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    const std::string &file_name_ref = calc_mode=="SCF" ?
        file_names.tc_energy_scf() : file_names.tc_energy_band();
    const std::vector<std::vector<std::vector<Complex> > > &eigenvalues_ref = calc_mode=="SCF" ?
        bloch_states.eigenvalues_scf() : bloch_states.eigenvalues_band();

    if (calc_mode=="SCF") 
    {
        *ost << " Read SCF orbital energies (" << file_name_ref << ")" << std::endl;
    }
    else if (calc_mode=="BAND")
    {
        *ost << " Read BAND orbital energies (" << file_name_ref << ")" << std::endl;
    }
    std::ifstream ifs(file_name_ref, std::ios::in | std::ios::binary); // binary
    if (ifs.fail()) { error_messages::cannot_open(file_name_ref); }

    // Read a size of array
    int num_independent_spins, num_irreducible_kpoints_scf, num_spinor;
    ifs.read(reinterpret_cast<char*>(&num_independent_spins), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&num_irreducible_kpoints_scf), sizeof(int));
    std::vector<int> num_bands_scf(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        ifs.read(reinterpret_cast<char*>(&num_bands_scf[ispin]), sizeof(int));
    }

    // check consistency
    if (num_independent_spins != eigenvalues_ref.size()) { error_messages::stop("num_independent_spins is inconsistent (read_energy)"); }
    if (num_irreducible_kpoints_scf != eigenvalues_ref[0].size()) { error_messages::stop("num_irreducible_kpoints_scf is inconsistent (read_energy)"); }
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        if (num_bands_scf[ispin] != eigenvalues_ref[ispin][0].size()) { error_messages::stop("num_bands_tc is inconsistent (read_energy)"); }
    }

    // dump eigenvalues_scf
    for (int ispin=0; ispin<eigenvalues_ref.size(); ispin++)
    {
        Eigen::VectorXcd eigenvalues_temp(num_bands_scf[ispin]);
        for (int ik=0; ik<eigenvalues_ref[ispin].size(); ik++)
        {
            for (int iband=0; iband<eigenvalues_ref[ispin][ik].size(); iband++)
            {
                double re, im;
                ifs.read(reinterpret_cast<char*>(&re), sizeof(double));
                ifs.read(reinterpret_cast<char*>(&im), sizeof(double));
                eigenvalues_temp[iband] = re + I*im;
            } // iband
            bloch_states.set_eigenvalues_from_tc_energy(ispin, ik, eigenvalues_temp, calc_mode);
        } // ik
    } // ispin

    ifs.close();
}
