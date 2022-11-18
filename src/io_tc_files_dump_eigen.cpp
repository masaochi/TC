// [namespace io_tc_files]
// Read input.in, read/write tc_wfc.dat (wave functions) and tc_energy.dat (orbital energies), etc.

// This file contains functions dumping "tc_wfc.dat" etc..

#include "include/header.hpp"

// dump eigenvalues & eigenvectors in each SCF loop
void io_tc_files::dump_eigen(const FileNames &file_names,
                             const Method &method,
                             const std::string &calc_mode,
                             const BlochStates &bloch_states,
                             std::ostream *ost)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    dump_wfc(file_names, method, calc_mode, bloch_states, ost);
    dump_energy(file_names, method, calc_mode, bloch_states, ost);
}

void io_tc_files::dump_wfc(const FileNames &file_names,
                           const Method &method,
                           const std::string &calc_mode,
                           const BlochStates &bloch_states,
                           std::ostream *ost)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    const std::string &file_name_ref = calc_mode=="SCF" ? 
        file_names.tc_wfc_scf() : file_names.tc_wfc_band();
    const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_ref = calc_mode=="SCF" ?
        bloch_states.phik_scf() : bloch_states.phik_band();    
    const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_left_ref = calc_mode=="SCF" ?
        bloch_states.phik_left_scf() : bloch_states.phik_left_band(); 

    if (calc_mode=="SCF")
    {
        *ost << " Dump SCF wave functions (" << file_name_ref << ")" << std::endl;
    }
    else if (calc_mode=="BAND")
    {
        *ost << " Dump BAND wave functions (" << file_name_ref << ")" << std::endl;
    }
    std::ofstream ofs(file_name_ref, std::ios::out | std::ios::binary); // binary
    if (ofs.fail()) { error_messages::cannot_open(file_name_ref); }

    // dump a size of array
    const int num_independent_spins = phik_ref.size();
    ofs.write(reinterpret_cast<const char*>(&num_independent_spins), sizeof(int));
    const int num_irreducible_kpoints_scf = phik_ref[0].size();
    ofs.write(reinterpret_cast<const char*>(&num_irreducible_kpoints_scf), sizeof(int));
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        const int num_bands_scf = phik_ref[ispin][0].size();
        ofs.write(reinterpret_cast<const char*>(&num_bands_scf), sizeof(int));
    }
    const int num_spinor = phik_ref[0][0][0].size();
    ofs.write(reinterpret_cast<const char*>(&num_spinor), sizeof(int));
    for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
    {
        const int num_G_at_k = phik_ref[0][ik][0][0].size();
        ofs.write(reinterpret_cast<const char*>(&num_G_at_k), sizeof(int)); 
    }

    // dump phik
    for (int ispin=0; ispin<phik_ref.size(); ispin++)
    {
        for (int ik=0; ik<phik_ref[ispin].size(); ik++)
        {
            for (int iband=0; iband<phik_ref[ispin][ik].size(); iband++)
            {
                for (int ispinor=0; ispinor<phik_ref[ispin][ik][iband].size(); ispinor++)
                {
                    for (int ipw_at_k=0; ipw_at_k<phik_ref[ispin][ik][iband][ispinor].size(); ipw_at_k++)
                    {
                        // not a clever way...
                        double re = phik_ref[ispin][ik][iband][ispinor](ipw_at_k).real();
                        double im = phik_ref[ispin][ik][iband][ispinor](ipw_at_k).imag();
                        ofs.write(reinterpret_cast<const char*>(&re), sizeof(double));
                        ofs.write(reinterpret_cast<const char*>(&im), sizeof(double));
                    } // ipw_at_k
                } // ispinor
            } // iband
        } // ik
    } // ispin

    if (method.calc_method()=="BITC")
    {
        // dump phik_left
        for (int ispin=0; ispin<phik_ref.size(); ispin++)
        {
            for (int ik=0; ik<phik_ref[ispin].size(); ik++)
            {
                for (int iband=0; iband<phik_ref[ispin][ik].size(); iband++)
                {
                    for (int ispinor=0; ispinor<phik_ref[ispin][ik][iband].size(); ispinor++)
                    {
                        for (int ipw_at_k=0; ipw_at_k<phik_ref[ispin][ik][iband][ispinor].size(); ipw_at_k++)
                        {
                            double re = phik_left_ref[ispin][ik][iband][ispinor](ipw_at_k).real();
                            double im = phik_left_ref[ispin][ik][iband][ispinor](ipw_at_k).imag();
                            ofs.write(reinterpret_cast<const char*>(&re), sizeof(double));
                            ofs.write(reinterpret_cast<const char*>(&im), sizeof(double));
                        } // ipw_at_k
                    } // ispinor
                } // iband
            } // ik
        } // ispin
    } // if (BITC)

    ofs.close();
}

void io_tc_files::dump_energy(const FileNames &file_names,
                              const Method &method,
                              const std::string &calc_mode,
                              const BlochStates &bloch_states,
                              std::ostream *ost)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    const std::string &file_name_ref = calc_mode=="SCF" ? 
        file_names.tc_energy_scf() : file_names.tc_energy_band();
    const std::vector<std::vector<std::vector<Complex> > > &eigenvalues_ref = calc_mode=="SCF" ?
        bloch_states.eigenvalues_scf() : bloch_states.eigenvalues_band();

    if (calc_mode=="SCF") 
    {
        *ost << " Dump SCF orbital energies (" << file_name_ref << ")" << std::endl;
    }
    else if (calc_mode=="BAND")
    {
        *ost << " Dump BAND orbital energies (" << file_name_ref << ")" << std::endl;
    }
    std::ofstream ofs(file_name_ref, std::ios::out | std::ios::binary); // binary
    if (ofs.fail()) { error_messages::cannot_open(file_name_ref); }

    // dump a size of array
    const int num_independent_spins = eigenvalues_ref.size();
    ofs.write(reinterpret_cast<const char*>(&num_independent_spins), sizeof(int));
    const int num_irreducible_kpoints_scf = eigenvalues_ref[0].size();
    ofs.write(reinterpret_cast<const char*>(&num_irreducible_kpoints_scf), sizeof(int));
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        const int num_bands_scf = eigenvalues_ref[ispin][0].size();
        ofs.write(reinterpret_cast<const char*>(&num_bands_scf), sizeof(int));
    }

    // dump eigenvalues_scf
    for (int ispin=0; ispin<eigenvalues_ref.size(); ispin++)
    {
        for (int ik=0; ik<eigenvalues_ref[ispin].size(); ik++)
        {
            for (int iband=0; iband<eigenvalues_ref[ispin][ik].size(); iband++)
            {
                double re = eigenvalues_ref[ispin][ik][iband].real();
                double im = eigenvalues_ref[ispin][ik][iband].imag();
                ofs.write(reinterpret_cast<const char*>(&re), sizeof(double));
                ofs.write(reinterpret_cast<const char*>(&im), sizeof(double));
            } // iband
        } // ik
    } // ispin

    ofs.close();
}
