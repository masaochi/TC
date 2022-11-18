// [namespace io_tc_files]
// Read input.in, read/write tc_wfc.dat (wave functions) and tc_energy.dat (orbital energies), etc.

#include "include/header.hpp"

// dump tc_bandplot.dat (BAND eigenvalues)
void io_tc_files::dump_bandplot(const FileNames &file_names,
                                const CrystalStructure &crystal_structure,
                                const Kpoints &kpoints,
                                const BlochStates &bloch_states,
                                std::ostream *ost)
{
    const int num_independent_spins = bloch_states.eigenvalues_band().size();
    if (num_independent_spins==1)
    {
        dump_bandplot_each(file_names.tc_bandplot(), 0,
                           crystal_structure, kpoints,
                           bloch_states, ost);
    }
    else
    {
        dump_bandplot_each(file_names.tc_bandplot_up(), 0,
                           crystal_structure, kpoints,
                           bloch_states, ost);
        dump_bandplot_each(file_names.tc_bandplot_dn(), 1,
                           crystal_structure, kpoints,
                           bloch_states, ost);
    }
}

void io_tc_files::dump_bandplot_each(const std::string &tc_bandplot,
                                     const int &ispin,
                                     const CrystalStructure &crystal_structure,
                                     const Kpoints &kpoints,
                                     const BlochStates &bloch_states,
                                     std::ostream *ost)
{
    *ost << " Dump BAND eigenvalues (" << tc_bandplot << ")" << std::endl;
    std::ofstream ofs(tc_bandplot, std::ios::out);
    if (ofs.fail()) { error_messages::cannot_open(tc_bandplot); }
    ofs << std::setprecision(10);

    ofs << "# ka, kb, kc in crystal coordinate, |k| (/bohr), real & imaginary parts of energy (eV)" << std::endl;
    ofs << "# Fermi energy = " << bloch_states.fermi_energy() * Ht_in_eV << " eV" << std::endl;
    for (int iband=0; iband<bloch_states.eigenvalues_band()[ispin][0].size(); iband++)
    {
        double norm_k = 0.0;
        for (int ik=0; ik<bloch_states.eigenvalues_band()[ispin].size(); ik++)
        {
            // k-point information
            for (int idim=0; idim<3; idim++)
            {
                ofs << kpoints.kvectors_band()[ik][0](idim) << " ";
            }
            if (ik!=0)
            {
                Eigen::Vector3d kdiff = 
                    kpoints.kvectors_band()[ik][0] - kpoints.kvectors_band()[ik-1][0];
                // norm in the cartesian coordinate
                norm_k += (crystal_structure.reciprocal_vectors().transpose() * kdiff).norm();
            }
            ofs << norm_k << " ";
            
            // BAND eigenvalue
            Complex ene_in_ev = bloch_states.eigenvalues_band()[ispin][ik][iband] * Ht_in_eV;
            ofs << ene_in_ev.real() << " " << ene_in_ev.imag() << std::endl;
        } // k-point
        ofs << std::endl;
    } // band

    ofs.close();
}
