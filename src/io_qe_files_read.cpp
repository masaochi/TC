// [namspace io_qe_files]
// read Quantum-Espresso (QE) files

#include "include/header.hpp"

// Read QE files
void io_qe_files::read(FileNames &file_names,
                       const Method &method,
                       CrystalStructure &crystal_structure,
                       Symmetry &symmetry,
                       Potentials &potentials,
                       Spin &spin,
                       Kpoints &kpoints, 
                       PlaneWaveBasis &plane_wave_basis, 
                       BlochStates &bloch_states,
                       TotalEnergy &total_energy,
                       const bool am_i_mpi_rank0,
                       std::ostream *ost)
{
    if (am_i_mpi_rank0) // open files
    {
        // read the xml file
        file_names.set_qe_xml_file_name();
        read_xml(file_names, method, crystal_structure, symmetry, potentials,
                 spin, kpoints, plane_wave_basis, bloch_states, total_energy, ost);
        
        // read the wfc files
        if (method.calc_mode()=="SCF")
        {
            file_names.set_qe_wfc_file_names(spin, kpoints.num_irreducible_kpoints_scf());
        }
        else if (method.calc_mode()=="BAND")
        {
            file_names.set_qe_wfc_file_names(spin, kpoints.num_irreducible_kpoints_band());
        }
        read_wfcs(file_names.qe_wfc(), file_names.reads_binary(), method.calc_mode(),
                  symmetry, spin, kpoints, plane_wave_basis, bloch_states, ost);
    }
    bcast(file_names, method, crystal_structure, symmetry, potentials, 
          spin, kpoints, plane_wave_basis, bloch_states, total_energy,
          am_i_mpi_rank0, ost);
}
 
// call read_qe_wfc() many times
void io_qe_files::read_wfcs(const std::vector<std::vector<std::string> > &qe_wfc,
                            const bool reads_binary,
                            const std::string &calc_mode,
                            const Symmetry &symmetry,
                            const Spin &spin,
                            const Kpoints &kpoints,
                            PlaneWaveBasis &plane_wave_basis, 
                            BlochStates &bloch_states,
                            std::ostream *ost)
{
    const int num_independent_spins = spin.num_independent_spins();

    // set an array size of phik (evc in QE)
    bloch_states.resize_phik_qe(plane_wave_basis, spin.is_spinor(), num_independent_spins, calc_mode);

    // Left-hand side = QE name, Right-hand side = TC++ name
    const std::vector<int> nbnd = bloch_states.num_bands_qe();   // (nspin==2) nbnd[0]=nbnd_up, nbnd[1]=nbnd_dw, (nspin==1, 4) nbnd[0]=nbnd
    const std::vector<int> npw = calc_mode=="SCF" ?
        plane_wave_basis.num_G_at_k_scf() : plane_wave_basis.num_G_at_k_band();  // npw[nks]: num. of plane waves at each k-point
    const int nks = calc_mode=="SCF" ?
        kpoints.num_irreducible_kpoints_scf() : kpoints.num_irreducible_kpoints_band(); // num. of irreducible k-points
    const int npol = spin.is_spinor() ? 2 : 1; // (non-collinear) 2 (otherwise) 1
    const bool gamma_only = plane_wave_basis.gamma_only();
    
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<nks; ik++)
        {
            // QE names
            std::vector<int> mill(3*npw[ik]); // mill(3,igwx) in fortran. G-vector at each k-point
            std::vector<Complex> evc(npol*npw[ik]*nbnd[ispin]); // evc(npol*igwx, nbnd[is]) in fortran. Wave function at each k-point
                
            *ost << " Read QE wave functions from " << qe_wfc[ispin][ik] << std::endl;
            if (reads_binary)
            {
                read_wfc(qe_wfc[ispin][ik].length(), qe_wfc[ispin][ik].c_str(),
                         ik+1, npw[ik], nbnd[ispin], npol, gamma_only, mill.data(), evc.data());
            }
            else // reads a non-binary format (for test calculation)
            {
                read_wfc_nonbin(qe_wfc[ispin][ik].length(), qe_wfc[ispin][ik].c_str(),
                                ik+1, npw[ik], nbnd[ispin], npol, gamma_only, mill.data(), evc.data());
            }

            // set TC++ variables (Gvector_at_k & phik) from QE variables (mill & evc)
            int zero_index; // used only if gamma_only==true
            plane_wave_basis.set_Gvector_at_k(symmetry, kpoints, ispin, ik, mill, zero_index, calc_mode); // If gamma_only, change num_G_at_k. (u(-G)=u(G)^*)
            bloch_states.set_phik_qe(ispin, ik, evc, calc_mode, gamma_only, zero_index); // If gamma_only, apply u(-G)=u(G)^*
        } // ik
    } // ispin
    *ost << std::endl;
}

void io_qe_files::bcast(const FileNames &file_names,
                        const Method &method,
                        CrystalStructure &crystal_structure,
                        Symmetry &symmetry,
                        Potentials &potentials,
                        Spin &spin,
                        Kpoints &kpoints, 
                        PlaneWaveBasis &plane_wave_basis, 
                        BlochStates &bloch_states,
                        TotalEnergy &total_energy,
                        const bool am_i_mpi_rank0,
                        std::ostream *ost)
{
    if (am_i_mpi_rank0) { *ost << " Bcast QE input variables" << std::endl; }

    crystal_structure.bcast(am_i_mpi_rank0);
    symmetry.bcast(am_i_mpi_rank0);
    spin.bcast(am_i_mpi_rank0);
    kpoints.bcast_qe_input(method.calc_mode(), am_i_mpi_rank0);
    plane_wave_basis.bcast(method.calc_mode(), am_i_mpi_rank0, true);

    // NOTE: spin.bcast() should be done before bloch_states.bcast_qe_input()
    bloch_states.bcast_qe_input(plane_wave_basis, spin.is_spinor(),
                                method.calc_mode(), method.calc_method(),
                                am_i_mpi_rank0);

    const int num_irreducible_kpoints = method.calc_mode()=="SCF" ?
        kpoints.num_irreducible_kpoints_scf() : kpoints.num_irreducible_kpoints_band();

    // NOTE: spin.bcast() & kpoints.bcast_qe_input() & bloch_states.bcast_qe_input()
    //       should be done before total_energy.bcast_qe_input()
    total_energy.bcast_qe_input(spin.num_independent_spins(),
                                num_irreducible_kpoints,
                                bloch_states.num_bands_scf(),
                                method.calc_mode(), method.calc_method(),
                                potentials.is_heg());

    // Jastrow parameters unnormalization using the cell volume and the num. of electrons
    if (method.calc_method()=="TC" || method.calc_method()=="BITC")
    {
        potentials.jastrow.unnormalize_A_long(crystal_structure.unit_cell_volume(),
                                              bloch_states.num_electrons(), spin);
        potentials.jastrow.impose_cusp(true); // determine F_long
    }

    if (am_i_mpi_rank0) { *ost << std::endl; }
}
