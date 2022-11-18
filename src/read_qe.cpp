// [namspace read_qe]
// read Quantum-Espresso (QE) files
// declared in include/read_qe.hpp,
// defined in  read_qe.cpp, read_qe_xml.cpp, read_qe_wfc.f90

#include "include/header.hpp"

// Read QE files
void read_qe::read_qe(FileNames &file_names,
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
        read_qe_xml(file_names, method, crystal_structure, symmetry, potentials,
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
        read_qe_wfcs(file_names.qe_wfc(), method.calc_mode(),
                     symmetry, spin, kpoints, plane_wave_basis, bloch_states, ost);
    }
    bcast_qe(file_names, method, crystal_structure, symmetry, potentials, 
             spin, kpoints, plane_wave_basis, bloch_states, total_energy,
             am_i_mpi_rank0, ost);
}
 
// call read_qe_wfc() many times
void read_qe::read_qe_wfcs(const std::vector<std::vector<std::string> > &qe_wfc,
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
    
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<nks; ik++)
        {
            // QE names
            std::vector<int> mill(3*npw[ik]); // mill(3,igwx) in fortran. G-vector at each k-point
            std::vector<Complex> evc(npol*npw[ik]*nbnd[ispin]); // evc(npol*igwx, nbnd[is]) in fortran. Wave function at each k-point
                
            *ost << " Read QE wave functions from " << qe_wfc[ispin][ik] << std::endl;
            read_qe_wfc(qe_wfc[ispin][ik].length(), qe_wfc[ispin][ik].c_str(),
                        ik+1, npw[ik], nbnd[ispin], npol, mill.data(), evc.data());
                
            // set TC++ variables (Gvector_at_k & phik) from QE variables (mill & evc)
            plane_wave_basis.set_Gvector_at_k(symmetry, kpoints, ispin, ik, mill, calc_mode);
            bloch_states.set_phik_qe(ispin, ik, evc, calc_mode);
        } // ik
    } // ispin
    *ost << std::endl;
}

void read_qe::bcast_qe(const FileNames &file_names,
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
    if (method.calc_method()=="TC" || method.calc_method()=="BITC") { potentials.jastrow.bcast(); }
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
    
    if (am_i_mpi_rank0) { *ost << std::endl; }
}
