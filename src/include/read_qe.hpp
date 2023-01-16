// [namspace read_qe]
// read Quantum-Espresso (QE) files
// declared in include/read_qe.hpp,
// defined in  read_qe.cpp, read_qe_xml.cpp, read_qe_wfc.f90

#ifndef TC_READ_QE_HPP
#define TC_READ_QE_HPP

namespace read_qe
{
    
// read QE files placed in the "file_names.qe_save_dir".    
void read_qe(FileNames &file_mames,
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
             std::ostream *ost);

void bcast_qe(const FileNames &file_mames,
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
              std::ostream *ost);

// call read_qe_wfc() many times in read_qe_wfcs()
void read_qe_wfcs(const std::vector<std::vector<std::string> > &qe_wfc,
                  const bool reads_binary,
                  const std::string &calc_mode,
                  const Symmetry &symmetry,
                  const Spin &spin,
                  const Kpoints &kpoints,
                  PlaneWaveBasis &plane_wave_basis, 
                  BlochStates &bloch_states,
                  std::ostream *ost);

extern "C"
{
    // read the QE wfc.dat files
    void read_qe_wfc(const int &nlength,
                     const char *qe_wfc_char,
                     const int &ik,
                     const int &igwx,
                     const int &nbnd,
                     const int &npol,
                     int *mill,
                     Complex *evc);

    // read the QE wfc.dat files (transformed to a non-binary format)
    void read_qe_wfc_nonbin(const int &nlength,
                            const char *qe_wfc_char,
                            const int &ik,
                            const int &igwx,
                            const int &nbnd,
                            const int &npol,
                            int *mill,
                            Complex *evc);
}

// read the QE xml file
void read_qe_xml(FileNames &file_names,
                 const Method &method,
                 CrystalStructure &crystal_structure,
                 Symmetry &symmetry,
                 Potentials &potentials,
                 Spin &spin, Kpoints &kpoints,
                 PlaneWaveBasis &plane_wave_basis,
                 BlochStates &bloch_states,
                 TotalEnergy &total_energy,
                 std::ostream *ost);

} // namespace read_qe

#endif // READ_QE_HPP
