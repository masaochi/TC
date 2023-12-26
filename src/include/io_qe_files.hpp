// [namspace io_qe_files]
// read Quantum-Espresso (QE) files
// declared in include/read_qe.hpp, include/io_qe_files_upf.hpp
// defined in  io_qe_files_read_qe.cpp, io_qe_files_read_qe_xml.cpp, io_qe_files_read_qe_wfc.f90, io_qe_files_read_upf.cpp

#ifndef TC_IO_QE_FILES_HPP
#define TC_IO_QE_FILES_HPP

namespace io_qe_files
{
    
// read QE files placed in the "file_names.qe_save_dir".    
void read(FileNames &file_mames,
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

void bcast(const FileNames &file_mames,
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

// call read_wfc() many times in read_wfcs()
void read_wfcs(const std::vector<std::vector<std::string> > &qe_wfc,
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
    void read_wfc(const int &nlength,
                  const char *qe_wfc_char,
                  const int &ik,
                  const int &igwx,
                  const int &nbnd,
                  const int &npol,
                  const bool &gamma_only,
                  int *mill,
                  Complex *evc);

    // read the QE wfc.dat files (transformed to a non-binary format)
    void read_wfc_nonbin(const int &nlength,
                         const char *qe_wfc_char,
                         const int &ik,
                         const int &igwx,
                         const int &nbnd,
                         const int &npol,
                         const bool &gamma_only,
                         int *mill,
                         Complex *evc);
}

// read the QE xml file
void read_xml(FileNames &file_names,
              const Method &method,
              CrystalStructure &crystal_structure,
              Symmetry &symmetry,
              Potentials &potentials,
              Spin &spin, Kpoints &kpoints,
              PlaneWaveBasis &plane_wave_basis,
              BlochStates &bloch_states,
              TotalEnergy &total_energy,
              std::ostream *ost);

} // namespace io_qe_files

#endif // IO_QE_FILES_HPP
