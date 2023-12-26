// [namespace io_tc_files]
// Read input.in, read/write tc_wfc.dat (wave functions) and tc_energy.dat (orbital energies), etc.

#ifndef TC_IO_TC_FILES_HPP
#define TC_IO_TC_FILES_HPP

namespace io_tc_files
{

// input.in
void read_input_in(FileNames &file_names,
                   Diagonalization &diagonalization,
                   Method &method, Potentials &potentials,
                   Kpoints &kpoints, BlochStates &bloch_states, 
                   const bool am_i_mpi_rank0,
                   std::ostream *ost);
void bcast_input_in(FileNames &file_names,
                    Diagonalization &diagonalization,
                    Method &method, Potentials &potentials,
                    Kpoints &kpoints, BlochStates &bloch_states,
                    const bool am_i_mpi_rank0,
                    std::ostream *ost);

// read tc_wfc.dat & tc_energy.dat (eigenenergies and eigenfunctions) in each SCF loop
void read_eigen(const FileNames &file_names,
                const Method &method,
                const std::string &calc_mode,
                BlochStates &bloch_states,
                const bool am_i_mpi_rank0,
                std::ostream *ost);
void read_wfc(const FileNames &file_names,
              const Method &method,
              const std::string &calc_mode,
              BlochStates &bloch_states, 
              std::ostream *ost);
void read_energy(const FileNames &file_names,
                 const Method &method,
                 const std::string &calc_mode,
                 BlochStates &bloch_states, 
                 std::ostream *ost);

// dump tc_wfc.dat & tc_energy.dat (eigenenergies and eigenfunctions) in each SCF loop
void dump_eigen(const FileNames &file_names,
                const Method &method,
                const std::string &calc_mode,
                const BlochStates &bloch_states, 
                std::ostream *ost);    
void dump_wfc(const FileNames &file_names,
              const Method &method,
              const std::string &calc_mode,
              const BlochStates &bloch_states, 
              std::ostream *ost);
void dump_energy(const FileNames &file_names,
                 const Method &method,
                 const std::string &calc_mode,
                 const BlochStates &bloch_states, 
                 std::ostream *ost);

// read & dump tc_scfinfo.dat (SCF information) that is necessary in BAND calculation
void set_scfinfo(const FileNames &file_names,
                 const Method &method,
                 const Spin &spin,
                 Kpoints &kpoints,
                 PlaneWaveBasis &plane_wave_basis,
                 BlochStates &bloch_states,
                 const bool am_i_mpi_rank0, 
                 std::ostream *ost);
void read_scfinfo(const FileNames &file_names,
                  const Spin &spin,
                  Kpoints &kpoints,
                  PlaneWaveBasis &plane_wave_basis,
                  BlochStates &bloch_states,
                  std::ostream *ost);    
void bcast_scfinfo(Kpoints &kpoints,
                   PlaneWaveBasis &plane_wave_basis,
                   BlochStates &bloch_states,
                   const bool is_spinor,
                   const std::string &calc_method,
                   const bool am_i_mpi_rank0);
void dump_scfinfo(const FileNames &file_names,
                  const Kpoints &kpoints,
                  const PlaneWaveBasis &plane_wave_basis,
                  const BlochStates &bloch_states,
                  std::ostream *ost);

// dump tc_bandplot.dat (BAND eigenvalues)
void dump_bandplot(const FileNames &file_names,
                   const CrystalStructure &crystal_structure,
                   const Kpoints &kpoints,
                   const BlochStates &bloch_states,
                   std::ostream *ost);
void dump_bandplot_each(const std::string &tc_bandplot,
                        const int &ispin,
                        const CrystalStructure &crystal_structure,
                        const Kpoints &kpoints,
                        const BlochStates &bloch_states,
                        std::ostream *ost);

// read & dump crystal structure (mainly for restarting structural opt., but can always be used)
void read_crystal_structure(const FileNames &file_names,
                            CrystalStructure &crystal_structure,
                            const bool am_i_mpi_rank0,
                            std::ostream *ost);
void dump_crystal_structure(const FileNames &file_names,
                            const CrystalStructure &crystal_structure,
                            std::ostream *ost);

} // namespace io_tc_files

#endif // TC_IO_TC_FILES_HPP
