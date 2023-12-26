// [namespace io_qmc_files]
// Read parameters.casl, write pwfn.data, write jastrow.plt

#ifndef TC_IO_QMC_FILES_HPP
#define TC_IO_QMC_FILES_HPP

namespace io_qmc_files
{

// dump jastrow.plt for plotting a Jastrow function using gnuplot
void dump_jastrow_plt(const FileNames &file_names,
                      const Jastrow &jastrow,
                      const bool dump_down_down,
                      std::ostream *ost);

// dump pwfn.data in SCF calculation (for CASINO)    
void dump_pwfn(const FileNames &file_names,
               const CrystalStructure &crystal_structure,
               const Kpoints &kpoints,
               const PlaneWaveBasis &plane_wave_basis,
               const BlochStates &bloch_states, 
               const TotalEnergy &total_energy,
               std::ostream *ost);

// read parameters.casl (Jastrow parameters) for CASINO
void read_casl(const FileNames &file_names,
               const Spin &spin,
               const CrystalStructure &crystal_structure,
               Jastrow &jastrow,
               const PlaneWaveBasis &plane_wave_basis,
               const bool am_i_mpi_rank0,
               std::ostream *ost);
void read_each_term_casl(std::ifstream &ifs, const int num_independent_spins, Jastrow &jastrow,
                         bool &is_RPA_read, bool &is_POLY_read, bool &cusp_poly);
bool read_RPA_casl(std::ifstream &ifs, const int num_independent_spins, Jastrow &jastrow);    
bool read_POLY_casl(std::ifstream &ifs, const int num_independent_spins, Jastrow &jastrow, bool &cusp_poly);

// write parameters.casl (Jastrow parameters) for CASINO
// when parameters.casl is not given and a default Jastrow is used.
void dump_casl(const FileNames &file_names,
               const Spin &spin,
               const Jastrow &jastrow,
               std::ostream *ost);    

} // namespace io_qmc_files

#endif // TC_IO_QMC_FILES_HPP
