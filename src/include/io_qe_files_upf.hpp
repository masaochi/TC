// [namespace io_qe_files::upf]
// Read upf (pseudopotential) files and MPI_Bcast

#ifndef TC_IO_QE_FILES_UPF_HPP
#define TC_IO_QE_FILES_UPF_HPP

namespace io_qe_files::upf
{

void read_upf(FileNames &file_names,
              const AtomicSpecies &atomic_species,
              CrystalStructure &crystal_structure,
              Potentials &potentials,
              PlaneWaveBasis &plane_wave_basis,
              const bool am_i_mpi_rank0,
              std::ostream *ost);

bool check_upf_ver2(std::ifstream &ifs, std::ostream *ost);
void read_upf_header(const AtomicSpecies &atomic_species,
                     std::ifstream &ifs, const std::string &file_name_upf,
                     const bool is_upf_ver2, int &atomic_number,
                     double &Z_valence, int &mesh_size, int &num_projectors,
                     std::ostream *ost);
void read_upf_mesh(std::ifstream &ifs, std::vector<double> &pseudo_rmesh,
                   std::vector<double> &pseudo_weight_rmesh);
void read_upf_local(std::ifstream &ifs, std::vector<double> &pseudo_local);
void read_upf_nonlocal(std::ifstream &ifs, const bool is_upf_ver2, 
                       std::vector<std::vector<double> > &pseudo_beta,
                       std::vector<int> &pseudo_lbeta,
                       Eigen::MatrixXd &pseudo_dij,
                       std::ostream *ost);
void read_upf_nlcc(std::ifstream &ifs, std::vector<double> &pseudo_nlcc);

} // namespace io_qe_files::upf

#endif // TC_IO_QE_FILES_UPF
