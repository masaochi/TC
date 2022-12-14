// [namespace read_upf]
// Read upf (pseudopotential) files and MPI_Bcast

#ifndef TC_READ_UPF_HPP
#define TC_READ_UPF_HPP

namespace read_upf
{

void read_upf(FileNames &file_names,
              CrystalStructure &crystal_structure,
              Potentials &potentials,
              PlaneWaveBasis &plane_wave_basis,
              const bool am_i_mpi_rank0,
              std::ostream *ost);
    
bool check_upf_ver2(std::ifstream &ifs, std::ostream *ost);    
void read_upf_header(std::ifstream &ifs, const bool is_upf_ver2, 
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

} // namespace read_upf

#endif // TC_READ_UPF
