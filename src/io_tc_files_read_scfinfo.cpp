// [namespace io_tc_files]
// Read input.in, read/write tc_wfc.dat (wave functions) and tc_energy.dat (orbital energies), etc.

#include "include/header.hpp"

// read SCF information for BAND calculation
void io_tc_files::set_scfinfo(const FileNames &file_names,
                              const Method &method,
                              const Spin &spin,
                              Kpoints &kpoints,
                              PlaneWaveBasis &plane_wave_basis,
                              BlochStates &bloch_states,
                              const bool am_i_mpi_rank0, 
                              std::ostream *ost)
{
    if (am_i_mpi_rank0) { read_scfinfo(file_names, spin, kpoints, plane_wave_basis, bloch_states, ost); }
    bcast_scfinfo(kpoints, plane_wave_basis, bloch_states, spin.is_spinor(), method.calc_method(), am_i_mpi_rank0);
}

void io_tc_files::read_scfinfo(const FileNames &file_names,
                               const Spin &spin,
                               Kpoints &kpoints,
                               PlaneWaveBasis &plane_wave_basis,
                               BlochStates &bloch_states,
                               std::ostream *ost)
{
    *ost << " Read SCF informations (" << file_names.tc_scfinfo() << ")" << std::endl;
    std::ifstream ifs(file_names.tc_scfinfo(), std::ios::in | std::ios::binary); // binary
    if (ifs.fail()) { error_messages::cannot_open(file_names.tc_scfinfo()); }

    // basic information
    int num_independent_spins, num_irreducible_kpoints_scf;
    ifs.read(reinterpret_cast<char*>(&num_independent_spins), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&num_irreducible_kpoints_scf), sizeof(int));
    if (num_independent_spins != spin.num_independent_spins()) { error_messages::stop("num_independent_spins is inconsistent (read_scfinfo)"); }

    std::vector<int> num_sym_at_k(num_irreducible_kpoints_scf);
    for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
    {
        ifs.read(reinterpret_cast<char*>(&num_sym_at_k[ik]), sizeof(int));
    }

    std::vector<int> size_FFT_grid_vec(3);
    ifs.read(reinterpret_cast<char*>(&size_FFT_grid_vec[0]), sizeof(int)*3);
    for (int idim=0; idim<3; idim++)
    {
        if (size_FFT_grid_vec[idim] != plane_wave_basis.size_FFT_grid_vec()[idim]) { error_messages::stop("inconsistent FFT grid between SCF and BAND."); }
    }

    // SCF variables in class BlochStates
    std::vector<int> num_bands_scf(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        ifs.read(reinterpret_cast<char*>(&num_bands_scf[ispin]), sizeof(int));
    }
    bloch_states.set_num_bands_tc(num_bands_scf, "SCF");

    // SCF variables in class Kpoints
    std::vector<double> kweight_scf(num_irreducible_kpoints_scf);
    for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
    {
        ifs.read(reinterpret_cast<char*>(&kweight_scf[ik]), sizeof(double)); 
    }
    std::vector<std::vector<Eigen::Vector3d> > kvectors_scf(num_irreducible_kpoints_scf);
    std::vector<std::vector<int> > index_of_rotation_at_k(num_irreducible_kpoints_scf);
    std::vector<std::vector<bool> > is_time_reversal_used_at_k(num_irreducible_kpoints_scf);
    for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
    {
        kvectors_scf[ik].resize(num_sym_at_k[ik]);
        index_of_rotation_at_k[ik].resize(num_sym_at_k[ik]);
        is_time_reversal_used_at_k[ik].resize(num_sym_at_k[ik]);
        for (int isym=0; isym<num_sym_at_k[ik]; isym++)
        {
            for (int idim=0; idim<3; idim++)
            {
                ifs.read(reinterpret_cast<char*>(&kvectors_scf[ik][isym](idim)), sizeof(double));
            }
            ifs.read(reinterpret_cast<char*>(&index_of_rotation_at_k[ik][isym]), sizeof(int));
            bool bool_tmp;
            ifs.read(reinterpret_cast<char*>(&bool_tmp), sizeof(bool));
            is_time_reversal_used_at_k[ik][isym] = bool_tmp;
        }
    }
    kpoints.set_scfinfo(kweight_scf, kvectors_scf, index_of_rotation_at_k, is_time_reversal_used_at_k);

    // SCF varialbes in class PlaneWaveBasis
    std::vector<int> num_G_at_k_scf(num_irreducible_kpoints_scf);
    for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
    {
        ifs.read(reinterpret_cast<char*>(&num_G_at_k_scf[ik]), sizeof(int)); 
    }
    std::vector<std::vector<std::vector<Eigen::VectorXi> > > Gindex_at_k_scf(num_independent_spins);
    std::vector<std::vector<std::vector<Eigen::VectorXcd> > > phase_factor_at_k(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        Gindex_at_k_scf[ispin].resize(num_irreducible_kpoints_scf);
        phase_factor_at_k[ispin].resize(num_irreducible_kpoints_scf);
        for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
        {
            Gindex_at_k_scf[ispin][ik].resize(num_sym_at_k[ik]);
            phase_factor_at_k[ispin][ik].resize(num_sym_at_k[ik]);
            for (int isym=0; isym<num_sym_at_k[ik]; isym++)
            {
                Gindex_at_k_scf[ispin][ik][isym].resize(num_G_at_k_scf[ik]);
                phase_factor_at_k[ispin][ik][isym].resize(num_G_at_k_scf[ik]);
                for (int ipw_at_k=0; ipw_at_k<num_G_at_k_scf[ik]; ipw_at_k++)
                {
                    ifs.read(reinterpret_cast<char*>(&Gindex_at_k_scf[ispin][ik][isym](ipw_at_k)), sizeof(int));
                    double re,im;
                    ifs.read(reinterpret_cast<char*>(&re), sizeof(double));
                    ifs.read(reinterpret_cast<char*>(&im), sizeof(double));
                    phase_factor_at_k[ispin][ik][isym](ipw_at_k) = re + I*im;
                }
            }
        }
    }
    plane_wave_basis.set_scfinfo(num_G_at_k_scf, Gindex_at_k_scf, phase_factor_at_k);

    ifs.close();
}

void io_tc_files::bcast_scfinfo(Kpoints &kpoints,
                                PlaneWaveBasis &plane_wave_basis,
                                BlochStates &bloch_states,
                                const bool is_spinor,
                                const std::string &calc_method,
                                const bool am_i_mpi_rank0)
{
    kpoints.bcast_qe_input("SCF", am_i_mpi_rank0);
    plane_wave_basis.bcast("SCF", am_i_mpi_rank0, false); // false = FFT grid was already setup
    bloch_states.bcast_scfinfo(plane_wave_basis, is_spinor, calc_method, am_i_mpi_rank0); // also allocate several varialbes such as phik_scf
}
