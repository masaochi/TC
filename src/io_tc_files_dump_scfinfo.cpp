// [namespace io_tc_files]
// Read input.in, read/write tc_wfc.dat (wave functions) and tc_energy.dat (orbital energies), etc.

#include "include/header.hpp"

// dump SCF information for BAND calculation
void io_tc_files::dump_scfinfo(const FileNames &file_names,
                               const Kpoints &kpoints,
                               const PlaneWaveBasis &plane_wave_basis,
                               const BlochStates &bloch_states,
                               std::ostream *ost)
{
    *ost << " Dump SCF informations (" << file_names.tc_scfinfo() << ")" << std::endl;
    std::ofstream ofs(file_names.tc_scfinfo(), std::ios::out | std::ios::binary); // binary
    if (ofs.fail()) { error_messages::cannot_open(file_names.tc_scfinfo()); }

    // basic information
    const int num_independent_spins = bloch_states.phik_scf().size();
    ofs.write(reinterpret_cast<const char*>(&num_independent_spins), sizeof(int));
    const int num_irreducible_kpoints_scf = kpoints.num_irreducible_kpoints_scf();
    ofs.write(reinterpret_cast<const char*>(&num_irreducible_kpoints_scf), sizeof(int));

    std::vector<int> num_sym_at_k(num_irreducible_kpoints_scf); // number of symmetrically equivalent k-points for each irreducible k-point
    for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
    {
        num_sym_at_k[ik] = kpoints.kvectors_scf()[ik].size();
        ofs.write(reinterpret_cast<const char*>(&num_sym_at_k[ik]), sizeof(int));
    }

    std::vector<int> size_FFT_grid_vec = plane_wave_basis.size_FFT_grid_vec(); // to check consistency of the FFT grid
    ofs.write(reinterpret_cast<const char*>(&size_FFT_grid_vec[0]), sizeof(int)*3);

    // SCF variables in class BlochStates
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        const int num_bands_scf = bloch_states.num_bands_scf()[ispin];
        ofs.write(reinterpret_cast<const char*>(&num_bands_scf), sizeof(int));
    }

    // SCF variables in class Kpoints
    for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
    {
        double kweight_scf = kpoints.kweight_scf()[ik];
        ofs.write(reinterpret_cast<const char*>(&kweight_scf), sizeof(double)); 
    }
    for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
    {
        for (int isym=0; isym<num_sym_at_k[ik]; isym++)
        {
            for (int idim=0; idim<3; idim++)
            {
                double kvectors_scf_i = kpoints.kvectors_scf()[ik][isym](idim);
                ofs.write(reinterpret_cast<const char*>(&kvectors_scf_i), sizeof(double));
            }
            int index_of_rotation_at_k = kpoints.index_of_rotation_at_k()[ik][isym];
            bool is_time_reversal_used_at_k = kpoints.is_time_reversal_used_at_k()[ik][isym];
            ofs.write(reinterpret_cast<const char*>(&index_of_rotation_at_k), sizeof(int));
            ofs.write(reinterpret_cast<const char*>(&is_time_reversal_used_at_k), sizeof(bool));
        }
    }

    // SCF varialbes in class PlaneWaveBasis
    std::vector<int> num_G_at_k_scf = plane_wave_basis.num_G_at_k_scf();
    for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
    {
        ofs.write(reinterpret_cast<const char*>(&num_G_at_k_scf[ik]), sizeof(int)); 
    }
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
        {
            for (int isym=0; isym<num_sym_at_k[ik]; isym++)
            {
                Eigen::VectorXi Gindex_at_k_scf = plane_wave_basis.Gindex_at_k_scf()[ispin][ik][isym];
                Eigen::VectorXcd phase_factor_at_k = plane_wave_basis.phase_factor_at_k()[ispin][ik][isym];
                for (int ipw_at_k=0; ipw_at_k<num_G_at_k_scf[ik]; ipw_at_k++)
                {
                    ofs.write(reinterpret_cast<const char*>(&Gindex_at_k_scf(ipw_at_k)), sizeof(int));

                    double re = phase_factor_at_k(ipw_at_k).real();
                    double im = phase_factor_at_k(ipw_at_k).imag();
                    ofs.write(reinterpret_cast<const char*>(&re), sizeof(double));
                    ofs.write(reinterpret_cast<const char*>(&im), sizeof(double));
                }
            }
        }
    }

    ofs.close();
}
