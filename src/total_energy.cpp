// [class TotalEnergy]
// total energy, ewald energy, etc.

#include "include/header.hpp"

// private 

void TotalEnergy::resize_energy(const int &num_independent_spins,
                                const int &num_irreducible_kpoints_scf,
                                const std::vector<int> &num_bands_scf,
                                const std::string &calc_method)
{
    assert(num_independent_spins > 0);
    assert(num_irreducible_kpoints_scf > 0);
    assert(num_bands_scf.size() > 0);

    energy_1body_.resize(num_independent_spins);
    energy_2body_.resize(num_independent_spins);
    if (calc_method=="TC" || calc_method=="BITC")
    {
        energy_3body_.resize(num_independent_spins);
    }
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        energy_1body_[ispin].resize(num_irreducible_kpoints_scf);
        energy_2body_[ispin].resize(num_irreducible_kpoints_scf);
        if (calc_method=="TC" || calc_method=="BITC")
        {
            energy_3body_[ispin].resize(num_irreducible_kpoints_scf);
        }
        for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
        {
            energy_1body_[ispin][ik].resize(num_bands_scf[ispin]);
            energy_2body_[ispin][ik].resize(num_bands_scf[ispin]);
            if (calc_method=="TC" || calc_method=="BITC")
            {
                energy_3body_[ispin][ik].resize(num_bands_scf[ispin]);
            }
        }
    }
}

// public

void TotalEnergy::set_ewald_energy(const bool is_heg, const double &ewald_energy)
{
    if (!is_heg) { ewald_energy_ = ewald_energy; } // ewald_energy = 0 for is_heg==true
}

void TotalEnergy::bcast_qe_input(const int &num_independent_spins,
                                 const int &num_irreducible_kpoints,
                                 const std::vector<int> &num_bands_scf,
                                 const std::string &calc_mode, const std::string &calc_method,
                                 const bool is_heg)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");

    if (calc_mode=="SCF") { resize_energy(num_independent_spins, num_irreducible_kpoints, num_bands_scf, calc_method); }
    if (!is_heg) { MPI_Bcast(&ewald_energy_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); }
}

void TotalEnergy::reset_energies(const int &ispin, const int &ik,
                                 const Eigen::MatrixXcd &energy_1body,
                                 const Eigen::MatrixXcd &energy_2body,
                                 const Eigen::MatrixXcd &energy_3body,
                                 const bool uses_3body)
{
    assert(ispin >= 0 && ispin < energy_1body_.size());
    assert(ik >= 0 && ik < energy_1body_[ispin].size());

    const int num_bands_tc = energy_1body_[ispin][ik].size();
    assert(energy_1body.cols() >= num_bands_tc);
    assert(energy_1body.rows() >= num_bands_tc);
    assert(energy_2body.cols() >= num_bands_tc);
    assert(energy_2body.rows() >= num_bands_tc);
    assert(!uses_3body || energy_3body.cols() >= num_bands_tc);
    assert(!uses_3body || energy_3body.rows() >= num_bands_tc);

    for (int iband=0; iband<num_bands_tc; iband++)
    {
        energy_1body_[ispin][ik][iband] = energy_1body(iband, iband);
        energy_2body_[ispin][ik][iband] = energy_2body(iband, iband);
        if (uses_3body) { energy_3body_[ispin][ik][iband] = energy_3body(iband, iband); }
    }
}

void TotalEnergy::bcast_energies(const bool uses_3body)
{
    for (int ispin=0; ispin<energy_1body_.size(); ispin++)
    {
        for (int ik=0; ik<energy_1body_[ispin].size(); ik++)
        {
            MPI_Bcast(&energy_1body_[ispin][ik][0], energy_1body_[ispin][ik].size(),
                      MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&energy_2body_[ispin][ik][0], energy_2body_[ispin][ik].size(),
                      MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            if (uses_3body)
            {
                MPI_Bcast(&energy_3body_[ispin][ik][0], energy_3body_[ispin][ik].size(),
                          MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            }
        }
    }
}

void TotalEnergy::calc_total_energy(const Spin &spin,
                                    const Kpoints &kpoints,
                                    const std::vector<std::vector<std::vector<double> > > &filling,
                                    const std::vector<std::vector<std::vector<Complex> > > &eigenvalues_scf,
                                    const double &fermi_energy,
                                    const bool uses_3body,
                                    const bool am_i_mpi_rank0,
                                    std::ostream *ost)
{
    const int num_independent_spins = filling.size();
    const int num_irreducible_kpoints_scf = filling[0].size();
    const double spin_factor = (num_independent_spins==1 && !spin.is_spinor()) ? 2.0 : 1.0; // 2 for no-spin

    Complex total_energy_sigma0 = 0.0; // Total energy with sigma=0

    total_energy_1body_ = 0.0;
    total_energy_2body_ = 0.0;
    total_energy_3body_ = 0.0;
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
        {
            if (kpoints.kweight_scf()[ik] < 1e-8) { continue; } // zero k-weight (e.g. band k-point)

            const int num_bands_tc = energy_1body_[ispin][ik].size();
            for (int isym=0; isym<kpoints.kvectors_scf()[ik].size(); isym++)
            {
                if (!kpoints.is_time_reversal_used_at_k()[ik][isym])
                {
                    for (int iband=0; iband<num_bands_tc; iband++)
                    {
                        total_energy_1body_ += filling[ispin][ik][iband] * energy_1body_[ispin][ik][iband];
                        total_energy_2body_ += filling[ispin][ik][iband] * energy_2body_[ispin][ik][iband] / 2.0;
                        if (uses_3body) { total_energy_3body_ += filling[ispin][ik][iband] * energy_3body_[ispin][ik][iband] / 3.0; }

                        double filling_sigma0 = eigenvalues_scf[ispin][ik][iband].real() < fermi_energy + 1e-8 ?
                                                                                           spin_factor : 0.0; // fully-occupied or empty
                        total_energy_sigma0 += filling_sigma0 * (energy_1body_[ispin][ik][iband] +
                                                                 energy_2body_[ispin][ik][iband] / 2.0);
                        if (uses_3body) { total_energy_sigma0 += filling_sigma0 * energy_3body_[ispin][ik][iband] / 3.0; }
                    }
                }
                else
                {
                    for (int iband=0; iband<num_bands_tc; iband++)
                    {
                        total_energy_1body_ += filling[ispin][ik][iband] * std::conj(energy_1body_[ispin][ik][iband]);
                        total_energy_2body_ += filling[ispin][ik][iband] * std::conj(energy_2body_[ispin][ik][iband]) / 2.0;
                        if (uses_3body) { total_energy_3body_ += filling[ispin][ik][iband] * std::conj(energy_3body_[ispin][ik][iband]) / 3.0; }

                        double filling_sigma0 = eigenvalues_scf[ispin][ik][iband].real() < fermi_energy + 1e-8 ?
                                                                                           spin_factor : 0.0; // fully-occupied or empty
                        total_energy_sigma0 += filling_sigma0 * (std::conj(energy_1body_[ispin][ik][iband]) +
                                                                 std::conj(energy_2body_[ispin][ik][iband]) / 2.0);
                        if (uses_3body) { total_energy_sigma0 += filling_sigma0 * std::conj(energy_3body_[ispin][ik][iband]) / 3.0; }                    
                    }
                }
            }
        }
    }
    total_energy_1body_ /= kpoints.num_kpoints();
    total_energy_2body_ /= kpoints.num_kpoints();
    total_energy_3body_ /= kpoints.num_kpoints();
    total_energy_sigma0 /= kpoints.num_kpoints();
    total_energy_sigma0 += ewald_energy_;

    total_energy_difference_ = 
        ewald_energy_ + total_energy_1body_ + total_energy_2body_ + total_energy_3body_ - total_energy_;
    total_energy_ = ewald_energy_ + total_energy_1body_ + total_energy_2body_ + total_energy_3body_;

    if (am_i_mpi_rank0) 
    {
        *ost << "   Total energy (Ewald energy) = " << ewald_energy_ << " Ht." << std::endl; 
        *ost << "   Total energy (1-body terms) = " << total_energy_1body_ << " Ht." << std::endl; 
        *ost << "   Total energy (2-body terms) = " << total_energy_2body_ << " Ht." << std::endl; 
        if (uses_3body) { *ost << "   Total energy (3-body terms) = " << total_energy_3body_ << " Ht." << std::endl; }
        *ost << "   Total energy = " << total_energy_ << " Ht." << std::endl; 
        if (kpoints.smearing_mode()=="gaussian") { *ost << "   Total energy (sigma->0) = " << (total_energy_ + total_energy_sigma0)/2.0 << " Ht." << std::endl; }
        *ost << "   Total energy difference from the previous loop = " << total_energy_difference_ << " Ht." << std::endl; 
    }
}
