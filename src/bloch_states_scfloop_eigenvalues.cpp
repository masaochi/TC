// [class BlochStates]
// wave functions, eigenvalues, filling, num. of bands, etc.

#include "include/header.hpp"

void BlochStates::print_band_energies(const Kpoints &kpoints, std::ostream *ost) const
{
    *ost << "   Print band energies" << std::endl;
    for (int ispin=0; ispin<eigenvalues_band_.size(); ispin++)
    {
        *ost << "   Spin " << ispin << std::endl;
        for (int ik=0; ik<eigenvalues_band_[ispin].size(); ik++)
        {
            *ost << "   k-vector " << kpoints.kvectors_band()[ik][0].transpose() << std::endl;
            *ost << "   (band index, energy (eV))" << std::endl;
            for (int iband=0; iband<eigenvalues_band_[ispin][ik].size(); iband++)
            {
                *ost << "   " << iband << " " << eigenvalues_band_[ispin][ik][iband].real() * Ht_in_eV << std::endl;
            }
            *ost << std::endl;
        }
    }
}

void BlochStates::reset_eigenvalues(const int &ispin, const int &ik,
                                    const Eigen::VectorXcd &eigenvalues,
                                    const std::string &calc_mode)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    std::vector<std::vector<std::vector<Complex> > > &eigenvalues_ref_ =
        calc_mode=="SCF" ? eigenvalues_scf_ : eigenvalues_band_;

    assert(ispin >= 0 && ispin < eigenvalues_ref_.size());
    assert(ik >= 0 && ik < eigenvalues_ref_[ispin].size());
    assert(eigenvalues.size() >= eigenvalues_ref_[ispin][ik].size());

    for (int iband=0; iband<eigenvalues_ref_[ispin][ik].size(); iband++)
    {
        eigenvalues_ref_[ispin][ik][iband] = eigenvalues(iband);
    }
}

void BlochStates::bcast_eigenvalues(const std::string &calc_mode)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    std::vector<std::vector<std::vector<Complex> > > &eigenvalues_ref_ =
        calc_mode=="SCF" ? eigenvalues_scf_ : eigenvalues_band_;

    for (int ispin=0; ispin<eigenvalues_ref_.size(); ispin++)
    {
        for (int ik=0; ik<eigenvalues_ref_[ispin].size(); ik++)
        {
            MPI_Bcast(&eigenvalues_ref_[ispin][ik][0], eigenvalues_ref_[ispin][ik].size(),
                      MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
        }
    }
}
