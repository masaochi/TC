// [class BlochStates]
// wave functions, eigenvalues, filling, num. of bands, etc.

#include "include/header.hpp"

// non-collinear calculation not supported
void BlochStates::reset_phik(const int &ispin, const int &ik,
                             const std::vector<Eigen::VectorXcd> &V,
                             const std::string &right_or_left,
                             const std::string &calc_mode)
{
    assert(right_or_left=="right" || right_or_left=="left");
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_ref_ =
        (right_or_left=="right" && calc_mode=="SCF") ? phik_scf_ :
        (right_or_left=="right" && calc_mode=="BAND") ? phik_band_ :
        (right_or_left=="left" && calc_mode=="SCF") ? phik_left_scf_ : phik_left_band_;

    assert(ispin >= 0 && ispin < phik_ref_.size());
    assert(ik >= 0 && ik < phik_ref_[ispin].size());

    const int num_spinor = phik_ref_[ispin][ik][0].size();
    assert(num_spinor==1);

    assert(phik_ref_[ispin][ik].size() <= V.size()); // band index
    assert(phik_ref_[ispin][ik][0][0].size() == V[0].size());

    for (int iband=0; iband<phik_ref_[ispin][ik].size(); iband++)
    {
        phik_ref_[ispin][ik][iband][0] = V[iband];
    }
}

void BlochStates::bcast_phik(const bool bcast_phik_left,
                             const std::string &calc_mode)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_ref_ =
        calc_mode=="SCF" ? phik_scf_ : phik_band_;
    std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_left_ref_ =
        calc_mode=="SCF" ? phik_left_scf_ : phik_left_band_;

    for (int ispin=0; ispin<phik_ref_.size(); ispin++)
    {
        for (int ik=0; ik<phik_ref_[ispin].size(); ik++)
        {
            for (int iband=0; iband<phik_ref_[ispin][ik].size(); iband++)
            {
                for (int ispinor=0; ispinor<phik_ref_[ispin][ik][iband].size(); ispinor++)
                {
                    MPI_Bcast(phik_ref_[ispin][ik][iband][ispinor].data(),
                              phik_ref_[ispin][ik][iband][ispinor].size(),
                              MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
                    if (bcast_phik_left)
                    {
                        MPI_Bcast(phik_left_ref_[ispin][ik][iband][ispinor].data(),
                                  phik_left_ref_[ispin][ik][iband][ispinor].size(),
                                  MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
                    }
                }
            }
        }
    }
}
