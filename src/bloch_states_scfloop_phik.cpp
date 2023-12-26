// [class BlochStates]
// wave functions, eigenvalues, filling, num. of bands, etc.

#include "include/header.hpp"

// private

void BlochStates::switch_left_right_orbitals_each(const std::string &calc_mode)
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
                    Eigen::VectorXcd &right_vec = phik_ref_[ispin][ik][iband][ispinor];
                    Eigen::VectorXcd &left_vec = phik_left_ref_[ispin][ik][iband][ispinor];

                    assert(right_vec.size() == left_vec.size());
                    Eigen::VectorXcd temp_vec = right_vec;
                    right_vec = left_vec;
                    left_vec = temp_vec;
                } // ispinor
            } // iband
        } // ik
    } // ispin

    // phik_old
    if (calc_mode=="SCF")
    {
        for (int ispin=0; ispin<phik_scf_old_.size(); ispin++)
        {
            for (int ik=0; ik<phik_scf_old_[ispin].size(); ik++)
            {
                for (int iband=0; iband<phik_scf_old_[ispin][ik].size(); iband++)
                {
                    for (int ispinor=0; ispinor<phik_scf_old_[ispin][ik][iband].size(); ispinor++)
                    {
                        Eigen::VectorXcd &right_vec = phik_scf_old_[ispin][ik][iband][ispinor];
                        Eigen::VectorXcd &left_vec = phik_left_scf_old_[ispin][ik][iband][ispinor];

                        assert(right_vec.size() == left_vec.size());
                        Eigen::VectorXcd temp_vec = right_vec;
                        right_vec = left_vec;
                        left_vec = temp_vec;
                    } // ispinor
                } // iband
            } // ik
        } // ispin
    }
}

// public

// non-collinear calculation not supported
void BlochStates::reset_phik(const int &ispin, const int &ik,
                             const std::vector<Eigen::VectorXcd> &V,
                             const std::string &right_or_left,
                             const bool mixes_density_matrix,
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

    // phik_scf_old_ for density-matrix mixing
    if (calc_mode=="SCF" && mixes_density_matrix)
    {
        std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phik_old_ref_ =
            right_or_left=="right" ? phik_scf_old_ : phik_left_scf_old_;

        phik_old_ref_[ispin][ik].resize(num_occupied_bands_[ispin][ik]);
        phik_old_ref_[ispin][ik].shrink_to_fit(); // memory saving
        for (int iband=0; iband<num_occupied_bands_[ispin][ik]; iband++)
        {
            phik_old_ref_[ispin][ik][iband].resize(1); // num_spinor
            phik_old_ref_[ispin][ik][iband][0] = phik_ref_[ispin][ik][iband][0];
        }
    } // SCF

    // phik_ref_
    for (int iband=0; iband<phik_ref_[ispin][ik].size(); iband++)
    {
        phik_ref_[ispin][ik][iband][0] = V[iband];
    }
}

void BlochStates::bcast_phik(const bool bcast_phik_left,
                             const bool mixes_density_matrix,
                             const std::string &calc_mode,
                             const bool am_i_mpi_rank0)
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

    if (calc_mode=="SCF" && mixes_density_matrix)
    {
        const int num_spinor = phik_ref_[0][0][0].size();
        for (int ispin=0; ispin<phik_scf_old_.size(); ispin++)
        {
            for (int ik=0; ik<phik_scf_old_[ispin].size(); ik++)
            {
                if (!am_i_mpi_rank0) // allocate before MPI_Bcast
                {
                    const int num_G_at_k = phik_ref_[ispin][ik][0][0].size();

                    phik_scf_old_[ispin][ik].resize(num_occupied_bands_[ispin][ik]);
                    phik_scf_old_[ispin][ik].shrink_to_fit(); // memory saving
                    for (int iband=0; iband<num_occupied_bands_[ispin][ik]; iband++)
                    {
                        phik_scf_old_[ispin][ik][iband].resize(num_spinor);
                        for (int ispinor=0; ispinor<num_spinor; ispinor++)
                        {
                            phik_scf_old_[ispin][ik][iband][ispinor].resize(num_G_at_k);
                        }
                    }

                    if (bcast_phik_left)
                    {
                        phik_left_scf_old_[ispin][ik].resize(num_occupied_bands_[ispin][ik]);
                        phik_left_scf_old_[ispin][ik].shrink_to_fit(); // memory saving
                        for (int iband=0; iband<num_occupied_bands_[ispin][ik]; iband++)
                        {
                            phik_left_scf_old_[ispin][ik][iband].resize(num_spinor);
                            for (int ispinor=0; ispinor<num_spinor; ispinor++)
                            {
                                phik_left_scf_old_[ispin][ik][iband][ispinor].resize(num_G_at_k);
                            }
                        }
                    } // bcast_phik_left
                } // !am_i_mpi_rank0
                
                for (int iband=0; iband<phik_scf_old_[ispin][ik].size(); iband++)
                {
                    for (int ispinor=0; ispinor<phik_scf_old_[ispin][ik][iband].size(); ispinor++)
                    {
                        MPI_Bcast(phik_scf_old_[ispin][ik][iband][ispinor].data(),
                                  phik_scf_old_[ispin][ik][iband][ispinor].size(),
                                  MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
                        if (bcast_phik_left)
                        {
                            MPI_Bcast(phik_left_scf_old_[ispin][ik][iband][ispinor].data(),
                                      phik_left_scf_old_[ispin][ik][iband][ispinor].size(),
                                      MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
                        }
                    } // ispinor
                } // iband
            } // ik
        } // ispin
    } // SCF
}

void BlochStates::switch_left_right_orbitals(const std::string &calc_mode)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");

    switch_left_right_orbitals_each("SCF"); // SCF orbitals are switched
    if (calc_mode=="BAND") { switch_left_right_orbitals_each("BAND"); } // BAND orbitals are switched
}
