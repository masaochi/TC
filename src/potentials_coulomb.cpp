// [class Potentials]
// pseudopotential, Jastrow, Ewald sum

// given in atomic unit (Hartree)

#include "include/header.hpp"

// For divergence correction (Coulomb potential)
void Potentials::set_sum_of_Vaux(const Parallelization &parallelization,
                                 const Method &method,
                                 const CrystalStructure &crystal_structure,
                                 const Kpoints &kpoints)
{
    if (method.calc_method()=="FREE") { return; } // free-electron mode: no interaction

    // not used in HF-BAND.
    if (!(method.calc_method()=="HF" && method.calc_mode()=="BAND"))
    {
        set_sum_of_Vaux_each(parallelization,
                             crystal_structure,
                             kpoints,
                             kpoints.kvectors_scf(),
                             sum_of_Vaux_scf_);

        // used only when (BI)TC && SCF
        if ((method.calc_method()=="TC" || method.calc_method()=="BITC") && method.calc_mode()=="SCF")
        {
            set_sum_of_kVaux_each(parallelization,
                                  crystal_structure,
                                  kpoints,
                                  kpoints.kvectors_scf(),
                                  sum_of_kVaux_scf_);
        }
    }
    if (method.calc_mode()=="BAND") 
    {
        set_sum_of_Vaux_each(parallelization,
                             crystal_structure,
                             kpoints,
                             kpoints.kvectors_band(),
                             sum_of_Vaux_band_);

        if (method.calc_method()=="TC" || method.calc_method()=="BITC")
        {
            // used only when (BI)TC && BAND
            set_sum_of_kVaux_each(parallelization,
                                  crystal_structure,
                                  kpoints,
                                  kpoints.kvectors_band(),
                                  sum_of_kVaux_band_);
        }
    }
}

// private, called in set_sum_of_Vaux
void Potentials::set_sum_of_Vaux_each(const Parallelization &parallelization,
                                      const CrystalStructure &crystal_structure,
                                      const Kpoints &kpoints,
                                      const std::vector<std::vector<Eigen::Vector3d> > &kvectors,
                                      std::vector<double> &sum_of_Vaux_)
{
    const int num_irreducible_kpoints = kvectors.size(); // num_irreducible_kpoints_(scf/band)
    const int num_irreducible_kpoints_scf = kpoints.num_irreducible_kpoints_scf();

    alpha_Vaux_ = std::cbrt(crystal_structure.unit_cell_volume())/(2.0*PI);
//    alpha_Vaux_ = (alpha_Vaux_ * alpha_Vaux_)/2.0; // alpha = 0.5*V^(2/3)/(4*PI^2)
    alpha_Vaux_ = (alpha_Vaux_ * alpha_Vaux_); // alpha = V^(2/3)/(4*PI^2)

    // set parallelization index (only used here) cf. Parallelization::set_orbital_parallelization()
    std::vector<bool> if_match(num_irreducible_kpoints);
    int ibegin, iend;
    set_parallelization_in_sum_Vaux(parallelization, kpoints, kvectors, if_match, ibegin, iend);

    // calculate
    bool if_match_kq;
    std::vector<double> sum_of_Vaux_local(num_irreducible_kpoints);
    for (int ik=0; ik<num_irreducible_kpoints; ik++)
    {
        sum_of_Vaux_local[ik] = 0.0;

        int icount = 0;
        for (int iq=0; iq<num_irreducible_kpoints_scf; iq++)
        {
            for (int isymq=0; isymq<kpoints.kvectors_scf()[iq].size(); isymq++)
            {
                if (kpoints.kweight_scf()[iq] < 1e-8) { continue; } // zero k-weight (i.e. band k-weight)
                Eigen::Vector3d kq = kvectors[ik][0] - kpoints.kvectors_scf()[iq][isymq];
                for (int idim=0; idim<3; idim++)
                {
                    while (kq(idim) < -1e-8) { kq(idim) += 1.0; }
                    while (kq(idim) > 1.0-1e-8) { kq(idim) -= 1.0; }
                }
                if (kq.squaredNorm()<1e-8) 
                {
                    if_match_kq = true;
                }
                else
                {
                    if_match_kq = false;
                }

                icount++;
                if ((icount-1)<ibegin || (icount-1)>iend) { continue; } // not assigned to this MPI process

                if (!if_match_kq) // nshell=0 contribution
                {
                    Eigen::Vector3d kqvect = crystal_structure.reciprocal_vectors().transpose()*kq;
                    double kq2 = kqvect.squaredNorm();
                    sum_of_Vaux_local[ik] += std::exp(-alpha_Vaux_*kq2)/kq2;
                }
                for (int nshell=1; ; nshell++)
                {
                    double sum_on_shell = 0.0;
                    for (int Gx=-nshell; Gx<=+nshell; Gx++)
                    { 
                        for (int Gy=-nshell; Gy<=+nshell; Gy++)
                        {
                            for (int Gz=-nshell; Gz<=+nshell; Gz++)
                            {
                                if (Gx!=nshell && -Gx!=nshell && Gy!=nshell && -Gy!=nshell
                                    && Gz!=nshell && -Gz!=nshell) { continue; }
                                Eigen::Vector3i tmpG = {Gx, Gy, Gz};
                                Eigen::Vector3d kqGvect = crystal_structure.reciprocal_vectors().transpose() 
                                    * (kq + tmpG.cast<double>());
                                double kqG2 = kqGvect.squaredNorm();
                                sum_on_shell += std::exp(-alpha_Vaux_*kqG2)/kqG2;
                            } // Gz
                        } // Gy
                    } // Gx
                    sum_of_Vaux_local[ik] += sum_on_shell;
                    if (sum_on_shell < 1e-20) { break; }
                } // int nshell
            } // isymq
        } // iq
        sum_of_Vaux_local[ik] *= FourPI;
    } // ik

    // allreduce
    sum_of_Vaux_.resize(num_irreducible_kpoints);
    MPI_Allreduce(&sum_of_Vaux_local[0],
                  &sum_of_Vaux_[0],
                  num_irreducible_kpoints,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int ik=0; ik<num_irreducible_kpoints; ik++)
    {
        if (if_match[ik]) // 4pi*exp[-alpha q^2]/q^2 ~ 4pi/q^2 - 4pi*alpha at q=0 
                          // (Note: if_match==true, both nshell=0 contribution in Vaux and 1/0^2 in coulomb are not included)
        {
            // This is not the correction for \sum_G 4pi*exp[-alpha q^2]/q^2,
            // rather the correction in \int 4pi/q^2 ~ \int 4pi*exp[-alpha q^2]/q^2 + 4pi*alpha for \int dq near q=0.
            // In calc_hamiltonian_two_body_x, we calculate \int - sum_of_Vaux. 
            // Thus, "-4pi*alpha" correction for sum_of_Vaux = "+4pi*alpha" correction for \int.
            sum_of_Vaux_[ik] -= FourPI*alpha_Vaux_;
        } 
    }
}

void Potentials::set_sum_of_kVaux_each(const Parallelization &parallelization,
                                       const CrystalStructure &crystal_structure,
                                       const Kpoints &kpoints,
                                       const std::vector<std::vector<Eigen::Vector3d> > &kvectors,
                                       std::vector<Eigen::Vector3d> &sum_of_kVaux_)
{
    const int num_irreducible_kpoints = kvectors.size(); // num_irreducible_kpoints_(scf/band)
    const int num_irreducible_kpoints_scf = kpoints.num_irreducible_kpoints_scf();

    // alpha_Vaux_ is already set in set_sum_of_Vaux_each()
    assert(std::abs(alpha_Vaux_) > 1e-8);

    // set parallelization index (only used here) cf. Parallelization::set_orbital_parallelization()
    std::vector<bool> if_match(num_irreducible_kpoints);
    int ibegin, iend;
    set_parallelization_in_sum_Vaux(parallelization, kpoints, kvectors, if_match, ibegin, iend);

    sum_of_kVaux_.resize(num_irreducible_kpoints);

    // calculate
    bool if_match_kq;
    Eigen::Vector3d sum_of_kVaux_local;
    for (int ik=0; ik<num_irreducible_kpoints; ik++)
    {
        if (if_match[ik]) // if_match=true then sum_of_kVaux = 0
        {
            sum_of_kVaux_[ik] = {0.0, 0.0, 0.0};
            continue; 
        }

        int icount = 0;
        sum_of_kVaux_local = {0.0, 0.0, 0.0};
        for (int iq=0; iq<num_irreducible_kpoints_scf; iq++)
        {
            for (int isymq=0; isymq<kpoints.kvectors_scf()[iq].size(); isymq++)
            {
                if (kpoints.kweight_scf()[iq] < 1e-8) { continue; } // zero k-weight (i.e. band k-weight)
                Eigen::Vector3d kq = kvectors[ik][0] - kpoints.kvectors_scf()[iq][isymq];
                for (int idim=0; idim<3; idim++)
                {
                    while (kq(idim) < -1e-8) { kq(idim) += 1.0; }
                    while (kq(idim) > 1.0-1e-8) { kq(idim) -= 1.0; }
                }
                if (kq.squaredNorm()<1e-8) 
                {
                    if_match_kq = true;
                }
                else
                {
                    if_match_kq = false;
                }

                icount++;
                if ((icount-1)<ibegin || (icount-1)>iend) { continue; } // not assigned to this MPI process

                if (!if_match_kq) // nshell=0 contribution
                {
                    Eigen::Vector3d kqvect = crystal_structure.reciprocal_vectors().transpose()*kq;
                    double kq2 = kqvect.squaredNorm();
                    sum_of_kVaux_local += (std::exp(-alpha_Vaux_*kq2)/kq2) * kqvect;
                }
                for (int nshell=1; ; nshell++)
                {
                    Eigen::Vector3d sum_on_shell = {0.0, 0.0, 0.0};
                    for (int Gx=-nshell; Gx<=+nshell; Gx++)
                    { 
                        for (int Gy=-nshell; Gy<=+nshell; Gy++)
                        {
                            for (int Gz=-nshell; Gz<=+nshell; Gz++)
                            {
                                if (Gx!=nshell && -Gx!=nshell && Gy!=nshell && -Gy!=nshell
                                    && Gz!=nshell && -Gz!=nshell) { continue; }
                                Eigen::Vector3i tmpG = {Gx, Gy, Gz};
                                Eigen::Vector3d kqGvect = crystal_structure.reciprocal_vectors().transpose() 
                                    * (kq + tmpG.cast<double>());
                                double kqG2 = kqGvect.squaredNorm();
                                sum_on_shell += (std::exp(-alpha_Vaux_*kqG2)/kqG2) * kqGvect;
                            } // Gz
                        } // Gy
                    } // Gx
                    sum_of_kVaux_local += sum_on_shell;
                    if (sum_on_shell.squaredNorm() < 1e-20) { break; }
                } // int nshell
            } // isymq
        } // iq
        sum_of_kVaux_local *= FourPI;

        MPI_Allreduce(sum_of_kVaux_local.data(),
                      sum_of_kVaux_[ik].data(),
                      3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    } // ik
}

// set parallelization index (used only in set_sum_of_Vaux_each() and set_sum_of_kVaux_band)
// cf. Parallelization::set_orbital_parallelization()
void Potentials::set_parallelization_in_sum_Vaux(const Parallelization &parallelization,
                                                 const Kpoints &kpoints,
                                                 const std::vector<std::vector<Eigen::Vector3d> > &kvectors,
                                                 std::vector<bool> &if_match,
                                                 int &ibegin, int &iend)
{
    const int num_irreducible_kpoints = kvectors.size(); // num_irreducible_kpoints_(scf/band)
    const int num_irreducible_kpoints_scf = kpoints.num_irreducible_kpoints_scf();
    const int num_mpi_processes = parallelization.num_mpi_processes();
    const int my_mpi_rank = parallelization.my_mpi_rank();

    for (int ik=0; ik<num_irreducible_kpoints; ik++)
    {
        if_match[ik] = false;
        for (int iq=0; iq<num_irreducible_kpoints_scf; iq++)
        {
            for (int isymq=0; isymq<kpoints.kvectors_scf()[iq].size(); isymq++)
            {
                if (kpoints.kweight_scf()[iq] < 1e-8) { continue; } // zero k-weight (i.e. band k-weight)
                Eigen::Vector3d kq = kvectors[ik][0] - kpoints.kvectors_scf()[iq][isymq];
                for (int idim=0; idim<3; idim++)
                {
                    while (kq(idim) < -1e-8) { kq(idim) += 1.0; }
                    while (kq(idim) > 1.0-1e-8) { kq(idim) -= 1.0; }
                }
                if (kq.squaredNorm()<1e-8) { if_match[ik] = true; } // always true for the SCF k-mesh
            }
        }
    }

    int icount = 0;
    for (int iq=0; iq<num_irreducible_kpoints_scf; iq++)
    {
        for (int isymq=0; isymq<kpoints.kvectors_scf()[iq].size(); isymq++)
        {
            if (kpoints.kweight_scf()[iq] < 1e-8) { continue; } // zero k-weight (i.e. band k-weight)
            icount++;
        }
    }

    std::vector<int> ncount_each(num_mpi_processes);
    std::fill(ncount_each.begin(), ncount_each.end(), icount/num_mpi_processes);
    for (int iproc=0; iproc<(icount%num_mpi_processes); iproc++) 
    {
        ncount_each[iproc]++;
    }
    ibegin = 0;
    for (int iproc=0; iproc<my_mpi_rank; iproc++)
    {
        ibegin += ncount_each[iproc];
    }
    iend = ibegin + ncount_each[my_mpi_rank] - 1;
}
