// [class Parallelization]
// MPI variables

#include "include/header.hpp"

void Parallelization::set_mpi()
{
    MPI_Comm_rank(MPI_COMM_WORLD,&my_mpi_rank_);
    MPI_Comm_size(MPI_COMM_WORLD,&num_mpi_processes_);
    am_i_mpi_rank0_ = (my_mpi_rank_==0) ? true : false;
}

void Parallelization::print(std::ostream *ost) const
{
    *ost << "Num. of MPI processes = " << num_mpi_processes_ << std::endl;

#ifdef _OPENMP
    *ost << "Num. of OpenMP threads = " << omp_get_max_threads() << std::endl;    
#endif
}

// set is_assigned_...

// is_assigned_all_kpoints_occupied_bands_[ispin][ik][isym][iband]
//   ispin = spin (1 for no-spin, 2 for spin-polarized), ik = irreducible k-point (on the SCF k-mesh)
//   isym = index of symmetrically-equivalent k-point, iband = occupied band

// is_assigned_all_kpoints_occupied_bands_old_[ispin][ik][isym][iband]
//   same as above but for old wave functions (i.e., those obtained in the previous SCF loop)
//   necessary for density-matrix mixing (mixes_density_matrix && calc_mode==SCF)

// is_assigned_irreducible_kpoints_all_bands_[ispin][ik][iband]
//   ispin,ik = same as above but on the SCF or BAND k-mesh, iband = every band
void Parallelization::set_orbital_parallelization(const BlochStates &bloch_states, const Kpoints &kpoints,
                                                  const std::string &calc_mode,
                                                  const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi)
{
    assert(calc_mode=="SCF" || calc_mode=="BAND");
    const int num_independent_spins = bloch_states.filling().size();
    const int num_irreducible_kpoints_scf = kpoints.num_irreducible_kpoints_scf();
    const int num_irreducible_kpoints_ref = calc_mode=="SCF"
        ? num_irreducible_kpoints_scf : kpoints.num_irreducible_kpoints_band();

    int icount, ibegin, iend;
    std::vector<int> ncount_at_each_process(num_mpi_processes_);

    // (1) set is_assigned_all_kpoints_occupied_bands_[ispin][ik][isym][iband] and is_assigned_all_kpoints_occupied_bands_old_[ispin][ik][isym][iband]
    is_assigned_all_kpoints_occupied_bands_.resize(num_independent_spins);
    is_assigned_all_kpoints_occupied_bands_old_.resize(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++) // MPI processes are assigned for each "ispin" for efficient parallelization.
    {
        icount = 0;
        for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
        {
            for (int isym=0; isym<kpoints.kvectors_scf()[ik].size(); isym++)
            {
                // For mixes_density_matrix==false, filling_old()[ispin][ik].size()==0.
                // For mixes_density_matrix==true, this size equals the num. of occupied bands for old wave functions.
                const int nbands_old = bloch_states.filling_old()[ispin][ik].size();
                for (int iband=-nbands_old; iband<bloch_states.num_occupied_bands()[ispin][ik]; iband++)
                {
                    icount++;
                }
            }
        }

        // e.g., icount = 13, num_mpi_processes = 5
        //       each process is assigned (3, 3, 3, 2, 2) orbitals.
        std::fill(ncount_at_each_process.begin(), ncount_at_each_process.end(), icount/num_mpi_processes_); // (2,2,2,2,2)
        for (int iproc=0; iproc<(icount%num_mpi_processes_); iproc++) 
        {
            ncount_at_each_process[iproc]++; // (3,3,3,2,2)
        }
        //       (ibegin, iend) = (0, 2) for my_mpi_rank_==0, (3, 5) for rank_1, (6, 8) for rank_2,
        //                        (9, 10) for rank_3, (11, 12) for rank_4.
        ibegin = 0;
        for (int iproc=0; iproc<my_mpi_rank_; iproc++)
        {
            ibegin += ncount_at_each_process[iproc];
        }
        iend = ibegin + ncount_at_each_process[my_mpi_rank_] - 1;
        
        // set 
        icount = 0;

        is_assigned_all_kpoints_occupied_bands_[ispin].resize(num_irreducible_kpoints_scf);
        is_assigned_all_kpoints_occupied_bands_old_[ispin].resize(num_irreducible_kpoints_scf);
        for (int ik=0; ik<num_irreducible_kpoints_scf; ik++)
        {
            is_assigned_all_kpoints_occupied_bands_[ispin][ik].resize(kpoints.kvectors_scf()[ik].size());
            is_assigned_all_kpoints_occupied_bands_old_[ispin][ik].resize(kpoints.kvectors_scf()[ik].size());
            for (int isym=0; isym<kpoints.kvectors_scf()[ik].size(); isym++)
            {
                is_assigned_all_kpoints_occupied_bands_[ispin][ik][isym].resize(bloch_states.num_occupied_bands()[ispin][ik]);
                const int nbands_old = bloch_states.filling_old()[ispin][ik].size();
                is_assigned_all_kpoints_occupied_bands_old_[ispin][ik][isym].resize(nbands_old);
                for (int iband=-nbands_old; iband<bloch_states.num_occupied_bands()[ispin][ik]; iband++)
                {
                    if (icount>=ibegin && icount<=iend)
                    {
                        if (iband>=0)
                        {
                            is_assigned_all_kpoints_occupied_bands_[ispin][ik][isym][iband] = true;
                        }
                        else
                        {
                            is_assigned_all_kpoints_occupied_bands_old_[ispin][ik][isym][-1-iband] = true;
                        }
                    }
                    else
                    {
                        if (iband>=0)
                        {
                            is_assigned_all_kpoints_occupied_bands_[ispin][ik][isym][iband] = false;
                        }
                        else
                        {
                            is_assigned_all_kpoints_occupied_bands_old_[ispin][ik][isym][-1-iband] = false;
                        }
                    }
                    icount++;
                }
            }
        }
    } // ispin

    // (2) set is_assigned_irreducible_kpoints_all_bands_[ispin][ik][iband]
    is_assigned_irreducible_kpoints_all_bands_.resize(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++) // MPI processes are assigned for each "ispin" for efficient parallelization.
    {
        icount = 0;
        for (int ik=0; ik<num_irreducible_kpoints_ref; ik++)
        {
            for (int iband=0; iband<phi[ispin][ik].size(); iband++)
            {
                icount++;
            }
        }
    
        std::fill(ncount_at_each_process.begin(), ncount_at_each_process.end(), icount/num_mpi_processes_);
        for (int iproc=0; iproc<(icount%num_mpi_processes_); iproc++) 
        {
            ncount_at_each_process[iproc]++;
        }
        ibegin = 0;
        for (int iproc=0; iproc<my_mpi_rank_; iproc++)
        {
            ibegin += ncount_at_each_process[iproc];
        }
        iend = ibegin + ncount_at_each_process[my_mpi_rank_] - 1;
        
        // set 
        icount = 0;

        is_assigned_irreducible_kpoints_all_bands_[ispin].resize(num_irreducible_kpoints_ref);
        for (int ik=0; ik<num_irreducible_kpoints_ref; ik++)
        {
            is_assigned_irreducible_kpoints_all_bands_[ispin][ik].resize(phi[ispin][ik].size());
            for (int iband=0; iband<phi[ispin][ik].size(); iband++)
            {
                if (icount>=ibegin && icount<=iend)
                {
                    is_assigned_irreducible_kpoints_all_bands_[ispin][ik][iband] = true;
                }
                else
                {
                    is_assigned_irreducible_kpoints_all_bands_[ispin][ik][iband] = false;
                }
                icount++;
            }
        }
    } // ispin
}
