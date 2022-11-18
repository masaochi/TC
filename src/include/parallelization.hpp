// [class Parallelization]
// MPI variables

#ifndef TC_PARALLELIZATION_HPP
#define TC_PARALLELIZATION_HPP

class Parallelization
{
private:
    int my_mpi_rank_; // MPI rank (0,1,2,...)
    int num_mpi_processes_; // no. of MPI Processes
    bool am_i_mpi_rank0_; // true for rank==0

    // orbital parallelization for SCF
    // [ispin][ik][isym][iband] ispin = spin (1 for no-spin, 2 for spin-polarized), ik = irreducible k-point
    //                          isym = index of symmetrically-equivalent k-point, iband = occupied band
    std::vector<std::vector<std::vector<std::vector<bool> > > > is_assigned_all_kpoints_occupied_bands_; 
    // [ispin][ik][iband] ispin,ik = same as above, iband = every band
    std::vector<std::vector<std::vector<bool> > > is_assigned_irreducible_kpoints_all_bands_;

    void set_mpi(); // called by constructor

public:
    int my_mpi_rank() const { return my_mpi_rank_; }
    int num_mpi_processes() const { return num_mpi_processes_; }
    bool am_i_mpi_rank0() const { return am_i_mpi_rank0_; }
    const std::vector<std::vector<std::vector<std::vector<bool> > > > &is_assigned_all_kpoints_occupied_bands() const
    { return is_assigned_all_kpoints_occupied_bands_; }
    const std::vector<std::vector<std::vector<bool> > > &is_assigned_irreducible_kpoints_all_bands() const
    { return is_assigned_irreducible_kpoints_all_bands_; }

    Parallelization() { set_mpi(); }

    void print(std::ostream *ost) const;
    // set is_assigned_...
    void set_orbital_parallelization(const BlochStates &bloch_states, const Kpoints &kpoints,
                                     const std::string &calc_mode,
                                     const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi);
};


#endif // TC_PARALLELIZATION_HPP

