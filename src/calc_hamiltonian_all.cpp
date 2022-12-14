// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

void calc_hamiltonian::all(const Parallelization &parallelization, 
                           MyClock &my_clock, const Method &method,
                           const CrystalStructure &crystal_structure,
                           const Symmetry &symmetry, const Potentials &potentials,
                           const Spin &spin, const Kpoints &kpoints,
                           PlaneWaveBasis &plane_wave_basis, 
                           const BlochStates &bloch_states, 
                           const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
                           std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H1phi,
                           std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H2phi,
                           std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H3phi,
                           std::ostream *ost)
{
    const bool am_i_mpi_rank0 = parallelization.am_i_mpi_rank0();
    const int num_independent_spins = spin.num_independent_spins();
    const int num_irreducible_kpoints = method.calc_mode()=="SCF" ? 
        kpoints.num_irreducible_kpoints_scf() : kpoints.num_irreducible_kpoints_band();
    const int num_spinor = (spin.is_spinor()==false) ? 1 : 2;

    if (am_i_mpi_rank0) { *ost << "  calc_hamiltonian" << std::endl; }
    if (am_i_mpi_rank0) { *ost << "  cf. dimension of subspace (up, ik=0) = " << phi[0][0].size() << std::endl; }

    // initialize Hphi (0.0)
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            const int num_G_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.num_G_at_k_scf()[ik] : plane_wave_basis.num_G_at_k_band()[ik];
 
            for (int iband=0; iband<H1phi[ispin][ik].size(); iband++)
            {
                for (int ispinor=0; ispinor<num_spinor; ispinor++)
                {
                    H1phi[ispin][ik][iband][ispinor] = Eigen::VectorXcd::Zero(num_G_at_k);
                    H2phi[ispin][ik][iband][ispinor] = Eigen::VectorXcd::Zero(num_G_at_k);
                    if (method.calc_method()=="TC" || method.calc_method()=="BITC")
                    {
                        H3phi[ispin][ik][iband][ispinor] = Eigen::VectorXcd::Zero(num_G_at_k);
                    }
                }
            }
        }
    }

    kinetic(parallelization,
            method,
            crystal_structure,
            spin, kpoints,
            plane_wave_basis,
            phi, H1phi, ost);
    if (!potentials.is_heg()) // pseudopotentials are not used for the homogeneous-electron-gas mode
    {
        pseudo(parallelization,
               method,
               crystal_structure,
               potentials, spin, kpoints,
               plane_wave_basis,
               phi, H1phi, ost);
    }
    if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "calculating the one-body terms"); }

    if (method.calc_method()=="HF")
    {
        // Hartree
        hf2h(parallelization,
             method,
             crystal_structure,
             potentials, spin, kpoints,
             plane_wave_basis, bloch_states,
             phi, H2phi, ost);
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "calculating the Hartree term"); }
     
        // exchange
        hf2x(parallelization,
             method,
             crystal_structure,
             potentials, spin, kpoints,
             plane_wave_basis, bloch_states,
             phi, H2phi, ost);
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "calculating the exchange term"); }
    }
    else if (method.calc_method()=="TC" || method.calc_method()=="BITC")
    {
        // two-body terms
        // Hartree
        tc2h(parallelization,
             method,
             crystal_structure,
             potentials, spin, kpoints,
             plane_wave_basis, bloch_states,
             phi, H2phi, ost);
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "calculating the tc-2body Hartree term"); }

        // exchange
        tc2x(parallelization,
             method,
             crystal_structure,
             potentials, spin, kpoints,
             plane_wave_basis, bloch_states,
             phi, H2phi, ost);
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "calculating the tc-2body exchange term"); }

        // three-body terms

        // tc3a1: -0.5*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |j,q1,q2>
        // tc3a2: 0.5*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |j,q2,q1>
        // tc3a3: 0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q1,j,q2>
        //   [equivalent then including: 3a6] 0.5*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q2,q1,j>
        // tc3a4: -0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q1,q2,j>
        //   [equivalent then including: 3a5] -0.5*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |q2,j,q1>

        // tc3b1: -0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |j,q1,q2>
        // tc3b2: 0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |j,q2,q1>
        // tc3b3: 0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |q1,j,q2>
        // tc3b4: -0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |q1,q2,j>
        // tc3b5: -0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |q2,j,q1>
        // tc3b6: 0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |q2,q1,j>
        // [all tc3b* (\nabla_2 u_21 \nabla_2 u_23) are equivalent to tc3c* (\nabla_3 u_31 \nabla_3 u_32)
        //  by considering x2<->x3 and q1<->q2. Thus, multiplied by 2.]

        tc3a1(parallelization,
              method,
              crystal_structure,
              potentials, spin, kpoints,
              plane_wave_basis, bloch_states,
              phi, H3phi, ost);
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "calculating the tc-3body term (3a1)"); }

        tc3b1(parallelization,
              method,
              crystal_structure,
              potentials, spin, kpoints,
              plane_wave_basis, bloch_states,
              phi, H3phi, ost);
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "calculating the tc-3body term (3b1)"); }

        tc3a2a4b2b5(parallelization,
                    method,
                    crystal_structure,
                    symmetry,
                    potentials, spin, kpoints,
                    plane_wave_basis, bloch_states,
                    phi, H3phi, ost);
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "calculating the tc-3body term (3a2a4b2b5)"); }

        tc3a3b3b4b6(parallelization,
                    method,
                    crystal_structure,
                    potentials, spin, kpoints,
                    plane_wave_basis, bloch_states,
                    phi, H3phi, ost);
        if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "calculating the tc-3body term (3a3b3b4b6)"); }
    }

    // reduce Hphi (not allreduce!)
    // It is better to allreduce <phi|H|phi> (not H|phi>!) from the viewpoint of the communication cost,
    // but we do not do so for keeping readability of our code...
    Eigen::VectorXcd Hphi_local;
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            for (int iband=0; iband<H1phi[ispin][ik].size(); iband++)
            {
                for (int ispinor=0; ispinor<num_spinor; ispinor++)
                {
                    Hphi_local = H1phi[ispin][ik][iband][ispinor];
                    MPI_Reduce(Hphi_local.data(),
                               H1phi[ispin][ik][iband][ispinor].data(),
                               H1phi[ispin][ik][iband][ispinor].size(),
                               MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

                    Hphi_local = H2phi[ispin][ik][iband][ispinor];
                    MPI_Reduce(Hphi_local.data(),
                               H2phi[ispin][ik][iband][ispinor].data(),
                               H2phi[ispin][ik][iband][ispinor].size(),
                               MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
                    
                    if (method.calc_method()=="TC" || method.calc_method()=="BITC")
                    {
                        Hphi_local = H3phi[ispin][ik][iband][ispinor];
                        MPI_Reduce(Hphi_local.data(),
                                   H3phi[ispin][ik][iband][ispinor].data(),
                                   H3phi[ispin][ik][iband][ispinor].size(),
                                   MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
                    }
                }
            }
        }
    }
    if (am_i_mpi_rank0) { my_clock.print_time_from_save(ost, "Reduce() Hamiltonian*phi"); }
}
