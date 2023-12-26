// [class Diagonalization]
// scf_loop, convergence check, control paramters, etc.

#include "include/header.hpp"

// Quasi-Newton (BFGS)
void Diagonalization::structural_optimization(const Method &method,
                                              CrystalStructure &crystal_structure,
                                              const Symmetry &symmetry,
                                              TotalEnergy &total_energy,
                                              const int &istep_ionic,
                                              const bool am_i_mpi_rank0,
                                              std::ostream *ost) const
{
    const int num_atoms = total_energy.force().size();
    const int ndegrees = 3 * num_atoms; // num. of degrees of freedom

    Eigen::VectorXd position(ndegrees); // atomic unit, Cartesian
    Eigen::VectorXd position_new(ndegrees); // atomic unit, Cartesian
    Eigen::VectorXd gradient(ndegrees); // atomic unit, Cartesian
    Eigen::MatrixXd Hessian_inv(ndegrees, ndegrees);

    if (am_i_mpi_rank0)
    {
        // set position & gradient
        for (int idegree=0; idegree<ndegrees; idegree++)
        {
            int iatom = idegree/3; // 0, 0, 0, 1, 1, 1, ...
            int idim = idegree%3; // 0, 1, 2, 0, 1, 2, ...
            position(idegree) = crystal_structure.atomic_position_cartesian()[iatom](idim);
            gradient(idegree) = total_energy.force()[iatom](idim).real(); // imaginary part neglected
        }
        
        // (1) construct an approximate Hessian (by BFGS)
        if (istep_ionic==0) // first stesp
        {
            // line search gives "(typical step size) * normalized gradient" as a trial displacement
            Hessian_inv = (size_ionic_displacement_ / gradient.norm())
                * Eigen::MatrixXd::Identity(ndegrees, ndegrees);
        }
        else
        {
            assert(crystal_structure.atomic_position_cartesian_old().size() == num_atoms);
            assert(total_energy.force_old().size() == num_atoms);

            Eigen::VectorXd position_difference(ndegrees); // position - position_old
            Eigen::VectorXd gradient_difference(ndegrees); // gradient - gradient_old
            for (int idegree=0; idegree<ndegrees; idegree++)
            {
                int iatom = idegree/3; // 0, 0, 0, 1, 1, 1, ...
                int idim = idegree%3; // 0, 1, 2, 0, 1, 2, ...
                position_difference(idegree) = 
                    position(idegree) - crystal_structure.atomic_position_cartesian_old()[iatom](idim);
                gradient_difference(idegree) = 
                    gradient(idegree) - total_energy.force_old()[iatom](idim).real(); // imaginary part neglected
            }

            const Eigen::MatrixXd &Hessian_inv_old = total_energy.Hessian_inv_old();
            assert(Hessian_inv_old.rows() == ndegrees && Hessian_inv_old.cols() == ndegrees);

            // BFGS formula
            double prod = (position_difference.transpose() * gradient_difference)(0,0);
            double coeff = 1.0 + (gradient_difference.transpose() * Hessian_inv_old * gradient_difference)(0,0) / prod;
            Hessian_inv = Hessian_inv_old +
                coeff / prod * position_difference * position_difference.transpose() -
                ((position_difference * gradient_difference.transpose() * Hessian_inv_old +
                  Hessian_inv_old * gradient_difference * position_difference.transpose()) / prod);
        }

        // (2) set a direction
        Eigen::VectorXd displacement = - Hessian_inv * gradient;

        // (3) symmetrize a displacement vector
        Eigen::Vector3d displacement_sub_before, displacement_sub_after;
        for (int isym=0; isym<symmetry.num_symmetry_operations(); isym++)
        {
            const int &order = symmetry.order()[isym]; // n such that g^n = identity for symmetry operation g
            const std::vector<int> &atom_index_after_sym = symmetry.atom_index_after_sym()[isym];

            std::vector<Eigen::VectorXd> displacement_sym(order);
            for (int iorder=0; iorder<order; iorder++) { displacement_sym[iorder].resize(ndegrees); }

            displacement_sym[0] = displacement;
            for (int iorder=0; iorder<order-1; iorder++)
            {
                for (int iatom=0; iatom<num_atoms; iatom++)
                {
                    for (int idim=0; idim<3; idim++) // ndegrees -> 3-dim vector
                    {
                        displacement_sub_before(idim) = displacement_sym[iorder](3*iatom + idim);
                    }

                    // cartesian -> crystal coord.: lattice_vectors.transpose().inverse()
                    // apply rotation (note! no translation for force)
                    // crystal -> cartesian coord.: lattice_vector.transpose()
                    displacement_sub_after = 
                        crystal_structure.lattice_vectors().transpose() *
                        symmetry.rotation()[isym].cast<double>() * 
                        crystal_structure.lattice_vectors().transpose().inverse() *
                        displacement_sub_before;
                    
                    for (int idim=0; idim<3; idim++) // 3-dim -> ndegrees vector
                    {
                        displacement_sym[iorder+1](3*atom_index_after_sym[iatom] + idim)
                            = displacement_sub_after(idim);
                    }
                } // iatom
            } // iorder

            // symmetrizes displacement
            displacement = Eigen::VectorXd::Zero(ndegrees);
            for (int iorder=0; iorder<order; iorder++)
            {
                displacement = displacement + displacement_sym[iorder] / static_cast<double>(order);
            }
        } // isym

        // (4) Line search 
        double coeff_step_size = 0.5; // temporary
        displacement *= coeff_step_size;

        // (5) set a new position
        position_new = position + displacement;
        *ost << "   Symmetrized displacement vector in structural opt. (in Bohr, cartesian coord.)" << std::endl
             << displacement << std::endl;
    } // am_i_mpi_rank0

    // Bcast
    MPI_Bcast(Hessian_inv.data(), ndegrees*ndegrees, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(position_new.data(), ndegrees, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // save force_old and Hessian_inv_old
    total_energy.set_force_old(); // force_old_ = force_
    total_energy.set_Hessian_inv_old(Hessian_inv); // Hessian_inv_old_ = Hessian

    // update & print atomic_position_cartesian, set atomic_position_cartesian_old
    crystal_structure.update_atomic_positions(position_new, am_i_mpi_rank0, ost);
}
