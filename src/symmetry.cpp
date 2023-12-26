// [class Symmetry]
// symmetry operations

#include "include/header.hpp"

void Symmetry::set(const Spin &spin, const CrystalStructure &crystal_structure,
                   const int num_symmetry_operations,
                   const std::vector<Eigen::Matrix3i> &rotation,
                   const std::vector<Eigen::Vector3d> &translation,
                   const std::string &no_t_rev, const std::string &noinv)
{
    assert(num_symmetry_operations >= 0);
    assert(rotation.size() == num_symmetry_operations);
    assert(translation.size() == num_symmetry_operations);

    num_symmetry_operations_ = num_symmetry_operations;
    rotation_ = rotation;
    translation_ = translation;

    // set atom_index_after_sym_
    const int num_atoms = crystal_structure.num_atoms();
    atom_index_after_sym_.resize(num_symmetry_operations);
    for (int isym=0; isym<num_symmetry_operations_; isym++)
    {
        atom_index_after_sym_[isym].resize(num_atoms);
    }

    for (int isym=0; isym<num_symmetry_operations_; isym++)
    {
        for (int iatom=0; iatom<num_atoms; iatom++)
        {
            // perform a symmetry operation (in crystal coord.)
            Eigen::Vector3d atomic_position_after_sym =
                rotation_[isym].cast<double>() * 
                (crystal_structure.atomic_position_crystal()[iatom] + translation_[isym]);

            // search jatom_match
            int jatom_match = -1;
            for (int jatom=0; jatom<num_atoms; jatom++)
            {
                Eigen::Vector3d diff =
                    atomic_position_after_sym - crystal_structure.atomic_position_crystal()[jatom];
                for (int idim=0; idim<3; idim++)
                {
                    while (diff(idim) < -0.5) { diff(idim) += 1.0; }
                    while (diff(idim) > 0.5) { diff(idim) -= 1.0; }
                }
                if (diff.squaredNorm() < 1e-6) 
                {
                    jatom_match = jatom;
                    break;
                }
            } // jatom
            if (jatom_match == -1) 
            {
                error_messages::stop("A symmetrically equivalent atom was not found for some symmetry operation.");
            }
            if (crystal_structure.index_of_atoms()[iatom] != crystal_structure.index_of_atoms()[jatom_match])
            {
                error_messages::stop("Something is wrong in symmetry: some atom moves to an inequivalent site.");
            }

            atom_index_after_sym_[isym][iatom] = jatom_match;
        } // iatom
    } // isym

    // set order_: n such that S^n = identity for each symmetry operation S
    order_.resize(num_symmetry_operations);
    Eigen::Matrix3d matrix_tmp, matrix_after;
    Eigen::Vector3d vector_tmp, vector_after;
    for (int isym=0; isym<num_symmetry_operations_; isym++)
    {
        // initial
        matrix_tmp = Eigen::Matrix3d::Identity();
        vector_tmp = Eigen::Vector3d::Zero();

        int iorder;
        for (iorder=0; iorder<num_order_max_; iorder++) 
        {
            // [Before sym. op.] (mat_ij x_j + vec_i) * e_i where i,j-sum is taken
            // [After sym. op.: r -> S(r+t)] (S_ik mat_kj x_j + S_ik (vec_k + t_k)) * e_i
            matrix_after = rotation_[isym].cast<double>() * matrix_tmp;
            vector_after = rotation_[isym].cast<double>() * (vector_tmp + translation_[isym]);
            
            // check matrix_after = id. && vector_after = (0,0,0) (mod. Z)
            bool is_identity = true;
            if ((matrix_after - Eigen::Matrix3d::Identity()).squaredNorm() > 1e-6) { is_identity = false; }
            for (int idim=0; idim<3; idim++)
            {
                while (vector_after(idim) < -0.5) { vector_after(idim) += 1.0; }
                while (vector_after(idim) > 0.5) { vector_after(idim) -= 1.0; }
            }
            if (vector_after.squaredNorm() > 1e-6) { is_identity = false; }
            if (is_identity) { break; }

            // update
            matrix_tmp = matrix_after;
            vector_tmp = vector_after;
        } // iorder
        if (iorder==num_order_max_) { error_messages::stop("cannot set order in symmetry. Try increasing num_order_max_ in symmetry.hpp"); }
        order_[isym] = iorder + 1;
    } // isym

    // use time-reversal symmetry to make symmetrically-equivalent k-points
    if (!(no_t_rev=="true" && noinv=="true") && (spin.is_spinor() || spin.num_independent_spins()==2) )
    {
        error_messages::stop("Please specify noinv = .true. && no_t_rev = .true. in QE unless nspin=1...");
    }
    is_time_reversal_symmetric_ = (noinv=="true") ? false : true; // Note!
}

void Symmetry::bcast(const bool am_i_mpi_rank0)
{
    MPI_Bcast(&num_symmetry_operations_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!am_i_mpi_rank0)
    {
        rotation_.resize(num_symmetry_operations_);
        translation_.resize(num_symmetry_operations_);
        atom_index_after_sym_.resize(num_symmetry_operations_);
        order_.resize(num_symmetry_operations_);
    }
    for (int isym=0; isym<num_symmetry_operations_; isym++)
    {
        MPI_Bcast(rotation_[isym].data(), 9, MPI_INT, 0, MPI_COMM_WORLD);        
        MPI_Bcast(translation_[isym].data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    int num_atoms = atom_index_after_sym_[0].size();
    MPI_Bcast(&num_atoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!am_i_mpi_rank0)
    {
        for (int isym=0; isym<num_symmetry_operations_; isym++)
        {
            atom_index_after_sym_[isym].resize(num_atoms);
        }
    }
    for (int isym=0; isym<num_symmetry_operations_; isym++)
    {
        MPI_Bcast(&atom_index_after_sym_[isym][0], num_atoms, MPI_INT, 0, MPI_COMM_WORLD);        
    }

    MPI_Bcast(&order_[0], num_symmetry_operations_, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&is_time_reversal_symmetric_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
}
