// [namespace io_tc_files]

#include "include/header.hpp"

namespace
{

std::vector<std::string> get_str(const std::string sline)
{
    std::stringstream ss{sline};
    std::vector<std::string> str(3);
    for (int i=0; i<3; i++)
    {
        str[i] = "";
        ss >> str[i];
    }
    return str;
}

} // namespace

void io_tc_files::read_crystal_structure(const FileNames &file_names,
                                         CrystalStructure &crystal_structure,
                                         const bool am_i_mpi_rank0,
                                         std::ostream *ost)
{
    const int num_atoms = crystal_structure.num_atoms();

    Eigen::Matrix3d lattice_vectors;
    std::vector<Eigen::Vector3d> atomic_position_cartesian(num_atoms);

    if (am_i_mpi_rank0) 
    {
        *ost << " Read a crystal structure (" << file_names.tc_crystal_structure() << ")" << std::endl;
        std::ifstream ifs(file_names.tc_crystal_structure());
        if (ifs.fail()) { error_messages::cannot_open(file_names.tc_crystal_structure()); }

        std::string sline;

        // comment line
        if (!std::getline(ifs, sline)) { error_messages::stop("Crystal structure file is incomplete."); }

        // lattice vectors in Bohr
        for (int idim=0; idim<3; idim++) // a1, a2, a3
        {
            if (!std::getline(ifs, sline)) { error_messages::stop("Crystal structure file is incomplete."); }

            std::vector<std::string> ai = get_str(sline);
            for (int jdim=0; jdim<3; jdim++) // x, y, z
            {
                lattice_vectors(idim, jdim) = boost::lexical_cast<double>(ai[jdim]);
            }
        }
        if (!crystal_structure.do_lattice_vectors_coincide(lattice_vectors)) // At present, lattice vectors are fixed during structural opt...
        {
            error_messages::stop("Lattice vectors are not consistent between QE and the crystal structure file.");
        }

        // comment line
        if (!std::getline(ifs, sline)) { error_messages::stop("Crystal structure file is incomplete."); }

        // atomic positions in cartesian coordinate
        for (int iatom=0; iatom<num_atoms; iatom++)
        {
            if (!std::getline(ifs, sline)) { error_messages::stop("Crystal structure file is incomplete."); }

            std::vector<std::string> pos = get_str(sline);
            for (int jdim=0; jdim<3; jdim++) // x, y, z
            {
                atomic_position_cartesian[iatom](jdim) = boost::lexical_cast<double>(pos[jdim]);
            }
        }
        ifs.close();
    } // am_i_mpi_rank0

    // Bcast
    for (int iatom=0; iatom<num_atoms; iatom++)
    {
        MPI_Bcast(atomic_position_cartesian[iatom].data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // update atomic positions
    crystal_structure.update_atomic_positions_by_file(atomic_position_cartesian, am_i_mpi_rank0, ost);
}
