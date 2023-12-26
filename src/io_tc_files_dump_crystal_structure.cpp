// [namespace io_tc_files]

#include "include/header.hpp"

void io_tc_files::dump_crystal_structure(const FileNames &file_names,
                                         const CrystalStructure &crystal_structure,
                                         std::ostream *ost)
{
    const int num_atoms = crystal_structure.num_atoms();

    *ost << " Dump a crystal structure (" << file_names.tc_crystal_structure() << ")" << std::endl;
    std::ofstream ofs(file_names.tc_crystal_structure());
    if (ofs.fail()) { error_messages::cannot_open(file_names.tc_crystal_structure()); }
    ofs << std::setprecision(10);

    ofs << " Lattice vectors in Bohr: a_i for i-th line" << std::endl;

    // lattice vectors in Bohr
    const Eigen::Matrix3d &lattice_vectors = crystal_structure.lattice_vectors();
    for (int idim=0; idim<3; idim++) // a1, a2, a3
    {
        ofs << " " << lattice_vectors(idim, 0)
            << " " << lattice_vectors(idim, 1)
            << " " << lattice_vectors(idim, 2) << std::endl;
    }

    ofs << " Atomic coordinates in cartesian coordinate" << std::endl;

    // atomic positions in cartesian coordinate
    const std::vector<Eigen::Vector3d> &atomic_position_cartesian = crystal_structure.atomic_position_cartesian();
    for (int iatom=0; iatom<num_atoms; iatom++)
    {
        ofs << " " << atomic_position_cartesian[iatom](0)
            << " " << atomic_position_cartesian[iatom](1)
            << " " << atomic_position_cartesian[iatom](2) << std::endl;
    }

    ofs.close();
}
