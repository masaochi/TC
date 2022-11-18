// [namspace read_qe]
// read Quantum-Espresso (QE) files
// declared in include/read_qe.hpp,
// defined in  read_qe.cpp, read_qe_xml.cpp, read_qe_wfc.f90

#include "include/header.hpp"

namespace
{

std::string get_element(const boost::property_tree::ptree &pt, const std::string &comment,
                        const std::string &element_name, std::ostream *ost)
{
    if (auto tmp = pt.get_optional<std::string>(element_name))
    {
        *ost << "  " << comment << element_name << " = " << tmp << std::endl;
        return tmp.get();
    }
    else
    {
        error_messages::not_found(element_name,"the QE xml file");
    }
}

// Keyword "element_name" is read from the QE xml file
std::string read_qe_xml_element(const boost::property_tree::ptree &pt, const std::string &element_name,
                                std::ostream *ost)
{
    return get_element(pt, "", element_name, ost);
}

// Keyword "element_name" is read in "tree_name" (child tree) from the QE xml file
std::string read_qe_xml_element_child(const boost::property_tree::ptree &pt, const std::string &tree_name,
                                      const std::string &element_name, std::ostream *ost)
{
    return get_element(pt, " " + tree_name + ".", element_name, ost); // tree_name is just used as comment
}

// transform string into std::vector<double>
void string_to_dvect(const std::string &str, std::vector<double> &dvect_tmp)
{
    std::string str_sub;
    std::stringstream iss(str);
    
    int idata = 0;
    while (iss >> str_sub)
    {
        dvect_tmp[idata] = boost::lexical_cast<double>(str_sub);
        idata++;
    }
    assert(idata == dvect_tmp.size());
}

} // namespace

// Read the QE xml file
void read_qe::read_qe_xml(FileNames &file_names,
                          const Method &method,
                          CrystalStructure &crystal_structure,
                          Symmetry &symmetry,
                          Potentials &potentials,
                          Spin &spin,
                          Kpoints &kpoints, 
                          PlaneWaveBasis &plane_wave_basis, 
                          BlochStates &bloch_states,
                          TotalEnergy &total_energy,
                          std::ostream *ost)
{
    // open the xml file
    const std::string qe_xml = file_names.qe_xml();
    *ost << " Read QE information from " << qe_xml << std::endl;
    boost::property_tree::ptree pt;
    boost::property_tree::xml_parser::read_xml(qe_xml, pt);
    
    std::vector<double> dvect_tmp1, dvect_tmp2;
    
    // check consistency of calc_mode between QE and TC++
    method.check_consistency_of_calc_mode(read_qe_xml_element(pt, "qes:espresso.input.control_variables.calculation", ost));
    
    // lattice and reciprocal vectors
    dvect_tmp1.resize(3);
    Eigen::Matrix3d lattice_vectors, reciprocal_vectors, reciprocal_vectors_normalized;

    string_to_dvect(read_qe_xml_element(pt, "qes:espresso.output.atomic_structure.cell.a1", ost), dvect_tmp1);
    for (int idim=0; idim<3; idim++) { lattice_vectors(0, idim) = dvect_tmp1[idim]; }
    string_to_dvect(read_qe_xml_element(pt, "qes:espresso.output.atomic_structure.cell.a2", ost), dvect_tmp1);
    for (int idim=0; idim<3; idim++) { lattice_vectors(1, idim) = dvect_tmp1[idim]; }
    string_to_dvect(read_qe_xml_element(pt, "qes:espresso.output.atomic_structure.cell.a3", ost), dvect_tmp1);
    for (int idim=0; idim<3; idim++) { lattice_vectors(2, idim) = dvect_tmp1[idim]; }

    string_to_dvect(read_qe_xml_element(pt, "qes:espresso.output.basis_set.reciprocal_lattice.b1", ost), dvect_tmp1);
    for (int idim=0; idim<3; idim++) { reciprocal_vectors_normalized(0, idim) = dvect_tmp1[idim]; }
    string_to_dvect(read_qe_xml_element(pt, "qes:espresso.output.basis_set.reciprocal_lattice.b2", ost), dvect_tmp1);
    for (int idim=0; idim<3; idim++) { reciprocal_vectors_normalized(1, idim) = dvect_tmp1[idim]; }
    string_to_dvect(read_qe_xml_element(pt, "qes:espresso.output.basis_set.reciprocal_lattice.b3", ost), dvect_tmp1);
    for (int idim=0; idim<3; idim++) { reciprocal_vectors_normalized(2, idim) = dvect_tmp1[idim]; }

    // alat: lattice parameter that is used as the normalized constant in QE
    const double alat = boost::lexical_cast<double>(read_qe_xml_element(pt, "qes:espresso.output.atomic_structure.<xmlattr>.alat", ost));
    reciprocal_vectors = (2*PI/alat)*reciprocal_vectors_normalized;

    // atom informations
    double num_atomic_species = boost::lexical_cast<int>(read_qe_xml_element(pt, "qes:espresso.output.atomic_species.<xmlattr>.ntyp", ost));
    std::vector<std::string> name_atomic_species(num_atomic_species);
    std::vector<std::string> pseudo_file_atomic_species(num_atomic_species);
    int ispecies = 0;
    BOOST_FOREACH (auto child, pt.get_child("qes:espresso.output.atomic_species")) 
    {
        if (child.first=="species")
        {
            *ost << "  Species: " << ispecies << std::endl;

            name_atomic_species[ispecies] 
                = read_qe_xml_element_child(child.second, "qes:espresso.output.atomic_species.species", "<xmlattr>.name", ost);
            pseudo_file_atomic_species[ispecies] 
                = read_qe_xml_element_child(child.second, "qes:espresso.output.atomic_species.species", "pseudo_file", ost);
            ispecies++;
        }
    }
    assert(ispecies == num_atomic_species);
    file_names.set_upf_file_names(pseudo_file_atomic_species);

    const int num_atoms = boost::lexical_cast<int>(read_qe_xml_element(pt, "qes:espresso.output.atomic_structure.<xmlattr>.nat", ost));
    std::vector<int> index_of_atoms(num_atoms);
    std::vector<Eigen::Vector3d> atomic_position_cartesian(num_atoms);

    dvect_tmp1.resize(3);
    int iatom = 0;
    BOOST_FOREACH (auto child, pt.get_child("qes:espresso.output.atomic_structure.atomic_positions")) 
    {
        if (child.first=="atom")
        {
            *ost << "  Atom: " << iatom << std::endl;
            std::string name 
                = read_qe_xml_element_child(child.second, "qes:espresso.output.atomic_structure.atomic_positions.atom", "<xmlattr>.name", ost);
            auto itr = std::find(name_atomic_species.begin(), name_atomic_species.end(), name);
            if (itr == name_atomic_species.end())
            {
                error_messages::stop("Name of atomic species is not found.");
            }
            else
            {
                index_of_atoms[iatom] = std::distance(name_atomic_species.begin(), itr); 
            }

            string_to_dvect(child.second.data(), dvect_tmp1);
            for (int idim=0; idim<3; idim++) { atomic_position_cartesian[iatom](idim) = dvect_tmp1[idim]; }
            *ost << "   qes.espresso.output.atomic_structure.atomic_positions.atom = " << atomic_position_cartesian[iatom].transpose() << std::endl;

            iatom++;
        }
    }
    assert(iatom == num_atoms);
    crystal_structure.set(lattice_vectors, reciprocal_vectors, num_atomic_species, index_of_atoms, atomic_position_cartesian);

    // spin    
    spin.set(read_qe_xml_element(pt, "qes:espresso.output.magnetization.noncolin", ost),
             read_qe_xml_element(pt, "qes:espresso.output.magnetization.lsda", ost));
    
    // symmetry operations
    const int num_symmetry_operations = boost::lexical_cast<int>(read_qe_xml_element(pt,"qes:espresso.output.symmetries.nsym",ost));
    std::vector<Eigen::Matrix3i> rotation(num_symmetry_operations);
    std::vector<Eigen::Vector3d> translation(num_symmetry_operations);
    
    dvect_tmp1.resize(9);
    dvect_tmp2.resize(3);
    int isym = 0;
    BOOST_FOREACH (auto child, pt.get_child("qes:espresso.output.symmetries")) 
    {
        if (child.first=="symmetry")
        {
            // cannot use this symmetry (cannot find "equivalent_atoms")
            if (!child.second.get_optional<std::string>("equivalent_atoms")) { continue; } 
            
            string_to_dvect(read_qe_xml_element_child(child.second, "qes:espresso.output.symmetries.symmetry",
                                                      "rotation", ost), dvect_tmp1);
            for (int idim=0; idim<9; idim++) { rotation[isym](idim/3, idim%3) = std::round(dvect_tmp1[idim]); }

            string_to_dvect(read_qe_xml_element_child(child.second, "qes:espresso.output.symmetries.symmetry",
                                                      "fractional_translation", ost), dvect_tmp2);
            for (int idim=0; idim<3; idim++) { translation[isym](idim) = dvect_tmp2[idim]; }                
            isym++;
        }
    }
    assert(isym == num_symmetry_operations);
    symmetry.set(spin, num_symmetry_operations, rotation, translation,
                 read_qe_xml_element(pt, "qes:espresso.input.symmetry_flags.no_t_rev", ost),
                 read_qe_xml_element(pt, "qes:espresso.input.symmetry_flags.noinv", ost));
    
    // read nbnd
    std::vector<int> nbnd(spin.num_independent_spins());
    if (spin.num_independent_spins()==2)
    {
        nbnd[0] = boost::lexical_cast<int>(read_qe_xml_element(pt, "qes:espresso.output.band_structure.nbnd_up", ost));
        nbnd[1] = boost::lexical_cast<int>(read_qe_xml_element(pt, "qes:espresso.output.band_structure.nbnd_dw", ost));
    }
    else
    {
        nbnd[0] = boost::lexical_cast<int>(read_qe_xml_element(pt, "qes:espresso.output.band_structure.nbnd", ost));
    }
    
    // read nks (no. of irreducible k-points)
    const int num_irreducible_kpoints = boost::lexical_cast<int>(read_qe_xml_element(pt, "qes:espresso.output.band_structure.nks", ost));
    
    // read k-point information, eigenvalues, and npw
    Eigen::Vector3d kvectors_irred_cart; // cartesian coordinate
    std::vector<Eigen::Vector3d> kvectors_irred_crys(num_irreducible_kpoints); // crystal coordinate
    std::vector<double> kweight(num_irreducible_kpoints);
    std::vector<int> num_G_at_k(num_irreducible_kpoints);
    std::vector<std::vector<double> > eigenvalues(num_irreducible_kpoints);
    
    dvect_tmp1.resize(3);
    int ik = 0;
    BOOST_FOREACH (auto child, pt.get_child("qes:espresso.output.band_structure")) 
    {
        if (child.first=="ks_energies")
        {
            *ost << "  Irreducible k-point: " << ik << std::endl;

            string_to_dvect(read_qe_xml_element_child(child.second, "qes:espresso.output.band_structure.ks_energies",
                                                      "k_point", ost), dvect_tmp1);
            for (int idim=0; idim<3; idim++) { kvectors_irred_cart(idim) = dvect_tmp1[idim]; }

            // cartesian (in 2pi/alat) => crystal coordinate. kvec_irred_cart = (b1 b2 b3)*kvec_irred_crys.
            kvectors_irred_crys[ik] = (reciprocal_vectors_normalized.transpose()).colPivHouseholderQr().solve(kvectors_irred_cart);
            
            kweight[ik] = boost::lexical_cast<double>(read_qe_xml_element_child(child.second, "qes:espresso.output.band_structure.ks_energies",
                                                                                "k_point.<xmlattr>.weight", ost));            
            num_G_at_k[ik] = boost::lexical_cast<int>(read_qe_xml_element_child(child.second, "qes:espresso.output.band_structure.ks_energies",
                                                                                "npw", ost));

            int nbands_tmp = boost::lexical_cast<double>(read_qe_xml_element_child(child.second, "qes::espresso.output.band_structure_ks_energies",
                                                                                   "eigenvalues.<xmlattr>.size",ost));
            eigenvalues[ik].resize(nbands_tmp); // nbands_tmp = nbnd[up] + nbnd[dn] for nspin=2, otherwise nbnd[0]
            string_to_dvect(read_qe_xml_element_child(child.second, "qes:espresso.output.band_structure.ks_energies",
                                                      "eigenvalues", ost), eigenvalues[ik]);
            ik++;
        }
    }
    assert(ik == num_irreducible_kpoints);

    // set these varaibles & make equivalent k-points using symmetry (in the SCF mode)
    kpoints.set_qe_input(spin, symmetry, num_irreducible_kpoints, kweight, kvectors_irred_crys, method.calc_mode(), ost);

    bloch_states.set_qe_xml(nbnd, eigenvalues, method.calc_mode(),
                            boost::lexical_cast<double>(read_qe_xml_element(pt, "qes:espresso.output.band_structure.nelec", ost)));
    plane_wave_basis.resize_G_at_k(spin.num_independent_spins(), num_irreducible_kpoints, method.calc_mode());
    plane_wave_basis.set_num_G_at_k(num_G_at_k, method.calc_mode());
    
    // set FFT grid
    const int nr1 = boost::lexical_cast<int>(read_qe_xml_element(pt, "qes:espresso.output.basis_set.fft_grid.<xmlattr>.nr1", ost));
    const int nr2 = boost::lexical_cast<int>(read_qe_xml_element(pt, "qes:espresso.output.basis_set.fft_grid.<xmlattr>.nr2", ost));
    const int nr3 = boost::lexical_cast<int>(read_qe_xml_element(pt, "qes:espresso.output.basis_set.fft_grid.<xmlattr>.nr3", ost));
    plane_wave_basis.setup_FFT(nr1, nr2, nr3);
    
    // Ewald energy
    total_energy.set_ewald_energy(potentials.is_heg(),
                                  boost::lexical_cast<double>
                                  (read_qe_xml_element(pt, "qes:espresso.output.total_energy.ewald", ost)));
    
    // Pseudopotential information
    if (read_qe_xml_element(pt, "qes:espresso.output.algorithmic_info.uspp", ost) == "true" ||
        read_qe_xml_element(pt, "qes:espresso.output.algorithmic_info.paw", ost) == "true")
    {
        error_messages::stop("Only norm-conserving pseudopotential is supported in TC++...");
    }

    // Jastrow parameters normalization using the cell volume and the num. of electrons
    if (method.calc_method()=="TC" || method.calc_method()=="BITC") 
    {
        potentials.jastrow.normalize_A_long(crystal_structure, bloch_states, spin); 
    }

    *ost << std::endl;
}
