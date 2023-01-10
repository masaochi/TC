// [namespace read_upf]
// Read upf (pseudopotential) files and MPI_Bcast

#include "include/header.hpp"

namespace
{

// e.g. find "PP_HEADER"    
bool find_tag(std::ifstream &ifs, const std::string tag_name)
{
    std::string sline;
    while (std::getline(ifs, sline))
    {
        if (sline.find(tag_name) != std::string::npos) // found
        {
            return true;
        }
    }
    return false;
}

// e.g. <PP_BETA.1 \n type="real" \n .... > then, getline() until reaching ">"
bool find_tag_until_end_of_tag(std::ifstream &ifs, const std::string tag_name)
{
    std::string sline;
    while (std::getline(ifs, sline))
    {
        if (sline.find(tag_name) != std::string::npos) // found
        {
            while (sline.find(">") == std::string::npos) // not found
            {
                std::getline(ifs, sline); // getline() until reaching ">"
            }
            return true;
        }
    }
    return false;
}

// e.g., get 6 from 'number_of_proj="6"'
std::string get_keyword(const std::string &sline)
{
    std::string s1 = sline.substr(sline.find("\"") + 1); // s1 = '6"'
    // s1.pop_back(); // remove "
    s1.erase(s1.end() - 1); // remove "
    return s1;
}

// almost std::getline() but return the first keyword
// e.g., "aaa bbb ccc" => return "aaa"
std::string getline_and_return_first(std::ifstream &ifs)
{
    std::string sline;
    std::getline(ifs, sline);
    std::stringstream ss{sline};

    std::string s;
    ss >> s;
    return s;
}
    
std::string getline_and_return_second(std::ifstream &ifs)
{
    std::string sline;
    std::getline(ifs, sline);
    std::stringstream ss{sline};

    std::string s1, s2;
    ss >> s1 >> s2;
    return s2;
}

// read array data and set rvec (e.g. PP_MESH/PP_R)    
void get_array(std::ifstream &ifs, std::vector<double> &rvec, const std::string tag_name)
{
    std::string sline, s1;

    if (!find_tag_until_end_of_tag(ifs, tag_name))
    {
        s1 = "Error in reading " + tag_name;
        error_messages::stop(s1);
    }
    int imesh = 0;
    while (std::getline(ifs, sline))
    {
        if (sline.find("/" + tag_name) != std::string::npos) { break; } // end reading

        std::stringstream ss{sline};
        while (ss >> s1)
        {
            rvec[imesh] = std::stod(s1);
            imesh++;
        }
    }
    assert(imesh == rvec.size());
}

std::string get_element(const boost::property_tree::ptree &pt, const std::string &comment,
                        const std::string &element_name, std::ostream *ost)
{
    if (auto tmp = pt.get_optional<std::string>(element_name))
    {
        *ost << "   " << comment << element_name << " = " << tmp << std::endl;
        return tmp.get();
    }
    else
    {
        error_messages::not_found(element_name,"upf files");
    }
}

// Keyword "element_name" is read from the upf file
std::string read_upf_element(const boost::property_tree::ptree &pt, const std::string &element_name,
                             std::ostream *ost)
{
    return get_element(pt, "", element_name, ost);
}

} // namespace

void read_upf::read_upf(FileNames &file_names,
                        CrystalStructure &crystal_structure,
                        Potentials &potentials,
                        PlaneWaveBasis &plane_wave_basis,
                        const bool am_i_mpi_rank0,
                        std::ostream *ost)
{
    if (am_i_mpi_rank0) // open upf files
    {
        *ost << " Read_upf" << std::endl;

        const int num_atomic_species = file_names.upf().size();

        // variables for reading upf files
        std::vector<double> Z_valence(num_atomic_species);
        //std::vector<int> lmax(num_atomic_species); // maximum angular momentum

        std::vector<std::vector<double> > pseudo_rmesh(num_atomic_species); // radial mesh
        // \int f(r) dr = \sum_i c_i f(i)*weight_rmesh(i). c_i = coefficients in Simpson's rule
        std::vector<std::vector<double> > pseudo_weight_rmesh(num_atomic_species);
        std::vector<std::vector<double> > pseudo_local(num_atomic_species); // local potential
        //std::vector<std::vector<double> > pseudo_nlcc(num_atomic_species); // non-linera core correction

        // Kleinman-Bylander projectors (pseudo_beta), their angular momenta (pseudo_lbeta),
        // and their coefficients (pseudo_dij)
        std::vector<std::vector<std::vector<double> > > pseudo_beta(num_atomic_species);
        std::vector<std::vector<int> > pseudo_lbeta(num_atomic_species);
        std::vector<Eigen::MatrixXd> pseudo_dij(num_atomic_species);

        for (int iatomic_species=0; iatomic_species<num_atomic_species; iatomic_species++)
        {
            *ost << "  Species = " << iatomic_species << std::endl;
            std::ifstream ifs(file_names.upf()[iatomic_species]);
            if (ifs.fail()) { error_messages::cannot_open(file_names.upf()[iatomic_species]); }

            // read upf
            bool is_upf_ver2 = check_upf_ver2(ifs, ost); // UPF ver.2 or ver.1 (XML schema not supported)
            int mesh_size; // num. of radial mesh grid
            int num_projectors; // num. of Kleinman-Bylander projectors
            read_upf_header(ifs, file_names.upf()[iatomic_species],
                            is_upf_ver2, Z_valence[iatomic_species], 
                            mesh_size, num_projectors, ost);

            // set some arrays
            pseudo_rmesh[iatomic_species].resize(mesh_size);
            pseudo_weight_rmesh[iatomic_species].resize(mesh_size);
            pseudo_local[iatomic_species].resize(mesh_size);
            pseudo_beta[iatomic_species].resize(num_projectors);
            for (int iproj=0; iproj<num_projectors; iproj++) { pseudo_beta[iatomic_species][iproj].resize(mesh_size); }
            pseudo_lbeta[iatomic_species].resize(num_projectors);
            pseudo_dij[iatomic_species].resize(num_projectors, num_projectors);
            //pseudo_nlcc[iatomic_species].resize(mesh_size);

            // read upf
            read_upf_mesh(ifs, pseudo_rmesh[iatomic_species], pseudo_weight_rmesh[iatomic_species]);
            read_upf_local(ifs, pseudo_local[iatomic_species]);
            read_upf_nonlocal(ifs, is_upf_ver2, pseudo_beta[iatomic_species], 
                              pseudo_lbeta[iatomic_species], pseudo_dij[iatomic_species], ost);
            //if (is_nlcc_included[iatomic_species]) { read_upf_nlcc(ifs, pseudo_nlcc[iatomic_species]); }
        }
        if (am_i_mpi_rank0) { *ost << " set_upf" << std::endl; }
        potentials.set_upf(Z_valence, pseudo_rmesh, pseudo_weight_rmesh,
                           pseudo_local, pseudo_beta, pseudo_lbeta, pseudo_dij,
                           crystal_structure, plane_wave_basis);
        *ost << std::endl;
    }
    if (am_i_mpi_rank0) { *ost << " Bcast_upf" << std::endl; }
    potentials.bcast_upf(am_i_mpi_rank0);
    if (am_i_mpi_rank0) { *ost << std::endl; }
}
    
bool read_upf::check_upf_ver2(std::ifstream &ifs, std::ostream *ost)
{
    std::string sline;
    std::getline(ifs, sline);
    if (sline.find("qe_pp:pseudo") != std::string::npos)
    {
        error_messages::stop("UPF with XML schema is not supported...");
    }
    else if (sline.find("UPF") != std::string::npos)
    {
        *ost << "   UPF ver.2" << std::endl;
        return true; // UPF ver.2
    }
    else
    {
        *ost << "   UPF ver.1" << std::endl;
        return false; // UPF ver.1
    }
}

// read PP_HEADER
void read_upf::read_upf_header(std::ifstream &ifs, const std::string &file_name_upf,
                               const bool is_upf_ver2,
                               double &Z_valence, int &mesh_size, int &num_projectors,
                               std::ostream *ost)
{
    std::string sline, s1;
    std::stringstream ss;

    if (!find_tag(ifs, "PP_HEADER")) { error_messages::stop("Error: cannot find PP_HEADER"); }

    if (is_upf_ver2) // UPF ver.2
    {
        // read using boost
        boost::property_tree::ptree pt;
        boost::property_tree::xml_parser::read_xml(file_name_upf, pt);

        std::string pseudo_type = read_upf_element(pt, "UPF.PP_HEADER.<xmlattr>.pseudo_type", ost);
        std::string is_ultrasoft = read_upf_element(pt, "UPF.PP_HEADER.<xmlattr>.is_ultrasoft", ost);
        std::string is_paw = read_upf_element(pt, "UPF.PP_HEADER.<xmlattr>.is_paw", ost);
        std::string is_coulomb = read_upf_element(pt, "UPF.PP_HEADER.<xmlattr>.is_coulomb", ost);

        // pseudo type (should be norm-conserving pseudopot.)
        if (pseudo_type=="US" || pseudo_type=="PAW" ||
            is_ultrasoft=="T" || is_paw=="T" || is_coulomb=="T")
        {
            error_messages::stop("Ultrasoft, PAW, 1/r (is_coulomb=\"T\") are not supported.");
        }

        std::string core_correction = read_upf_element(pt, "UPF.PP_HEADER.<xmlattr>.core_correction", ost);
        // nonlinear core correction
        if (core_correction=="T")
        {
            error_messages::stop("Non-linear core correction is not supported...");
        }
        else if (core_correction=="F")
        {
            *ost << "   Non-linear core correction is not included in pseudopot." << std::endl;
        }
        
        std::string Z_valence_string = read_upf_element(pt, "UPF.PP_HEADER.<xmlattr>.z_valence", ost);
        Z_valence = std::stod(Z_valence_string);
        *ost << "   Z_valence = " << Z_valence << std::endl;

        std::string mesh_size_string = read_upf_element(pt, "UPF.PP_HEADER.<xmlattr>.mesh_size", ost);
        mesh_size = std::stoi(mesh_size_string);
        *ost << "   mesh_size = " << mesh_size << std::endl;

        std::string number_of_proj_string = read_upf_element(pt, "UPF.PP_HEADER.<xmlattr>.number_of_proj", ost);
        num_projectors = std::stoi(number_of_proj_string);
        *ost << "   num. of projectors = " << num_projectors << std::endl;
    }
    else // UPF ver.1
    {
        std::getline(ifs, sline); // 1st line not read
        std::getline(ifs, sline); // 2nd line not read

        s1 = getline_and_return_first(ifs); // 3rd line: type of pseudopot
        if (s1!="NC") { error_messages::stop("Ultrasoft, PAW, 1/r are not supported."); }

        s1 = getline_and_return_first(ifs); // 4th line: nonlinear core correction
        if (s1=="T")
        {
            error_messages::stop("Non-linear core correction is not supported...");
        }
        else
        {
            *ost << "   Non-linear core correction is not included in pseudopot." << std::endl;
        }

        std::getline(ifs, sline); // 5th line not read

        Z_valence = std::stod(getline_and_return_first(ifs)); // 6th line: Z_valence
        *ost << "   Z_valence = " << Z_valence << std::endl;

        std::getline(ifs, sline); // 7th line not read
        std::getline(ifs, sline); // 8th line not read

        std::getline(ifs, sline); // 9th line not read (lmax)
        //lmax = std::stoi(getline_and_return_first(ifs)); // 9th line: maximum angular momentum (lmax)
        //*ost << "   lmax: " << lmax << std::endl;

        mesh_size = std::stoi(getline_and_return_first(ifs)); // 10th line: number of points in the radial grid
        *ost << "   mesh_size = " << mesh_size << std::endl;

        num_projectors = std::stoi(getline_and_return_second(ifs)); // 11th line: (second argument) num. of projectors
        *ost << "   num. of projectors = " << num_projectors << std::endl;
    }
}

// read PP_MESH
void read_upf::read_upf_mesh(std::ifstream &ifs, std::vector<double> &pseudo_rmesh,
                             std::vector<double> &pseudo_weight_rmesh)
{
    ifs.clear();
    ifs.seekg(0, std::ios_base::beg);

    get_array(ifs, pseudo_rmesh, "PP_R");
    get_array(ifs, pseudo_weight_rmesh, "PP_RAB");
}

// read PP_LOCAL
void read_upf::read_upf_local(std::ifstream &ifs, std::vector<double> &pseudo_local)
{
    ifs.clear();
    ifs.seekg(0, std::ios_base::beg);

    get_array(ifs, pseudo_local, "PP_LOCAL");
}

// read PP_NONLOCAL    
void read_upf::read_upf_nonlocal(std::ifstream &ifs, const bool is_upf_ver2,
                                 std::vector<std::vector<double> > &pseudo_beta,
                                 std::vector<int> &pseudo_lbeta,
                                 Eigen::MatrixXd &pseudo_dij,
                                 std::ostream *ost)
{
    ifs.clear();
    ifs.seekg(0, std::ios_base::beg);

    std::string sline, s1;
    int num_projectors = pseudo_beta.size();

    if (is_upf_ver2) // UPF ver.2
    {
        bool is_beta_found;
        for (int iproj=0; iproj<num_projectors; iproj++)
        {
            while (std::getline(ifs, sline))
            {
                if (sline.find("PP_BETA") != std::string::npos) // found
                {
                    while (sline.find("angular_momentum") == std::string::npos) // not found
                    {
                        std::getline(ifs, sline);
                    }
                    s1 = sline.substr(sline.find("angular_momentum"));
                    s1 = s1.substr(s1.find("\"") + 1);
                    s1 = s1.erase(s1.find("\""));
                    pseudo_lbeta[iproj] = std::stoi(s1); // set pseudo_lbeta
                    *ost << "    angular momentum = " << pseudo_lbeta[iproj] << std::endl;

                    while (sline.find(">") == std::string::npos) // not found
                    {
                        std::getline(ifs, sline); // getline() until reaching ">"
                    }
                    is_beta_found = true;
                    break;
                }
                is_beta_found = false;
            }
            if (!is_beta_found) { error_messages::stop("Error: PP_BETA not found."); }

            int imesh = 0;
            while (std::getline(ifs, sline))
            {
                if (sline.find("/PP_BETA") != std::string::npos) { break; } // end reading
                
                std::stringstream ss{sline};
                while (ss >> s1)
                {
                    pseudo_beta[iproj][imesh] = std::stod(s1);
                    imesh++;
                }
            }
            assert(imesh == pseudo_beta[iproj].size());
        } // int iproj

        std::vector<double> tmp_dij(num_projectors * num_projectors);
        get_array(ifs, tmp_dij, "PP_DIJ");
        for (int iproj=0; iproj<num_projectors; iproj++)
        {
            for (int jproj=0; jproj<num_projectors; jproj++)
            {
                pseudo_dij(iproj, jproj) = tmp_dij[iproj*num_projectors + jproj];
            }
        }
    }
    else // UPF ver.1
    {
        for (int iproj=0; iproj<num_projectors; iproj++)
        {
            if (!find_tag(ifs, "PP_BETA")) { error_messages::stop("Error: PP_BETA not found."); }

            pseudo_lbeta[iproj] = std::stoi(getline_and_return_second(ifs)); // angular momentum
            *ost << "    angular momentum = " << pseudo_lbeta[iproj] << std::endl;
            std::getline(ifs, sline); // not used

            int imesh = 0;
            while (std::getline(ifs, sline))
            {
                if (sline.find("/PP_BETA") != std::string::npos) { break; } // end reading

                std::stringstream ss{sline};
                while (ss >> s1)
                {
                    pseudo_beta[iproj][imesh] = std::stod(s1);
                    imesh++;
                }
            }
            assert(imesh == pseudo_beta[iproj].size());
        }

        if (!find_tag(ifs, "PP_DIJ")) { error_messages::stop("Error: PP_DIJ not found"); }
        int num_elements = std::stoi(getline_and_return_first(ifs)); // num. of non-zero dij elements
        pseudo_dij = Eigen::MatrixXd::Zero(num_projectors, num_projectors);

        int icheck = 0;
        while (std::getline(ifs, sline))
        {
            if (sline.find("/PP_DIJ") != std::string::npos) { break; } // end reading

            std::stringstream ss{sline};
            int i1, i2;
            double dij;
            ss >> i1 >> i2 >> dij;
            if (i1 > num_projectors || i2 > num_projectors)
            {
                error_messages::stop("Error in reading PP_DIJ");
            }
            pseudo_dij(i1-1, i2-1) = dij;
            icheck++;
        }
        assert(icheck == num_elements);
    }
}

// read PP_NLCC (not used at present)
void read_upf::read_upf_nlcc(std::ifstream &ifs, std::vector<double> &pseudo_nlcc)
{
    ifs.clear();
    ifs.seekg(0, std::ios_base::beg);

    get_array(ifs, pseudo_nlcc, "PP_NLCC");
}
