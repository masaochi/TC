// [namespace io_tc_files]
// Read input.in, read/write tc_wfc.dat (wave functions) and tc_energy.dat (orbital energies)

// This file contains functions regarding "input.in"

#include "include/header.hpp"

namespace
{
    
bool read_bool(const std::string &keyword, const std::string &s2)
{
    if (s2=="true" || s2=="True" || s2=="TRUE" || s2=="T" || s2=="t"
        || s2==".true." || s2==".True." || s2==".TRUE.")
    {
        return true;
    }
    else if (s2=="false" || s2=="False" || s2=="FALSE" || s2=="F" || s2=="f"
             || s2==".false." || s2==".False." || s2==".FALSE.")
    {
        return false;
    }
    else
    {
        error_messages::inappropriate_argument(keyword, s2, "should be true or false");
    }
}

std::string bool_to_string(const bool b)
{
    if (b) 
    { 
        return "true";
    }
    else
    {
        return "false";
    }
}

void already_set_error(const std::string &s1, std::ostream *ost)
{
    *ost << "The keyword " << s1 << " is set twice in input.in!" << std::endl;
    error_messages::stop("Error in reading input.in (a keyword is specified twice). See output.out for details.");
}

} // namespace

void io_tc_files::read_input_in(FileNames &file_names,
                                Diagonalization &diagonalization,
                                Method &method, Potentials &potentials,
                                Kpoints &kpoints, BlochStates &bloch_states, 
                                const bool am_i_mpi_rank0, std::ostream *ost)
{
    if (am_i_mpi_rank0) // open input.in
    {
        *ost << " Read TC input (" << file_names.tc_input() << ")" << std::endl;
        std::ifstream ifs(file_names.tc_input());
        if (ifs.fail()) { error_messages::cannot_open(file_names.tc_input()); }

        // if not overwritten by reading input.in,
        // default values (set by the constructor) are used
        std::string smearing_mode = kpoints.smearing_mode();
        double smearing_width = kpoints.smearing_width(); 

        bool restarts = diagonalization.restarts();
        bool includes_div_correction = potentials.includes_div_correction();

        double energy_tolerance = diagonalization.energy_tolerance();
        double charge_tolerance = diagonalization.charge_tolerance();
        int max_num_iterations = diagonalization.max_num_iterations();
        bool mixes_density_matrix = diagonalization.mixes_density_matrix();
        double mixing_beta = diagonalization.mixing_beta();
        int num_refresh_david = diagonalization.num_refresh_david();
        int max_num_blocks_david = diagonalization.max_num_blocks_david();

        bool is_heg = potentials.is_heg();

        std::vector<std::vector<double> > A_long = potentials.jastrow.A_long();

        std::string sline, s1, s2;
        std::vector<bool> is_already_set(20, false); // 20 = total num. of the input keywords
        while (std::getline(ifs, sline))
        {
            std::stringstream ss{sline};
            s1 = s2 = ""; // initialization is necessary considering the case when sline=="" (a void line)
            ss >> s1 >> s2; // Note that, when sline=="", this operation is not performed.
            if (s1=="calc_method") 
            {
                if (is_already_set[0]) { already_set_error(s1, ost); }
                method.set_calc_method(s2);
                is_already_set[0] = true;
            }
            else if (s1=="calc_mode") 
            {
                if (is_already_set[1]) { already_set_error(s1, ost); }
                method.set_calc_mode(s2); 
                is_already_set[1] = true;
            }
            else if (s1=="qe_save_dir") 
            { 
                if (is_already_set[2]) { already_set_error(s1, ost); }
                file_names.set_qe_save_dir(s2); 
                is_already_set[2] = true;
            }
            else if (s1=="pseudo_dir") 
            { 
                if (is_already_set[3]) { already_set_error(s1, ost); }
                file_names.set_pseudo_dir(s2); 
                is_already_set[3] = true;
            }
            else if (s1=="smearing_mode") 
            {
                if (is_already_set[5]) { already_set_error(s1, ost); }
                smearing_mode = s2;
                is_already_set[5] = true;
            }
            else if (s1=="smearing_width") 
            {
                if (is_already_set[6]) { already_set_error(s1, ost); }
                smearing_width = boost::lexical_cast<double>(s2); 
                is_already_set[6] = true;
            }
            else if (s1=="restarts") 
            {
                if (is_already_set[7]) { already_set_error(s1, ost); }
                restarts = read_bool(s1,s2);
                is_already_set[7] = true;
            }
            else if (s1=="includes_div_correction") 
            {
                if (is_already_set[8]) { already_set_error(s1, ost); }
                includes_div_correction = read_bool(s1,s2); 
                is_already_set[8] = true;
            }
            else if (s1=="energy_tolerance") 
            { 
                if (is_already_set[9]) { already_set_error(s1, ost); }
                energy_tolerance = boost::lexical_cast<double>(s2); 
                is_already_set[9] = true;
            }
            else if (s1=="charge_tolerance") 
            {
                if (is_already_set[10]) { already_set_error(s1, ost); }
                charge_tolerance = boost::lexical_cast<double>(s2); 
                is_already_set[10] = true;
            }
            else if (s1=="max_num_iterations") 
            { 
                if (is_already_set[11]) { already_set_error(s1, ost); }
                max_num_iterations = boost::lexical_cast<int>(s2); 
                is_already_set[11] = true;
            }
            else if (s1=="mixing_beta") 
            { 
                if (is_already_set[12]) { already_set_error(s1, ost); }
                mixing_beta = boost::lexical_cast<double>(s2);
                is_already_set[12] = true;
            }
            else if (s1=="num_refresh_david") 
            { 
                if (is_already_set[13]) { already_set_error(s1, ost); }
                num_refresh_david = boost::lexical_cast<int>(s2); 
                is_already_set[13] = true;
            }
            else if (s1=="max_num_blocks_david") 
            { 
                if (is_already_set[14]) { already_set_error(s1, ost); }
                max_num_blocks_david = boost::lexical_cast<int>(s2); 
                is_already_set[14] = true;
            }
            else if (s1=="is_heg")
            {
                if (is_already_set[15]) { already_set_error(s1, ost); }
                is_heg = read_bool(s1,s2); 
                is_already_set[15] = true;
            }
            else if (s1=="A_up_up") 
            { 
                if (is_already_set[16]) { already_set_error(s1, ost); }
                A_long[0][0] = boost::lexical_cast<double>(s2);
                is_already_set[16] = true;
            }
            else if (s1=="A_up_dn") 
            { 
                if (is_already_set[17]) { already_set_error(s1, ost); }
                A_long[0][1] = boost::lexical_cast<double>(s2); // = A_up_dn
                A_long[1][0] = boost::lexical_cast<double>(s2); // = A_dn_up
                is_already_set[17] = true;
            }
            else if (s1=="A_dn_dn") 
            { 
                if (is_already_set[18]) { already_set_error(s1, ost); }
                A_long[1][1] = boost::lexical_cast<double>(s2);
                is_already_set[18] = true;
            }
            else if (s1=="mixes_density_matrix")
            {
                if (is_already_set[19]) { already_set_error(s1, ost); }
                mixes_density_matrix = read_bool(s1,s2); 
                is_already_set[19] = true;
            }
            else if (s1!="" && s1.c_str()[0]!='#' && s1!="num_bands_tc") // neither a void line nor a comment line
            {
                *ost << "Error in reading input.in. Your keyword " << s1 << " is not valid. " << std::endl;
                *ost << "(If you would like write to a comment, please note that a comment line should begin with #)" << std::endl;
                error_messages::stop("Error in reading input.in (invalid keyword). See output.out for details.");
            }
        }
        if (!is_already_set[0]) { error_messages::stop("The keyword calc_method should be specified in input.in"); }
        if (!is_already_set[1]) { error_messages::stop("The keyword calc_mode should be specified in input.in"); }
        if (!is_already_set[2]) { error_messages::stop("The keyword qe_save_dir should be specified in input.in"); }
        if (!is_already_set[3]) { error_messages::stop("The keyword pseudo_dir should be specified in input.in"); }
        if (!((is_already_set[16] && is_already_set[17] && is_already_set[18]) ||
              (!is_already_set[16] && !is_already_set[17] && !is_already_set[18])))
        {
            error_messages::stop("Jastrow A parameters (A_up_up, A_up_dn, A_dn_dn) are only partially given in input.in");
        }

        ifs.clear();
        ifs.seekg(0, std::ios_base::beg);
        while (std::getline(ifs, sline)) // read num_bands_tc (can be num_bands_scf/band, need calc_mode!)
        {
            std::stringstream ss{sline};
            s1 = s2 = ""; // initialization is necessary considering the case when sline=="" (a void line)
            ss >> s1 >> s2; // Note that, when sline=="", this operation is not performed.
            if (s1=="num_bands_tc") 
            {
                if (is_already_set[4]) { already_set_error(s1, ost); }
                bloch_states.set_num_bands_tc({boost::lexical_cast<int>(s2)}, method.calc_mode());
                is_already_set[4] = true;
            }
        }

        kpoints.set_tc_input(smearing_mode, smearing_width);
        diagonalization.set(restarts,
                            energy_tolerance, charge_tolerance,
                            max_num_iterations, 
                            mixes_density_matrix, mixing_beta,
                            num_refresh_david, max_num_blocks_david,
                            method.calc_mode(), is_heg);
        potentials.set_includes_div_correction(includes_div_correction);
        potentials.set_is_heg(is_heg);
        potentials.jastrow.set(A_long);

        // Show parameters set here (if not specified in the input file, a default value is shown)
        *ost << "  calc_method = " << method.calc_method() << std::endl;
        *ost << "  calc_mode = " << method.calc_mode() << std::endl;
        *ost << "  qe_save_dir = " << file_names.qe_save_dir() << std::endl;
        *ost << "  pseudo_dir = " << file_names.pseudo_dir() << std::endl;

        const std::vector<int> &num_bands_tc = method.calc_mode()=="SCF"
            ? bloch_states.num_bands_scf() : bloch_states.num_bands_band();
        if (num_bands_tc[0]>=0) { *ost << "  num_bands_tc = " << num_bands_tc[0] << std::endl; }

        *ost << "  smearing_mode = " << kpoints.smearing_mode() << std::endl;
        *ost << "  smearing_width = " << kpoints.smearing_width() << std::endl;

        *ost << "  restarts = " << bool_to_string(diagonalization.restarts()) << std::endl;
        *ost << "  includes_div_correction = " << bool_to_string(potentials.includes_div_correction()) << std::endl;
        *ost << "  energy_tolerance = " << diagonalization.energy_tolerance() << std::endl;
        *ost << "  charge_tolerance = " << diagonalization.charge_tolerance() << std::endl;
        *ost << "  max_num_iterations = " << diagonalization.max_num_iterations() << std::endl;
        *ost << "  mixes_density_matrix = " << bool_to_string(diagonalization.mixes_density_matrix()) << std::endl;
        *ost << "  mixing_beta = " << diagonalization.mixing_beta() << std::endl;
        *ost << "  num_refresh_david = " << diagonalization.num_refresh_david() << std::endl;
        *ost << "  max_num_blocks_david = " << diagonalization.max_num_blocks_david() << std::endl;
        *ost << "  is_heg = " << bool_to_string(potentials.is_heg()) << std::endl;
        if (method.calc_method()=="TC" || method.calc_method()=="BITC")
        {
            *ost << "  A_up_up (unit: sqrt(V/(4piN)) ) = " << potentials.jastrow.A_long()[0][0] << std::endl;
            *ost << "  A_up_dn (unit: sqrt(V/(4piN)) ) = " << potentials.jastrow.A_long()[0][1] << std::endl;
            *ost << "  A_dn_up (unit: sqrt(V/(4piN)) ) = " << potentials.jastrow.A_long()[1][0] << std::endl;
            *ost << "  A_dn_dn (unit: sqrt(V/(4piN)) ) = " << potentials.jastrow.A_long()[1][1] << std::endl;
        }
    
        *ost << std::endl;
    } // am_i_mpi_rank0

    bcast_input_in(file_names, diagonalization,
                   method, potentials, kpoints,
                   bloch_states, am_i_mpi_rank0, ost);
}

void io_tc_files::bcast_input_in(FileNames &file_names,
                                 Diagonalization &diagonalization,
                                 Method &method, Potentials &potentials,
                                 Kpoints &kpoints, BlochStates &bloch_states,
                                 const bool am_i_mpi_rank0, std::ostream *ost)
{
    if (am_i_mpi_rank0) { *ost << " Bcast TC input (" << file_names.tc_input() << ")" << std::endl; }
    
    diagonalization.bcast(am_i_mpi_rank0);
    method.bcast(am_i_mpi_rank0);
    kpoints.bcast_tc_input(am_i_mpi_rank0);
    // Note: bloch_states.num_bands_tc_ will be Bcast in read_qe().
    // because this variable is possibly initialized by num_band_qe_ in read_qe_xml().
    potentials.bcast_tc_input();
    // Jatrow parameters will be bcast later (after normalization using QE variables)

    if (am_i_mpi_rank0) { *ost << std::endl; }
}
