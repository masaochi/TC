// [namespace io_qmc_files]
// Read/write parameters.casl, write pwfn.data, write jastrow.plt

#include "include/header.hpp"

namespace
{
    
std::vector<std::string> read_a_line_casl(std::ifstream &ifs)
{
    std::string sline;
    std::vector<std::string> keywords(2);

    std::getline(ifs, sline);
    boost::trim(sline);

    int colon_pos = sline.find(":"); // e.g., sline=="Type: RPA" -> colon_pos = 4
    if (colon_pos == std::string::npos)
    {
        error_messages::stop(": is not found. Invalid format of casl");
    }
    keywords[0] = sline.substr(0, colon_pos); // e.g., sline=="Type: RPA" -> keywords[0] = "Type"
    keywords[1] = colon_pos+2<=sline.size() ? sline.substr(colon_pos+2) : ""; // e.g., sline=="Type: RPA" -> keywords[1] = "RPA"
    
    if (keywords[1][0]=='[') // e.g., sline=="c_2: [ 0.0, fixed ]"
    {
        int comma_pos = keywords[1].find(",");
        if (comma_pos == std::string::npos) // "," was not found. e.g., sline=="Rules: [ 1-1=2-2 ]"
        {
            comma_pos = keywords[1].find("]");
            if (comma_pos == std::string::npos) // "]" was not found
            {
                error_messages::stop("Neither , nor ] are found. Invalid format of casl");
            }
            comma_pos -= 1;
        }

        // e.g., sline=="c_2: [ 0.0, ...." -> keywords[1] = "0.0"
        // e.g., sline=="Rules: [ 1-1=2-2 ]" -> keywords[1] = "1-1=2-2"
        keywords[1] = keywords[1].substr(2, comma_pos-2); 
    }
    return keywords;
}

} // namespace

void io_qmc_files::read_casl(const FileNames &file_names,
                             const Spin &spin,
                             const CrystalStructure &crystal_structure,
                             Jastrow &jastrow,
                             const PlaneWaveBasis &plane_wave_basis,
                             const bool am_i_mpi_rank0,
                             std::ostream *ost)
{
    bool is_RPA_read = false;
    bool is_POLY_read = false;
    bool cusp_poly; // whether the cusp condition is imposed on polynomial terms
    if (am_i_mpi_rank0)
    {
        *ost << " Read Jastrow parameters from " << file_names.tc_casl() << std::endl;
        std::ifstream ifs(file_names.tc_casl());
        if (ifs.fail()) 
        {
            *ost << "  The casl file was not found. Please make sure that the file should be placed in the current directory" << std::endl;
            *ost << "  if you would like to specify Jastrow parameters by " << file_names.tc_casl() << "." << std::endl;
            *ost << "  If you would like to specify Jastrow parameters by " << file_names.tc_input() << " or use the default Jastrow (A=1 RPA-type), do not mind this warning." << std::endl;
        }
        else
        {
            if (jastrow.is_A_given_in_input_in())
            {
                error_messages::stop("Do not specify Jastrow parameters in two files: input.in and parameters.casl. You can use one of them.");
            }

            jastrow.set_A_long({{0.0, 0.0}, {0.0, 0.0}}); // reset A parameters: default values in jastrow.hpp are non-zero.
            const int num_independent_spins = spin.num_independent_spins();
            std::string sline;

            // start reading!
            if (read_a_line_casl(ifs)[0]!="JASTROW") { error_messages::stop("The 1st line should be `JASTROW:' (read_casl)"); }
            if (read_a_line_casl(ifs)[0]!="Title") { error_messages::stop("The 2nd line should begin with `Title:' (read_casl)"); }
            if (read_a_line_casl(ifs)[0]!="TERM 1") { error_messages::stop("The 3rd line should be `TERM 1:' (read_casl)"); }
            read_each_term_casl(ifs, num_independent_spins, jastrow, is_RPA_read, is_POLY_read, cusp_poly);
            if (!ifs.eof()) { read_each_term_casl(ifs, num_independent_spins, jastrow, is_RPA_read, is_POLY_read, cusp_poly); } // "TERM 2" exists
        } // if (ifs.fail())
    } // if (am_i_mpi_rank0)
    MPI_Bcast(&is_RPA_read, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&is_POLY_read, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cusp_poly, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);

    // set & bcast & print RPA parameters
    if (is_RPA_read || is_POLY_read) // when casl was read: RPA parameters were overwritten
    {
        jastrow.bcast_A_long();
        jastrow.bcast_F_long();
    }
    else
    {
        if (am_i_mpi_rank0) { dump_casl(file_names, spin, jastrow, ost); }
    }
    if (am_i_mpi_rank0) { jastrow.print_RPA_parameters(ost); }
    jastrow.set_derived_RPA_parameters();

    // bcast & print polynomial parameters
    if (is_POLY_read) 
    {
        // consistency check for the cusp condition
        if (is_RPA_read && cusp_poly) { error_messages::stop("Cusp condition should not be imposed on polynomial terms when RPA terms are used (read_casl)"); }
//        if (!is_RPA_read && !cusp_poly) { error_messages::stop("Cusp condition should be imposed on polynomial terms when RPA terms are not used"); }
        jastrow.bcast_polynomial_parameters(am_i_mpi_rank0);
        if (am_i_mpi_rank0) { jastrow.print_polynomial_parameters(ost); }

        jastrow.set_derived_polynomial_parameters();
        if (is_RPA_read) { jastrow.set_derived_RPA_polynomial_parameters(crystal_structure, plane_wave_basis); }
        if (is_RPA_read && am_i_mpi_rank0) { jastrow.print_RPA_polynomial_parameters(ost); }
    }
}

void io_qmc_files::read_each_term_casl(std::ifstream &ifs, const int num_independent_spins, Jastrow &jastrow,
                                       bool &is_RPA_read, bool &is_POLY_read, bool &cusp_poly)
{
    if (read_a_line_casl(ifs)[0]!="Rank") { error_messages::stop("Rank is not given (read_casl)"); }
    if (read_a_line_casl(ifs)[0]!="Rules") { error_messages::stop("Rules is not given (read_casl)"); }
    if (read_a_line_casl(ifs)[0]!="e-e basis") { error_messages::stop("e-e-basis is not given (read_casl)"); }
    std::string jastrow_type = read_a_line_casl(ifs)[1];
    if (jastrow_type=="RPA")
    {
        is_RPA_read = read_RPA_casl(ifs, num_independent_spins, jastrow);
    }
    else if (jastrow_type=="natural power")
    {
        is_POLY_read = read_POLY_casl(ifs, num_independent_spins, jastrow, cusp_poly);
    }
    else
    {
        error_messages::stop("Type is unrecognized: should be `RPA' or `natural power' (read_casl)");
    }
}

bool io_qmc_files::read_RPA_casl(std::ifstream &ifs, const int num_independent_spins, Jastrow &jastrow)
{
    std::string sline;
    std::vector<std::string> keywords;
    std::vector<std::vector<double> > F_long(2);
    for (int is1=0; is1<2; is1++) { F_long[is1].resize(2); }

    keywords = read_a_line_casl(ifs);
    if (keywords[0]!="Order" || keywords[1]!="1") { error_messages::stop("Order should be given as `Order: 1` for RPA (read_casl)"); }

    // read F_long
    if (read_a_line_casl(ifs)[0]!="Parameters") { error_messages::stop("Parameters is not given for RPA (read_casl)"); }
    if (read_a_line_casl(ifs)[0]!="Channel 1-1") { error_messages::stop("Channel 1-1 is not given for RPA (read_casl)"); }
    keywords = read_a_line_casl(ifs);
    if (keywords[0]!="F") { error_messages::stop("F is not given for RPA (read_casl)"); }
    F_long[0][0] = boost::lexical_cast<double>(keywords[1]);
    std::getline(ifs, sline); // this line is not used

    if (read_a_line_casl(ifs)[0]!="Channel 1-2") { error_messages::stop("Channel 1-2 is not given for RPA (read_casl)"); }
    keywords = read_a_line_casl(ifs);
    if (keywords[0]!="F") { error_messages::stop("F is not given for RPA (read_casl)"); }
    F_long[0][1] = F_long[1][0] = boost::lexical_cast<double>(keywords[1]);
    std::getline(ifs, sline); // this line is not used

    if (num_independent_spins==1)
    {
        F_long[1][1] = F_long[0][0];
    }
    else
    {
        if (read_a_line_casl(ifs)[0]!="Channel 2-2") { error_messages::stop("Channel 2-2 is not given for RPA (read_casl)"); }
        keywords = read_a_line_casl(ifs);
        if (keywords[0]!="F") { error_messages::stop("F is not given for RPA (read_casl)"); }
        F_long[1][1] = boost::lexical_cast<double>(keywords[1]);
        std::getline(ifs, sline); // this line is not used
    }

    // set jastrow
    jastrow.set_F_long(F_long);
    jastrow.impose_cusp(false); // determine A_long

    // cutoff function is not read since TC++ does not use it.
    // So skip remaining lines before "TERM 2:" (if exists)
    while (std::getline(ifs, sline))
    {
        if (sline.find("TERM 2") != std::string::npos) { return true; } // "TERM 2" is found
    } // else, end of file
    return true;
}

bool io_qmc_files::read_POLY_casl(std::ifstream &ifs, const int num_independent_spins, Jastrow &jastrow, bool &cusp_poly)
{
    // polynomial parameters to be read
    int deg_poly; // degree of a polynomial
    std::vector<std::vector<double> > L_poly(2); // cutoff length
    for (int is1=0; is1<2; is1++) { L_poly[is1].resize(2); }
    std::vector<std::vector<std::vector<double> > > c_poly_input(2); // coefficients of a polynomial c[i] r^i (i=0,1,2...)
    for (int is1=0; is1<2; is1++) { c_poly_input[is1].resize(2); }

    std::string sline;
    std::vector<std::string> keywords;

    keywords = read_a_line_casl(ifs);
    if (keywords[0]!="Order") { error_messages::stop("Order is not given for polynomial terms (read_casl)"); }
    deg_poly = boost::lexical_cast<int>(keywords[1]);
    if (deg_poly<0) { error_messages::inappropriate_argument("Order: for polynomial terms in casl", deg_poly, ">=0"); }
    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            c_poly_input[is1][is2] = std::vector<double>(deg_poly, 0.0);
        }
    }

    if (read_a_line_casl(ifs)[0]!="e-e cutoff") { error_messages::stop("e-e cutoff is not given for polynomial terms (read_casl)"); }
    keywords = read_a_line_casl(ifs);
    if (keywords[0]!="Type" || keywords[1]!="polynomial") { error_messages::stop("Type: polynomial should be given for e-e cutoff of polynomial terms (read_casl)"); }
    if (read_a_line_casl(ifs)[0]!="Constants") { error_messages::stop("Constants are not given for e-e cutoff of polynomial terms (read_casl)"); }
    keywords = read_a_line_casl(ifs);
    if (keywords[0]!="C" || keywords[1]!="3") { error_messages::stop("C: 3 is not given for e-e cutoff of polynomial terms (read_casl)"); }
    if (read_a_line_casl(ifs)[0]!="Parameters") { error_messages::stop("Parameters are not given for e-e cutoff of polynomial terms (read_casl)"); }

    if (read_a_line_casl(ifs)[0]!="Channel 1-1") { error_messages::stop("Channel 1-1 is not given for polynomial cutoff (read_casl)"); }
    keywords = read_a_line_casl(ifs);
    if (keywords[0]!="L") { error_messages::stop("L is not given for polynomial cutoff (read_casl)"); }
    L_poly[0][0] = boost::lexical_cast<double>(keywords[1]);
    std::getline(ifs, sline); // this line is not used

    if (read_a_line_casl(ifs)[0]!="Channel 1-2") { error_messages::stop("Channel 1-2 is not given for polynomial cutoff (read_casl)"); }
    keywords = read_a_line_casl(ifs);
    if (keywords[0]!="L") { error_messages::stop("L is not given for polynomial cutoff (read_casl)"); }
    L_poly[0][1] = L_poly[1][0] = boost::lexical_cast<double>(keywords[1]);
    std::getline(ifs, sline); // this line is not used

    if (num_independent_spins==1)
    {
        L_poly[1][1] = L_poly[0][0];
    }
    else
    {
        if (read_a_line_casl(ifs)[0]!="Channel 2-2") { error_messages::stop("Channel 2-2 is not given for polynomial cutoff (read_casl)"); }
        keywords = read_a_line_casl(ifs);
        if (keywords[0]!="L") { error_messages::stop("L is not given for polynomial cutoff (read_casl)"); }
        L_poly[1][1] = boost::lexical_cast<double>(keywords[1]);
        std::getline(ifs, sline); // this line is not used
    }

    // cusp condition
    keywords = read_a_line_casl(ifs);
    if (keywords[0]!="e-e cusp") { error_messages::stop("e-e cusp is not given for polynomial terms (read_casl)"); }
    if (keywords[1]=="F")
    {
        cusp_poly = false;
    }
    else if (keywords[1]=="T")
    {
        cusp_poly = true;
    }
    else
    {
        error_messages::inappropriate_argument("e-e cusp: for polynomial terms in casl", keywords[1], "F or T");
    }


    // read linear parameters
    if (read_a_line_casl(ifs)[0]!="Linear parameters") { error_messages::stop("Linear parameters: is not given for polynomial terms (read_casl)"); }

    if (deg_poly>1)
    {
        if (read_a_line_casl(ifs)[0]!="Channel 1-1") { error_messages::stop("Channel 1-1 is not given for polynomial terms (read_casl)"); }
        for (int ideg=1; ideg<deg_poly; ideg++)
        {
            keywords = read_a_line_casl(ifs);
            std::string cname = "c_" + std::to_string(ideg+1); // c_2, c_3, ...
            if (keywords[0]!=cname) { error_messages::stop("An appropriate name for linear parameters is not given for polynomial terms (read_casl)"); }
            c_poly_input[0][0][ideg] = boost::lexical_cast<double>(keywords[1]);
        }
        
        if (read_a_line_casl(ifs)[0]!="Channel 1-2") { error_messages::stop("Channel 1-2 is not given for polynomial terms (read_casl)"); }
        for (int ideg=1; ideg<deg_poly; ideg++)
        {
            keywords = read_a_line_casl(ifs);
            std::string cname = "c_" + std::to_string(ideg+1); // c_2, c_3, ...
            if (keywords[0]!=cname) { error_messages::stop("An appropriate name for linear parameters is not given for polynomial terms (read_casl)"); }
            c_poly_input[0][1][ideg] = c_poly_input[1][0][ideg] = boost::lexical_cast<double>(keywords[1]);
        }
        
        if (num_independent_spins==1)
        {
            c_poly_input[1][1] = c_poly_input[0][0];
        }
        else
        {
            if (read_a_line_casl(ifs)[0]!="Channel 2-2") { error_messages::stop("Channel 2-2 is not given for polynomial terms (read_casl)"); }
            for (int ideg=1; ideg<deg_poly; ideg++)
            {
                keywords = read_a_line_casl(ifs);
                std::string cname = "c_" + std::to_string(ideg+1); // c_2, c_3, ...
                if (keywords[0]!=cname) { error_messages::stop("An appropriate name for linear parameters is not given for polynomial terms (read_casl)"); }
                c_poly_input[1][1][ideg] = boost::lexical_cast<double>(keywords[1]);
            }
        }
    } // if (deg_poly>1)

    // set polynomial-Jastrow parameters
    jastrow.set_polynomial_parameters(deg_poly, L_poly, c_poly_input, cusp_poly);

    // check whether "TERM 2" exists
    while (std::getline(ifs, sline))
    {
        if (sline.find("TERM 2") != std::string::npos) { return true; } // "TERM 2" is found
    } // else, end of file
    return true;
}
