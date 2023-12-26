// [class FileNames]
// file names for input & output

#include "include/header.hpp"

void FileNames::set_qe_save_dir(const std::string &qe_save_dir)
{
    if (qe_save_dir=="") { error_messages::inappropriate_argument("qe_save_dir", qe_save_dir, "should not be blank"); }
    qe_save_dir_ = qe_save_dir;
}

void FileNames::set_qe_xml_file_name()
{
    if (qe_save_dir_=="") { error_messages::stop("Error: set_qe_xml_file_name should be called after set_qe_save_dir"); }
    qe_xml_ = qe_save_dir_ + "/data-file-schema.xml";
}

void FileNames::set_qe_wfc_file_names(const Spin &spin, const int &num_irreducible_kpoints)
{
    if (qe_save_dir_=="") { error_messages::stop("Error: set_qe_wfc_file_names should be called after set_qe_save_dir"); }

    const int num_independent_spins = spin.num_independent_spins();
    qe_wfc_.resize(num_independent_spins);
    for (int ispin=0; ispin<num_independent_spins; ispin++) 
    {
        qe_wfc_[ispin].resize(num_irreducible_kpoints);
    }
    
    if (num_independent_spins!=2) 
    {
        for (int ik=0; ik<num_irreducible_kpoints; ik++) 
        {
	    // qe_wfc_[0][ik] = qe_save_dir_ + "/wfc" + std::to_string(ik+1) + ".dat"; // ik+1: start from 1
	    std::ostringstream oss;
	    oss << ik+1;
            qe_wfc_[0][ik] = qe_save_dir_ + "/wfc" + oss.str() + ".dat"; // ik+1: start from 1
        }
    }
    else
    {
        for (int ik=0; ik<num_irreducible_kpoints; ik++) 
        {
	    // qe_wfc_[0][ik] = qe_save_dir_ + "/wfcup" + std::to_string(ik+1) + ".dat";
            // qe_wfc_[1][ik] = qe_save_dir_ + "/wfcdw" + std::to_string(ik+1) + ".dat";
	    std::ostringstream oss;
	    oss << ik+1;
            qe_wfc_[0][ik] = qe_save_dir_ + "/wfcup" + oss.str() + ".dat";
            qe_wfc_[1][ik] = qe_save_dir_ + "/wfcdw" + oss.str() + ".dat";
        }        
    }
}

void FileNames::set_pseudo_dir(const std::string &pseudo_dir)
{
    if (pseudo_dir=="") { error_messages::inappropriate_argument("pseudo_dir", pseudo_dir, "should not be blank"); }
    pseudo_dir_ = pseudo_dir;
}

void FileNames::set_upf_file_names(const std::vector<std::string> &upf)
{
    upf_.resize(upf.size());
    if (pseudo_dir_=="") { error_messages::stop("Error: set_upf_file_names should be called after set_pseudo_dir"); }
    for (int iatomic_species=0; iatomic_species<upf.size(); iatomic_species++)
    {
        upf_[iatomic_species] = pseudo_dir_ + "/" + upf[iatomic_species];
    }
}

void FileNames::set_reads_binary_false(std::ostream *ost)
{
    reads_binary_ = false;
    *ost << "  WARNING: reads_binary = false was set. This option should not be specified unless you perform test calculations in test folder." << std::endl;
}

