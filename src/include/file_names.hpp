// [class FileNames]
// file names for input & output

#ifndef TC_FILE_NAMES_HPP
#define TC_FILE_NAMES_HPP

class FileNames
{
private:
    // "reads_binary_" is applied to qe_wfc. When reads_binary_==false, read non-binary files. (default: true)
    // This variables is prepared only for test calculation where environment-dependent binary files are not appropriate...
    bool reads_binary_;

    const std::string tc_input_;  // TC++ input file name (input.in)
    const std::string tc_output_;  // TC++ output file name (output.out)

    const std::string tc_wfc_scf_; // TC++ output (SCF wave functions) binary file name (tc_wfc_scf.dat)
    const std::string tc_wfc_band_; // TC++ output (BAND wave functions) binary file name (tc_wfc_band.dat)
    const std::string tc_energy_scf_; // TC++ output (SCF orbital energies) binary file name (tc_energy_scf.dat)
    const std::string tc_energy_band_; // TC++ output (BAND orbital energies) binary file name (tc_energy_band.dat)
    const std::string tc_scfinfo_; // TC++ output. SCF information that will be used in BAND calculation
    // TC++ output. BAND eigenvalues for band plot.
    const std::string tc_bandplot_; // for no-spin or spin-non-collinear 
    const std::string tc_bandplot_up_; // for up-band
    const std::string tc_bandplot_dn_; // for down-band

    // Address of the QE (quantum-espresso) "save" directory. e.g. "/home/hoge/QE/Si/prefix.save"
    // xml and wfc files are placed in the QE "save" directory specified here.
    std::string qe_save_dir_;

    std::string qe_xml_; //  e.g. "/home/hoge/QE/Si/prefix.save/data-file-schema.xml"

    // (spin-polarized, collinaer) qe_wfc_[is][ik]  is = 0 (up) or 1 (down), ik = irreducible k-point
    // e.g. "/home/hoge/QE/Si/prefix.save/wfcup3.dat" (is=0, ik=3) ".../wfcdw4.dat" (is=1, ik=4)
    // (no-spin or non-collinear) qe_wfc_[is][ik]  is = 0, ik = irreducible k-point
    // e.g. "/home/hoge/QE/Si/prefix.save/wfc3.dat" (is=0, ik=3)
    std::vector<std::vector<std::string> > qe_wfc_;

    // upf (pseudopotential) files are placed in "pseudo_dir_" directory.
    std::string pseudo_dir_;
    std::vector<std::string> upf_; // e.g. pseudo_dir_ + "Si.upf". upf_[index of atomic species]

public:
    bool reads_binary() const { return reads_binary_; }
    const std::string &tc_input() const { return tc_input_; }
    const std::string &tc_output() const { return tc_output_; }
    const std::string &tc_wfc_scf() const { return tc_wfc_scf_; }
    const std::string &tc_wfc_band() const { return tc_wfc_band_; }
    const std::string &tc_energy_scf() const { return tc_energy_scf_; }
    const std::string &tc_energy_band() const { return tc_energy_band_; }
    const std::string &tc_scfinfo() const { return tc_scfinfo_; }
    const std::string &tc_bandplot() const { return tc_bandplot_; }
    const std::string &tc_bandplot_up() const { return tc_bandplot_up_; }
    const std::string &tc_bandplot_dn() const { return tc_bandplot_dn_; }
    const std::string &qe_save_dir() const { return qe_save_dir_; }
    const std::string &qe_xml() const { return qe_xml_; }
    const std::vector<std::vector<std::string> > &qe_wfc() const { return qe_wfc_; }
    const std::string &pseudo_dir() const {return pseudo_dir_; }
    const std::vector<std::string> &upf() const { return upf_; }

    FileNames() : 
        reads_binary_(true),
        tc_input_("input.in"), 
        tc_output_("output.out"), 
        tc_wfc_scf_("tc_wfc_scf.dat"),
        tc_wfc_band_("tc_wfc_band.dat"),
        tc_energy_scf_("tc_energy_scf.dat"),
        tc_energy_band_("tc_energy_band.dat"),
        tc_scfinfo_("tc_scfinfo.dat"),
        tc_bandplot_("tc_bandplot.dat"),
        tc_bandplot_up_("tc_bandplot_up.dat"),
        tc_bandplot_dn_("tc_bandplot_dn.dat"),
        qe_save_dir_(""),
        pseudo_dir_("") { }

    void set_qe_save_dir(const std::string &qe_save_dir);
    void set_qe_xml_file_name(); // NOTE! should be called after set_qe_save_dir
    void set_qe_wfc_file_names(const Spin &spin, const int &num_irreducible_kpoints); // NOTE! should be called after set_qe_save_dir
    void set_pseudo_dir(const std::string &pseudo_dir);
    void set_upf_file_names(const std::vector<std::string> &upf); // NOTE! should be called after set_pseudo_dir
    void set_reads_binary_false(std::ostream *ost);
};

#endif // TC_FILE_NAMES_HPP
