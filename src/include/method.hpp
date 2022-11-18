// [class Method]
// method (HF, TC, BITC, FREE) and mode (SCF, BAND)

#ifndef TC_METHOD_HPP
#define TC_METHOD_HPP

class Method
{
private:
    std::string calc_method_;  // HF or TC or BITC or FREE
    std::string calc_mode_;    // SCF or BAND

public:
    const std::string &calc_method() const { return calc_method_; }
    const std::string &calc_mode() const { return calc_mode_; }

    Method() : calc_method_(""), calc_mode_("") {}
    void set_calc_method(const std::string &calc_method);
    void set_calc_mode(const std::string &calc_mode);
    void check_consistency_of_calc_mode(const std::string &calc_mode_QE) const;
    void bcast(const bool am_i_mpi_rank0);
};

#endif // TC_METHOD_HPP
