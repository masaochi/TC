// [class Spin]
// spin setting (no spin, noncollinear etc.), magnetization

#ifndef TC_SPIN_HPP
#define TC_SPIN_HPP

class Spin
{
private:

    // cf. in Quantum Espresso, nspin = 1 for no-spin, 2 for spin-polarized collinear, 4 for non-collinear
    // To summarize, (is_spinor_, num_independent_spins_) =
    // [no-spin: nspin = 1 in QE] (false, 1)
    // [spin-polarized collinear: nspin = 2 in QE] (false, 2)
    // [non-collinear: nspin = 4 in QE] (true, 1)
    bool is_spinor_;
    int num_independent_spins_; // (no-spin or non-collinear) 1, (spin-polarized spin-collinear) 2

// magnetization
// constraint?

public:
    bool is_spinor() const { return is_spinor_; }
    int num_independent_spins() const { return num_independent_spins_; }

    Spin() : is_spinor_(false), num_independent_spins_(-1)  {}

    // set spin variables from QE input (xml)
    void set(const std::string &noncolin, const std::string &lsda);
    void bcast(const bool am_i_mpi_rank0);
};


#endif // TC_SPIN_HPP

