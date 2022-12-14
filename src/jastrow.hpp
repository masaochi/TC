// [class Jastrow]
// jastrow parameters. "Jastrow jastrow" is a member of "Potential potential".

#ifndef TC_JASTROW_HPP
#define TC_JASTROW_HPP

class Jastrow
{
private:
    // Jastrow factor is defined as J = exp[-\sum_{i>j} u(x_i,x_j)] 

    // Jastrow parameters for u_long(|r-r'|) = (A_long/|r-r'|)*[1 - exp(-|r-r'|/F_long)]. 

    // In "input.in", "A" parameters are specified in the unit of sqrt(V/(4piN)), V:volume, N:no. of valence electrons,
    // while A_long in this class is in Hartree unit. i.e., sqrt(V/(4piN)) is multiplied with the values specified in input.in.

    // F_long is automatically determined so as to satisfy the cusp condition.
    // i.e., F = \sqrt(2A) for parallel spins, F = sqrt(A) for anti-parallel spins
    std::vector<std::vector<double> > A_long_; // e.g. A_long_[0][0]: up-up. [0][1]: up-down.
    std::vector<std::vector<double> > F_long_;

    std::vector<std::vector<double> > FourPI_A_long_; // 4*PI*A_long_
    std::vector<std::vector<double> > F2_long_; // F_long_*F_long_
    std::vector<std::vector<double> > Finv2_long_; // 1/(F_long_*F_long_)
    std::vector<std::vector<double> > PI_AAFinv2_long_; // PI*(A_long_/F_long)^2
    // bool uses_RPA_jastrow_; // false then u_long calc will be skipped

    // initialization, called in set()
    void set_A_long_unnormalized(const std::vector<std::vector<double> > &A_long);

public:
    const std::vector<std::vector<double> > &A_long() const { return A_long_; }
    const std::vector<std::vector<double> > &F_long() const { return F_long_; }

    Jastrow() : A_long_{{1.0, 1.0}, {1.0, 1.0}}, F_long_{{1.0, 1.0}, {1.0, 1.0}},
        FourPI_A_long_{{0.0, 0.0}, {0.0, 0.0}}, F2_long_{{1.0, 1.0}, {1.0, 1.0}},
        Finv2_long_{{1.0, 1.0}, {1.0, 1.0}}, PI_AAFinv2_long_{{0.0, 0.0}, {0.0, 0.0}} {}

    // initialization
    void set(const std::vector<std::vector<double> > &A_long);
    void normalize_A_long(const CrystalStructure &crystal_structure,
                          const BlochStates &bloch_states,
                          const Spin &spin);
    void bcast();

    // returns a value of the Jastrow function
    double uk(const Eigen::Vector3d &G, const int is1, const int is2) const;
    double tc_2body(const Eigen::Vector3d &G, const int is1, const int is2) const; // two-body pot.
};


#endif // TC_JASTROW_HPP
