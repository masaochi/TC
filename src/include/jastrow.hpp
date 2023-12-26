// [class Jastrow]
// jastrow parameters. "Jastrow jastrow" is a member of "Potential potential".

#ifndef TC_JASTROW_HPP
#define TC_JASTROW_HPP

class Jastrow
{
private:
    // Jastrow factor is defined as J = exp[-\sum_{i>j} u(x_i,x_j)] 
    // Note: "-" sign...

    // (1) When A parameters are specified in "input.in": RPA-type with specified A values in the unit of sqrt(V/(4piN)).
    // (2) When Jastrow parameters are specified in "parameters.casl": RPA+polynomial-type with specified values in the atomic unit.
    // (3) Both (1) and (2) are given: stop with error.
    // (4) Neither (1) nor (2) are given: RPA-type with A = 1 in the unit of sqrt(V/(4piN)) is set as default.
    // In all cases, cusp condition is imposed on RPA terms unless no RPA terms are given in case (2) (i.e., polynomial Jastrow is used.)

    // **
    // RPA jastrow: u_long(|r-r'|) = (A_long/|r-r'|)*[1 - exp(-|r-r'|/F_long)]. 
    // **

    // In "input.in", "A" parameters are specified in the unit of sqrt(V/(4piN)), V:volume, N:no. of valence electrons,
    // while A_long in this class is in Hartree unit. i.e., sqrt(V/(4piN)) is multiplied with the values specified in input.in.

    // F_long is automatically determined so as to satisfy the cusp condition.
    // i.e., F = \sqrt(2A) for parallel spins, F = sqrt(A) for anti-parallel spins
    std::vector<std::vector<double> > A_long_; // e.g. A_long_[0][0]: up-up. [0][1]: up-down.
    std::vector<std::vector<double> > F_long_;

    std::vector<std::vector<double> > FourPI_A_long_; // 4*PI*A_long_
    std::vector<std::vector<double> > Finv_long_; // 1/F_long_
    std::vector<std::vector<double> > F2_long_; // F_long_*F_long_
    std::vector<std::vector<double> > Finv2_long_; // 1/(F_long_*F_long_)
    std::vector<std::vector<double> > PI_AAFinv2_long_; // PI*(A_long_/F_long)^2

    bool is_A_zero_; // whether all A_long_[is1][is2] are zero
    bool is_A_given_in_input_in_; // whether A is given in "input.in"

    // **
    // POLY (polynomial) jastrow: u_poly(|r-r'|) = -\sum_{i=1}^deg_poly c_poly_input,i |r-r'|^{i-1} * (1 - |r-r'|/L_poly)^3 * step_function(L_poly - |r-r'|)
    //                                           = \sum_{i=1}^{deg_poly+3} c_poly_internal,i |r-r'|^{i-1} * step_function(L_poly - |r-r'|)
    // NOTE! "-" sign in the first line and that in J=exp[-u] are canceled. Thus, coefficients are consistent with CASINO definition.
    // **
    int deg_poly_; // degree of a polynomial
    std::vector<std::vector<double> > L_poly_; // cutoff length: L_poly[up(0),dn(1)][up(0),dn(1)]
    std::vector<std::vector<std::vector<double> > > c_poly_input_; // input coefficients of a polynomial: c_poly[up(0),dn(1)][up(0),dn(1)][0 -- deg_poly-1]
    std::vector<std::vector<std::vector<double> > > c_poly_internal_; // including a damping function and a minus sign (see above)

    std::vector<std::vector<std::vector<double> > > FourPI_power_L_; // FourPI_power_L_[is1][is2][i] = 4 * PI * L_poly[is1][is2]^i

    const int max_deg_poly_; // maximum degree of a polynomial
    // contants for evaluation of J_integ[n] = \int_0^L r^n exp[-iG \dot r] d^3r = (4pi/G)*Im[\int r^{n+1} exp(iGr) dr]
    const double threshold_Taylor_; // threshold for Taylor expansion
    const int max_deg_Taylor_; // maximum degree of Taylor expansion

    // RPA+POLY parameters
    std::vector<std::vector<double> > exp_LF_; // exp(-L_poly/F_long)
    std::vector<std::vector<std::vector<double> > > M_integ_; // see below, used for evaluating J2_integ

    // constants for evaluating M_integ[n] = L^{-n} \int_0^L exp(-r/F)*r^l dr
    const double threshold_Taylor2_; // threshold for Taylor expansion
    const int max_deg_Taylor2_; // maximum degree of Taylor expansion

    // variables for evaulating J3_integ
    int num_rmesh_J3_integ_;
    std::vector<std::vector<Eigen::VectorXd> > rmesh_J3_integ_;
    std::vector<std::vector<Eigen::VectorXd> > integrand_stored_J3_integ_;

    // RPA
    double uk_RPA(const Eigen::Vector3d &G, const int is1, const int is2) const;
    double tc_2body_RPA(const Eigen::Vector3d &G, const int is1, const int is2) const; // two-body pot.

    // polynomial
    double uk_POLY(const Eigen::Vector3d &G, const int is1, const int is2) const;
    double tc_2body_POLY(const Eigen::Vector3d &G, const int is1, const int is2) const; // two-body pot.
    void calculate_J_integ(const double G2, const int is1, const int is2, std::vector<double> &J_integ, const bool calculates_zeroth) const;

    // RPA+polynomial
    double tc_2body_RPAPOLY(const Eigen::Vector3d &G, const int is1, const int is2) const; // two-body pot.
    void calculate_J2_integ(const double G1, const int is1, const int is2, std::vector<double> &J2_integ) const;
    double J3_integ(const double G1, const int is1, const int is2) const;
    void prepare_for_J3_integ(const CrystalStructure &crystal_structure,
                              const PlaneWaveBasis &plane_wave_basis);
    // following functions are called in prepare_for_J3_integ()
    void test_J3_integ(const int num_rmesh_J3_integ, const std::vector<double> &Glist,
                       std::vector<std::vector<std::vector<double> > > & integral);
    double test_error_J3_integ(const std::vector<std::vector<std::vector<double> > > &integral1,
                               const std::vector<std::vector<std::vector<double> > > &integral2) const;
    void set_arrays_J3_integ(const int num_rmesh_J3_integ);

public:
    const std::vector<std::vector<double> > &A_long() const { return A_long_; }
    const std::vector<std::vector<double> > &F_long() const { return F_long_; }
    int deg_poly() const { return deg_poly_; }
    const std::vector<std::vector<double> > &L_poly() const { return L_poly_; }
    const std::vector<std::vector<std::vector<double> > > &c_poly_input() const { return c_poly_input_; }
    bool is_A_zero() const { return is_A_zero_; }
    bool is_A_given_in_input_in() const { return is_A_given_in_input_in_; }

    Jastrow() : 
        A_long_{{1.0, 1.0}, {1.0, 1.0}}, 
        F_long_{{1.0, 1.0}, {1.0, 1.0}},
        FourPI_A_long_{{0.0, 0.0}, {0.0, 0.0}},
        Finv_long_{{1.0, 1.0}, {1.0, 1.0}},
        F2_long_{{1.0, 1.0}, {1.0, 1.0}},
        Finv2_long_{{1.0, 1.0}, {1.0, 1.0}},
        PI_AAFinv2_long_{{0.0, 0.0}, {0.0, 0.0}},
        is_A_zero_(false),
        is_A_given_in_input_in_(false),
        deg_poly_(0),
        L_poly_{{1.0, 1.0}, {1.0, 1.0}},
        FourPI_power_L_{{{1.0}, {1.0}}, {{1.0}, {1.0}}},
        max_deg_poly_(16),
        threshold_Taylor_(1000.0), // can be 800, 900 etc. for cheaper calc.
        max_deg_Taylor_(50),
        exp_LF_{{1.0, 1.0}, {1.0, 1.0}},
        threshold_Taylor2_(1.5),
        max_deg_Taylor2_(30) {}

    // initialization (RPA)
    void set_A_long(const std::vector<std::vector<double> > &A_long);
    void unnormalize_A_long(const double &volume,
                            const double &num_electrons,
                            const Spin &spin);
    void set_F_long(const std::vector<std::vector<double> > &F_long);
    void impose_cusp(const bool A_to_F);
    void set_derived_RPA_parameters();
    void set_is_A_given_in_input_in(const bool is_A_given_in_input_in);
    void bcast_A_long();
    void bcast_F_long();
    void print_RPA_parameters(std::ostream *ost) const;

    // initialization (polynomial parameters)
    void set_polynomial_parameters(const int &deg_poly,
                                   const std::vector<std::vector<double> > &L_poly,
                                   const std::vector<std::vector<std::vector<double> > > &c_poly_input,
                                   const bool cusp_poly);
    void set_derived_polynomial_parameters();
    void set_derived_RPA_polynomial_parameters(const CrystalStructure &crystal_structure,
                                               const PlaneWaveBasis &plane_wave_basis);
    void bcast_polynomial_parameters(const bool am_i_mpi_rank0);
    void print_polynomial_parameters(std::ostream *ost) const;

    // initialization (RPA + polynomial)
    void print_RPA_polynomial_parameters(std::ostream *ost) const;

    // returns a value of the Jastrow function
    double uk(const Eigen::Vector3d &G, const int is1, const int is2) const;
    double tc_2body(const Eigen::Vector3d &G, const int is1, const int is2) const; // two-body pot.

    // Jastrow function u -> -u
    void take_Hermitian_conj();
};


#endif // TC_JASTROW_HPP
