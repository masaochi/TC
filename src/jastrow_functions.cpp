// [class Jastrow]
// jastrow parameters. "Jastrow jastrow" is a member of "Potential potential".

#include "include/header.hpp"

double Jastrow::uk(const Eigen::Vector3d &G, const int is1, const int is2) const
{
    double val = 0.0;
    if (!is_A_zero_)
    {
        val += uk_RPA(G, is1, is2);
    }
    if (deg_poly_>0)
    {
        val += uk_POLY(G, is1, is2);
    }
    return val;
}

double Jastrow::tc_2body(const Eigen::Vector3d &G, const int is1, const int is2) const
{
    if (is_A_zero_) // no RPA
    {
        if (deg_poly_>0) // polynomial only
        {
            return tc_2body_POLY(G, is1, is2);
        }
        else // no Jastrow (Coulomb only)
        {
            double G2 = G.squaredNorm();
            if (G2>1e-12) 
            {
                return FourPI/G2;
            }
            else
            {
                return 0.0;
            }
        }
    }
    else // with RPA
    {
        if (deg_poly_>1) // RPA + polynomial // NOTE! deg_poly_=1 then polynomial should be zero to be cuspless.
        {
            return tc_2body_RPAPOLY(G, is1, is2);
        }
        else // RPA only
        {
            return tc_2body_RPA(G, is1, is2);
        }
    }
}

// for RPA Jastrow

// u = 4piA*[1/G2 - 1/(G2 + (1/F2))] = 4piA * (1/F2) / G2(G2 + (1/F2))
// = 4piA / (G2*(F2G2+1))
double Jastrow::uk_RPA(const Eigen::Vector3d &G, const int is1, const int is2) const
{
    assert(is1>=0 && is1<=1 && is2>=0 && is2<=1);

    double G2 = G.squaredNorm();
    if (G2<1e-12)
    {
        // subtract 4piA/G2 then = -4piA/(G2 + (1/F2)) -> -4piA*F*F (G=0)
        return -FourPI_A_long_[is1][is2]*F2_long_[is1][is2];
    }
    else
    {
        return FourPI_A_long_[is1][is2]/(G2*(G2*F2_long_[is1][is2]+1));
    }
}

double Jastrow::tc_2body_RPA(const Eigen::Vector3d &G, const int is1, const int is2) const
{
    assert(is1>=0 && is1<=1 && is2>=0 && is2<=1);

    double G2 = G.squaredNorm();
    if (G2<1e-12)
    {
        // contribution from "\nabla^2 u(G=0)" is always canceled with 
        // G=0 contribution of "\nabla u \cdot \nabla".
        // Thus, here we only consider "-(\nabla u)^2".
        return -2*PI_AAFinv2_long_[is1][is2]*F_long_[is1][is2]; // -2*PI*A^2/F
    }
    else
    {
        double G1 = std::sqrt(G2);
        double FG1 = F_long_[is1][is2]*G1;
        double F2G2 = F2_long_[is1][is2]*G2;
        return FourPI/G2 // coulomb
            - FourPI_A_long_[is1][is2]/(1+F2G2) // + \nabla^2 u (= -G^2*uk)
            + (PI_AAFinv2_long_[is1][is2]/G1)*(PI*F2G2 // - (\nabla u)^2 (these three lines)
                                              +(4+2*F2G2)*std::atan(FG1/2)
                                              -(4+4*F2G2)*std::atan(FG1));
    }
}

// for polynomial Jastrow

// calculate J_integ[0:N_degree-1] where N_degree = J_integ.size()
// J_integ[n] = \int r^{n-1} exp[-iG \dot r] d^3r = (4pi/G)*Im[\int r^n exp(iGr) dr]
// When calculates_zeroth==false, J_integ[0] is NOT calculated.
void Jastrow::calculate_J_integ(const double G2, const int is1, const int is2, std::vector<double> &J_integ, const bool calculates_zeroth) const
{
    assert(is1>=0 && is1<=1 && is2>=0 && is2<=1);

    const int N_degree = J_integ.size();
    const double L = L_poly_[is1][is2];

    assert(G2>-1e-10); // not negative
    assert(L>1e-10); // positive
    assert(N_degree<37); // If you would like to use a larger N_degree, max_deg_Taylor_ and threshold_Taylor_ should be increased.

    double G2L2 = G2*L*L;
    int ipoly_start = calculates_zeroth ? 0 : 1;
    if (G2L2 < threshold_Taylor_)
    {
        // Taylor expansion for J_integ[n]
        // = 1/(n+2) - (GL)^2/(3!(n+4)) + (GL)^4/(5!(n+6)) - (GL)^6/(7!(n+8)) +...
        // = 1/(n+2) - (GL)^2/(3*2) * [ 1/(n+4) - (GL)^2/(5*4) * [ 1/(n+6) - (GL)^2/(7*6) * ...
        for (int ipoly=ipoly_start; ipoly<N_degree; ipoly++)
        {
            boost::multiprecision::cpp_dec_float_50 temporary_Taylor = 0;
            for (int jpoly=0; jpoly<max_deg_Taylor_; jpoly++)
            {
                int deg = 2*(max_deg_Taylor_ - jpoly);
                temporary_Taylor
                    = -temporary_Taylor * G2L2/((deg+1)*deg)
                    + 1.0/static_cast<boost::multiprecision::cpp_dec_float_50>(ipoly+deg);
            }
            // not a clever way to convert cpp_dec_float_50 to double...
            std::string temporary_Taylor_str = temporary_Taylor.str(0, std::ios_base::scientific);
            J_integ[ipoly] = std::stod(temporary_Taylor_str);
/*
            double temporary_Taylord = 0;
            for (int jpoly=0; jpoly<max_deg_Taylor_; jpoly++)
            {
                int deg = 2*(max_deg_Taylor_ - jpoly);
                temporary_Taylord
                    = -temporary_Taylord * G2L2/((deg+1)*deg)
                    + 1.0/static_cast<double>(ipoly+deg);
            }
            J_integ[ipoly] = temporary_Taylord;
*/
        } // ipoly
    }
    else
    {
        double GL = std::sqrt(G2L2);
        double cos = std::cos(GL);
        double sin = std::sin(GL);
        for (int ipoly=ipoly_start; ipoly<N_degree; ipoly++)
        {
            // 1) cosine term
            double d_tmp = -1.0/(G2L2);
            int n_tmp = ipoly;
            double coeff = d_tmp; // 0 for the case when the first order term is not required (canceled each other in the cutoff function)
            for (int jpoly=0; jpoly<ipoly/2; jpoly++) // ipoly = 0, 1, 2, 3,...; ipoly/2 = 0, 0, 1, 1,...
            {
                // e.g.) jpoly = 0, 1,...; numerator = n(n-1), (n-2)(n-3),... for ipoly=n.
                d_tmp *= -n_tmp*(n_tmp - 1)/G2L2;
                n_tmp -= 2;
                coeff += d_tmp;
            }
            J_integ[ipoly] = coeff * cos;
          
            // 2) const. term
            if (ipoly%2 == 0) { J_integ[ipoly] -= d_tmp; }

            // (3) sine term
            if (ipoly>0)
            {
                d_tmp = ipoly/(G2L2*GL);
                n_tmp = ipoly-1;
                coeff = d_tmp; // 0 for the case when the first order term is not required (canceled each other in the cutoff function)
                for (int jpoly=0; jpoly<(ipoly-1)/2; jpoly++) // ipoly = 0, 1, 2, 3,...; (ipoly-1)/2 = *, 0, 0, 1, 1,...
                {
                    // e.g.) jpoly = 0, 1, 2,...; numerator = (n-1)(n-2), (n-3)(n-4),... for ipoly=n.
                    d_tmp *= -n_tmp*(n_tmp - 1)/G2L2;
                    n_tmp -= 2;
                    coeff += d_tmp;
                }
                J_integ[ipoly] += coeff * sin;
            }
        } // ipoly
    } // if

    // multiply 4pi*L^{n+2}
    for (int ipoly=ipoly_start; ipoly<N_degree; ipoly++)
    {
        J_integ[ipoly] *= FourPI_power_L_[is1][is2][ipoly+2];
    }
}

double Jastrow::uk_POLY(const Eigen::Vector3d &G, const int is1, const int is2) const
{
    assert(is1>=0 && is1<=1 && is2>=0 && is2<=1);

    double G2 = G.squaredNorm();
    std::vector<double> J_integ(deg_poly_+4); // use J_integ[0:N+3] 
    calculate_J_integ(G2, is1, is2, J_integ, false); // false: J_integ[0] is not calculated

    double return_val = 0.0;
    for (int ipoly=0; ipoly<c_poly_internal_[is1][is2].size(); ipoly++)
    {
        return_val += c_poly_internal_[is1][is2][ipoly] * J_integ[ipoly+1];
    }
    return return_val;
/*
    // test for long-double
    boost::multiprecision::cpp_dec_float_50 return_val = 0;
    for (int ipoly=0; ipoly<c_poly_internal_[is1][is2].size(); ipoly++)
    {
        return_val += c_poly_internal_[is1][is2][ipoly] * J_integ[ipoly+1];
    }
    return std::stod(return_val.str(0, std::ios_base::scientific));
*/
}

double Jastrow::tc_2body_POLY(const Eigen::Vector3d &G, const int is1, const int is2) const
{
    assert(is1>=0 && is1<=1 && is2>=0 && is2<=1);

    double G2 = G.squaredNorm();
    std::vector<double> J_integ(2*(deg_poly_+C_damp_poly_)-2); // use J_integ[0:2N+3] for C_damp_poly=3, [0:2N+1] for C_damp_poly=2.
    calculate_J_integ(G2, is1, is2, J_integ, false); // false: J_integ[0] is not calculated

    // Coulomb
    double return_val_sub1 = 0.0;
    if (G2>1e-12) { return_val_sub1 += FourPI/G2; }

    // +\nabla^2 u (= -G^2 uk)
    double return_val_sub2 = 0.0;
    for (int ipoly=0; ipoly<c_poly_internal_[is1][is2].size(); ipoly++) 
    {
        return_val_sub2 -= G2 * c_poly_internal_[is1][is2][ipoly] * J_integ[ipoly+1];
    }

    // -(\nabla u)^2
/*
    double return_val_sub3 = 0;
    for (int ipoly=1; ipoly<c_poly_internal_[is1][is2].size(); ipoly++) // ipoly=0: no contribution
    {
        double return_val_sub4 = 0;
        for (int jpoly=1; jpoly<c_poly_internal_[is1][is2].size(); jpoly++) // jpoly=0: no contribution
        {
            return_val_sub4 -= ipoly * jpoly
                * c_poly_internal_[is1][is2][ipoly]
                * c_poly_internal_[is1][is2][jpoly]
                * J_integ[ipoly+jpoly-1];
        }
        return_val_sub3 += return_val_sub4;
    }
    return return_val_sub1 + return_val_sub2 + return_val_sub3;
*/
    boost::multiprecision::cpp_dec_float_50 return_val_sub3 = 0;
    for (int ipoly=1; ipoly<c_poly_internal_[is1][is2].size(); ipoly++) // ipoly=0: no contribution
    {
        for (int jpoly=1; jpoly<c_poly_internal_[is1][is2].size(); jpoly++) // jpoly=0: no contribution
        {
            return_val_sub3 -= ipoly * jpoly
                * c_poly_internal_[is1][is2][ipoly]
                * c_poly_internal_[is1][is2][jpoly]
                * J_integ[ipoly+jpoly-1];
        }
    }
    return return_val_sub1 + return_val_sub2 + std::stod(return_val_sub3.str(0, std::ios_base::scientific));
}

// for RPA+polynomial
// we assume c_poly_internal_[*][*][0,1] = 0 since J_integ[-1] is difficult to evaluate.
double Jastrow::tc_2body_RPAPOLY(const Eigen::Vector3d &G, const int is1, const int is2) const
{
    assert(is1>=0 && is1<=1 && is2>=0 && is2<=1);
    assert(deg_poly_>1); // Note: for deg_poly_=1, polynomial Jastrow = 0 (to satisfy the cuspless condition)

    double G2 = G.squaredNorm();
    std::vector<double> J_integ(2*(deg_poly_+C_damp_poly_)-2); // use J_integ[0:2N+3] for C_damp_poly=3, [0:2N+1] for C_damp_poly=2.
    calculate_J_integ(G2, is1, is2, J_integ, true); // true: J_integ[0] is calculated

    // Coulomb + \nabla^2 u (RPA) - (\nabla u)^2 (RPA x RPA)
    double return_val_sub1 = 0.0;
    return_val_sub1 += tc_2body_RPA(G, is1, is2);

    // +\nabla^2 u (= -G^2 uk): polynomial
    double return_val_sub2 = 0.0;
    for (int ipoly=0; ipoly<c_poly_internal_[is1][is2].size(); ipoly++) 
    {
        return_val_sub2 -= G2 * c_poly_internal_[is1][is2][ipoly] * J_integ[ipoly+1];
    }

    // -(\nabla u)^2: polynomial x polynomial
    boost::multiprecision::cpp_dec_float_50 return_val_sub3 = 0;
    for (int ipoly=1; ipoly<c_poly_internal_[is1][is2].size(); ipoly++) // ipoly=0: no contribution
    {
        for (int jpoly=1; jpoly<c_poly_internal_[is1][is2].size(); jpoly++) // jpoly=0: no contribution
        {
            return_val_sub3 -= ipoly * jpoly
                * c_poly_internal_[is1][is2][ipoly]
                * c_poly_internal_[is1][is2][jpoly]
                * J_integ[ipoly+jpoly-1];
        }
    }

    // -(\nabla u)^2: RPA x polynomial
    double G1 = std::sqrt(G2);
    std::vector<double> J2_integ(deg_poly_+C_damp_poly_-1); // use J2_integ[0:N+1] for C_damp_poly=3, [0:N] for C_damp_poly=2
    calculate_J2_integ(G1, is1, is2, J2_integ);

    double return_val_sub4 = 0.0;
    for (int ipoly=2; ipoly<c_poly_internal_[is1][is2].size(); ipoly++)
    {
        return_val_sub4 += 2 * ipoly * c_poly_internal_[is1][is2][ipoly]
            * A_long_[is1][is2] * (J_integ[ipoly-2]
                                   -  (J2_integ[ipoly-1]/F_long_[is1][is2] + J2_integ[ipoly-2]));
    }
    // -(\nabla u)^2: RPA x polynomial (ipoly=1). 
    return_val_sub4 -= 2 * c_poly_internal_[is1][is2][1]
        * A_long_[is1][is2] * J3_integ(G1, is1, is2);
    // Note! ipoly=0 has no contribution for "RPA x polynomial" of -(\nabla u)^2. (since \nabla r^0 = 0)

    return return_val_sub1 + return_val_sub2 + std::stod(return_val_sub3.str(0, std::ios_base::scientific)) + return_val_sub4;
}

// calculate J2_integ[0:N_degree-1] where N_degree = J2_integ.size()
// J2_integ[n] = \int r^{n-1} exp[-iG \dot r] exp[-r/F] d^3r = (4pi/G)*Im[\int r^n exp(iGr) exp(-r/F) dr]
void Jastrow::calculate_J2_integ(const double G1, const int is1, const int is2, std::vector<double> &J2_integ) const
{
    assert(is1>=0 && is1<=1 && is2>=0 && is2<=1);

    const int N_degree = J2_integ.size();
    const double L = L_poly_[is1][is2];

    assert(G1>-1e-10); // not negative
    assert(L>1e-10); // positive

    double G2L2 = (G1*L)*(G1*L);
    if (G2L2 < threshold_Taylor2_) // use Taylor expansion for a small G
    {
        for (int ipoly=0; ipoly<N_degree; ipoly++)
        {
            J2_integ[ipoly] = 0.0;
            double coeff = 1.0;
            for (int jpoly=0; jpoly<max_deg_Taylor2_/2; jpoly++)
            {
                J2_integ[ipoly] += coeff * M_integ_[is1][is2][ipoly + 2*jpoly + 1];
                coeff *= (-G2L2)/static_cast<double>((2*jpoly+3) * (2*jpoly+2));
            }
            J2_integ[ipoly] *= FourPI_power_L_[is1][is2][ipoly+1];
/*
            boost::multiprecision::cpp_dec_float_50 coeff = 1.0;
            boost::multiprecision::cpp_dec_float_50 temporary_Taylor = 0.0;
            for (int jpoly=0; jpoly<max_deg_Taylor2_/2; jpoly++)
            {
                temporary_Taylor += coeff * M_integ_[is1][is2][ipoly + 2*jpoly + 1];
                coeff *= (-G2L2)/static_cast<boost::multiprecision::cpp_dec_float_50>((2*jpoly+3) * (2*jpoly+2));
            }
            temporary_Taylor *= FourPI_power_L_[is1][is2][ipoly+1];
            J2_integ[ipoly] = std::stod( temporary_Taylor.str(0, std::ios_base::scientific) );
*/
        }
    }
    else
    {
        Complex ctmp1 = 1.0/(Finv_long_[is1][is2] - I*G1); // 1/(1/F - iG)
        Complex exp2 = exp_LF_[is1][is2] * std::exp(I*G1*L); // exp(-L/F) * exp(iGL)
        
        std::vector<Complex> J2_integ_complex(N_degree);
        J2_integ_complex[0] = FourPI * ctmp1 * (1.0 - exp2);
        for (int ipoly=1; ipoly<N_degree; ipoly++)
        {
            J2_integ_complex[ipoly]
                = ctmp1 * (ipoly * J2_integ_complex[ipoly-1]
                           - exp2 * FourPI_power_L_[is1][is2][ipoly]);
        }
        for (int ipoly=0; ipoly<N_degree; ipoly++)
        {
            J2_integ[ipoly] = J2_integ_complex[ipoly].imag() / G1;
        }
    } // if (G2L2 < threshold_Taylor2_)
}

// J3_integ = (4pi/G)*\int_0^{L/F} dr sin(FGr)*[exp(-r) + (exp(-r) - 1)/r]
// sin(FGr)/G is evaluated here but others are stored as integrand_stored_J3_integ_.
double Jastrow::J3_integ(const double G1, const int is1, const int is2) const
{
    double FG = F_long_[is1][is2] * G1;

    double return_val = 0.0;
    for (int ir=0; ir<rmesh_J3_integ_[is1][is2].size(); ir++)
    {
        double FGr = FG * rmesh_J3_integ_[is1][is2](ir);
        if (FGr<0.1) // use Taylor expansion
        {
            // sin(FGr)/G = Fr - (Fr)^3 G^2/3! + (Fr)^5 G^4/5! - ...
            double sin_val = 0.0;
            double each_term = F_long_[is1][is2] * rmesh_J3_integ_[is1][is2](ir); 
            for (int i=0; i<3; i++)
            {
                sin_val += each_term;
                each_term *= -FGr * FGr / static_cast<double>((2*i+3)*(2*i+2));
            }
            return_val += integrand_stored_J3_integ_[is1][is2](ir) * sin_val;
        }
        else
        {
            return_val += integrand_stored_J3_integ_[is1][is2](ir) 
                * std::sin(FGr) / G1;
        }
    }
    return return_val;
}

// Take Hermitian conjugate of the Jastrow factor. To say, u -> -u.
// Should be modified if you add additional Jastrow parameters!
void Jastrow::take_Hermitian_conj()
{
    assert(A_long_.size()==2 && A_long_[0].size()==2 && A_long_[1].size()==2);
    assert(FourPI_A_long_.size()==2 && FourPI_A_long_[0].size()==2 
           && FourPI_A_long_[1].size()==2);
    assert(c_poly_internal_.size()==2 && c_poly_internal_[0].size()==2
           && c_poly_internal_[1].size()==2);

    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            // RPA terms
            A_long_[is1][is2] *= -1;
            FourPI_A_long_[is1][is2] *= -1;

            // polynomial terms
            for (int ideg=0; ideg<c_poly_internal_[is1][is2].size(); ideg++)
            {
                c_poly_internal_[is1][is2][ideg] *= -1;
            }
        }
    }
}
