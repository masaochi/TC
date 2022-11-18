// [class Jastrow]
// jastrow parameters. "Jastrow jastrow" is a member of "Potential potential".

#include "include/header.hpp"

// u = 4piA*[1/G2 - 1/(G2 + (1/F2))] = 4piA * (1/F2) / G2(G2 + (1/F2))
// = 4piA / (G2*(F2G2+1))
double Jastrow::uk(const Eigen::Vector3d &G, const int is1, const int is2) const
{
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

double Jastrow::tc_2body(const Eigen::Vector3d &G, const int is1, const int is2) const
{
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
