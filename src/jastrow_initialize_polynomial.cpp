// [class Jastrow]
// jastrow parameters. "Jastrow jastrow" is a member of "Potential potential".

#include "include/header.hpp"

void Jastrow::set_polynomial_parameters(const int &deg_poly,
                                        const int &C_damp_poly,
                                        const std::vector<std::vector<double> > &L_poly,
                                        const std::vector<std::vector<std::vector<double> > > &c_poly_input,
                                        const bool cusp_poly)
{
    assert(L_poly.size()==2 && L_poly[0].size()==2 && L_poly[1].size()==2);
    assert(c_poly_input.size()==2 && c_poly_input[0].size()==2 && c_poly_input[1].size()==2);

    if (deg_poly<0) { error_messages::inappropriate_argument("deg_poly (a degree of a polynomial Jastrow)", deg_poly, "should be non-negative"); }
    if (deg_poly>max_deg_poly_) { error_messages::inappropriate_argument("deg_poly (a degree of a polynomial Jastrow)", deg_poly, "should be <= max_deg_poly_(16) in include/jastrow.hpp"); } 
    if (C_damp_poly!=2 && C_damp_poly!=3) { error_messages::inappropriate_argument("C_damp_poly (C in the damping function for the polynomial Jastrow)", C_damp_poly, "should be 2 or 3"); }
    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            if (L_poly[is1][is2]<1e-8) 
            {
                error_messages::inappropriate_argument("L_poly (cutoff length of a polynomial Jastrow)", L_poly[is1][is2], "should be positive"); 
            }

            // check consistency betwee c_poly_input[is1][is2].size() and deg_poly
            if (c_poly_input[is1][is2].size() != deg_poly)
            {
                error_messages::stop("the num. of coefficients of a polynomial Jastrow and its degree (deg_poly) is inconsistent."); 
            }
        }
    }

    // set input parameters
    deg_poly_ = deg_poly;
    C_damp_poly_ = C_damp_poly;
    L_poly_ = L_poly;
    c_poly_input_ = c_poly_input;

    // set the lowest order coefficients based on the cusp condition
    if (deg_poly==0)
    {
        if (cusp_poly) { error_messages::inappropriate_argument("deg_poly (a degree of polynomial Jastrow)", deg_poly, " should be >=1 when the cusp condition is imposed"); }
    }
    else // deg_poly>0
    {
        for (int is1=0; is1<2; is1++)
        {
            for (int is2=0; is2<2; is2++)
            {
                double L = L_poly_[is1][is2];
                double cusp_const = is1==is2 ? 0.25 : 0.5;
                if (!cusp_poly) { cusp_const = 0.0; } // cuspless
                double higher_order_coeff = deg_poly==1 ? 0.0 : c_poly_input_[is1][is2][1];

                c_poly_input_[is1][is2][0] = (L/static_cast<double>(C_damp_poly))*(higher_order_coeff - cusp_const);
            }
        }
    }

    // set c_poly_internal_: polynomial coefficients including a damping function (1-r/L)^3 and a minus sign
    // c_poly_internal(i) = -(c(i) - (3/L)c_{i-1} + (3/L^2)c_{i-2} - (1/L^3)c_{i-3})
    c_poly_internal_.resize(2);
    for (int is1=0; is1<2; is1++)
    {
        c_poly_internal_[is1].resize(2);
        for (int is2=0; is2<2; is2++)
        {
            double L = L_poly_[is1][is2];

            c_poly_internal_[is1][is2].resize(deg_poly+C_damp_poly);
            for (int ideg=0; ideg<deg_poly+C_damp_poly; ideg++)
            {
                c_poly_internal_[is1][is2][ideg] = 0.0;
            }

            for (int ideg=0; ideg<deg_poly; ideg++)
            {
                if (C_damp_poly==3)
                {
                    c_poly_internal_[is1][is2][ideg] -= c_poly_input_[is1][is2][ideg]; // c_i
                    c_poly_internal_[is1][is2][ideg+1] += 3.0*c_poly_input_[is1][is2][ideg]/L; // (3/L)*c_{i-1}
                    c_poly_internal_[is1][is2][ideg+2] -= 3.0*c_poly_input_[is1][is2][ideg]/(L*L); // (3/L^2)*c_{i-2}
                    c_poly_internal_[is1][is2][ideg+3] += c_poly_input_[is1][is2][ideg]/(L*L*L); // (1/L^3)*c_{i-3}
                }
                else if (C_damp_poly==2)
                {
                    c_poly_internal_[is1][is2][ideg] -= c_poly_input_[is1][is2][ideg]; // c_i
                    c_poly_internal_[is1][is2][ideg+1] += 2.0*c_poly_input_[is1][is2][ideg]/L; // (2/L)*c_{i-1}
                    c_poly_internal_[is1][is2][ideg+2] -= c_poly_input_[is1][is2][ideg]/(L*L); // (1/L^2)*c_{i-2}
                }
                else
                {
                    error_messages::inappropriate_argument("C_damp_poly (C in the damping function for the polynomial Jastrow)", C_damp_poly, "should be 2 or 3");
                }
            } // ideg
        } // is2
    } // is1
}

void Jastrow::set_derived_polynomial_parameters()
{
    assert(L_poly_.size()==2 && L_poly_[0].size()==2 && L_poly_[1].size()==2);
    assert(FourPI_power_L_.size()==2 && FourPI_power_L_[0].size()==2 && FourPI_power_L_[1].size()==2);
    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            FourPI_power_L_[is1][is2].resize(2*deg_poly_+6); // Fourpi*(L^i). +6 is sufficient both for C_damp_poly = 2 and 3.
            FourPI_power_L_[is1][is2][0] = FourPI;
            for (int ipoly=1; ipoly<FourPI_power_L_[is1][is2].size(); ipoly++)
            {
                FourPI_power_L_[is1][is2][ipoly] = L_poly_[is1][is2] * FourPI_power_L_[is1][is2][ipoly-1];
            }
        }
    }
}

// used only when RPA+poly Jastrow is used
void Jastrow::set_derived_RPA_polynomial_parameters(const CrystalStructure &crystal_structure,
                                                    const PlaneWaveBasis &plane_wave_basis)
{
    assert(F_long_.size()==2 && F_long_[0].size()==2 && F_long_[1].size()==2);
    assert(L_poly_.size()==2 && L_poly_[0].size()==2 && L_poly_[1].size()==2);
    assert(exp_LF_.size()==2 && exp_LF_[0].size()==2 && exp_LF_[1].size()==2);
    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            assert(std::abs(F_long_[is1][is2])>1e-8); // guaranteed in unnormalize_A_long()
            exp_LF_[is1][is2] = std::exp(-L_poly_[is1][is2] / F_long_[is1][is2]);
        }
    }

    // calculate M_integ[n] = L^{-n} \int_0^L exp(-r/F)*r^l dr (for J2_integ)
    M_integ_.resize(2);
    for (int is1=0; is1<2; is1++)
    {
        M_integ_[is1].resize(2);
        for (int is2=0; is2<2; is2++)
        {
            M_integ_[is1][is2].resize(max_deg_Taylor2_ + deg_poly_ + 1); // sufficient both for C_damp_poly = 2 and 3.
            M_integ_[is1][is2][0] = F_long_[is1][is2] * (1.0 - exp_LF_[is1][is2]);
            for (int ideg=1; ideg<max_deg_Taylor2_; ideg++)
            {
                M_integ_[is1][is2][ideg] = F_long_[is1][is2]
                    * (-exp_LF_[is1][is2] + ideg * M_integ_[is1][is2][ideg-1]/L_poly_[is1][is2]);
            }
        }
    }

    // Note: for deg_poly_=1, polynomial Jastrow = 0 (to satisfy the cuspless condition)
    if (deg_poly_>1) { prepare_for_J3_integ(crystal_structure, plane_wave_basis); } // for J3_integ 
}

void Jastrow::prepare_for_J3_integ(const CrystalStructure &crystal_structure,
                                   const PlaneWaveBasis &plane_wave_basis)
{
    // set Glist
    int num_G = plane_wave_basis.size_FFT_grid();
    Eigen::Vector3d Gvect;
    std::vector<double> Glist(num_G);
    for (int ipw=0; ipw<num_G; ipw++)
    {
        Gvect = crystal_structure.reciprocal_vectors().transpose()
            * plane_wave_basis.get_Gvector(ipw).cast<double>();
        Glist[ipw] = std::sqrt(Gvect.squaredNorm()); // |G|
    }

    // test error for various num_rmesh_J3_integ
    std::vector<std::vector<std::vector<double> > > integral(2);
    for (int is1=0; is1<2; is1++)
    {
        integral[is1].resize(2);
        for (int is2=0; is2<2; is2++)
        {
            integral[is1][is2].resize(num_G);
        }
    }
    std::vector<std::vector<std::vector<double> > > integral_new(2);
    for (int is1=0; is1<2; is1++)
    {
        integral_new[is1].resize(2);
        for (int is2=0; is2<2; is2++)
        {
            integral_new[is1][is2].resize(num_G);
        }
    }
    int num_rmesh_J3_integ = 101;
    test_J3_integ(num_rmesh_J3_integ, Glist, integral);
    bool is_error_small = false;
    for (int i=0; i<15; i++)
    {
        int num_rmesh_J3_integ_new = 2*num_rmesh_J3_integ - 1;
        test_J3_integ(num_rmesh_J3_integ_new, Glist, integral_new);

        double error = test_error_J3_integ(integral, integral_new);
        if (error < 1e-6) // accuracy of J3_integ
        {
            is_error_small = true;
            break;
        }
        else // for the next loop
        {
            num_rmesh_J3_integ = num_rmesh_J3_integ_new;
            integral = integral_new;
        }
    }
    if (!is_error_small) { error_messages::stop("Something wrong in prepare_for_J3_integ. Error is large."); }

    // set final variables
    MPI_Bcast(&num_rmesh_J3_integ, 1, MPI_INT, 0, MPI_COMM_WORLD);
    num_rmesh_J3_integ_ = num_rmesh_J3_integ;
    set_arrays_J3_integ(num_rmesh_J3_integ); // set rmesh_J3_integ_ & integrand_stored_J3_integ_
}
    
void Jastrow::test_J3_integ(const int num_rmesh_J3_integ, const std::vector<double> &Glist,
                            std::vector<std::vector<std::vector<double> > > &integral)
{
    int num_G = Glist.size();
    assert(integral.size() == 2);
    for (int is1=0; is1<2; is1++)
    {
        assert(integral[is1].size() == 2);
        for (int is2=0; is2<2; is2++)
        {
            assert(integral[is1][is2].size() == num_G);
        }
    }

    // set integrand_stored_J3_integ_ and rmesh_J3_integ_
    set_arrays_J3_integ(num_rmesh_J3_integ);

    // calculate J3_integ
    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            for (int ig=0; ig<num_G; ig++)
            {
                integral[is1][is2][ig] = J3_integ(Glist[ig], is1, is2);
            }
        } // is2
    } // is1
}

double Jastrow::test_error_J3_integ(const std::vector<std::vector<std::vector<double> > > &integral1,
                                    const std::vector<std::vector<std::vector<double> > > &integral2) const
{
    assert(integral1.size() == 2 && integral2.size() == 2);
    for (int is1=0; is1<2; is1++)
    {
        assert(integral1[is1].size() == 2 && integral2[is1].size() == 2);
    }
    const int num_G = integral1[0][0].size();
    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            assert(integral1[is1][is2].size() == num_G);
            assert(integral2[is1][is2].size() == num_G);
        }
    }

    double error_max = 0.0;
    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            for (int ig=0; ig<num_G; ig++)
            {
                double error_at_ig
                    = std::abs(integral1[is1][is2][ig] - integral2[is1][is2][ig]);
                if (error_max < error_at_ig) { error_max = error_at_ig; }
            }
        }
    }
    return error_max;
}
    
void Jastrow::set_arrays_J3_integ(const int num_rmesh_J3_integ)
{
    assert(F_long_.size()==2 && F_long_[0].size()==2 && F_long_[1].size()==2);
    assert(L_poly_.size()==2 && L_poly_[0].size()==2 && L_poly_[1].size()==2);
    assert(num_rmesh_J3_integ > 1);
    assert(num_rmesh_J3_integ%2==1);

    std::vector<std::vector<Eigen::VectorXd> > coeff(2);
    for (int is1=0; is1<2; is1++)
    {
        coeff[is1].resize(2);
        for (int is2=0; is2<2; is2++)
        {
            coeff[is1][is2].resize(num_rmesh_J3_integ);
        }
    }
    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            double dr = L_poly_[is1][is2] / F_long_[is1][is2]
                / static_cast<double>(num_rmesh_J3_integ - 1);
            for (int ir=0; ir<num_rmesh_J3_integ; ir++)
            {
                // coefficient for the composite Simpson's rule
                if (ir==0 || ir==num_rmesh_J3_integ-1)
                {
                    coeff[is1][is2](ir) = 1.0;
                }
                else if (ir%2==0)
                {
                    coeff[is1][is2](ir) = 2.0;
                }
                else
                {
                    coeff[is1][is2](ir) = 4.0;
                }
            }
            coeff[is1][is2] *= (dr/3.0)*FourPI;
        }
    }

    rmesh_J3_integ_.resize(2);
    integrand_stored_J3_integ_.resize(2);
    for (int is1=0; is1<2; is1++)
    {
        rmesh_J3_integ_[is1].resize(2);
        integrand_stored_J3_integ_[is1].resize(2);
        for (int is2=0; is2<2; is2++)
        {
            double dr = L_poly_[is1][is2] / F_long_[is1][is2]
                / static_cast<double>(num_rmesh_J3_integ - 1);

            rmesh_J3_integ_[is1][is2].resize(num_rmesh_J3_integ);
            for (int ir=0; ir<num_rmesh_J3_integ; ir++)
            {
                rmesh_J3_integ_[is1][is2](ir) = ir*dr;
            }

            integrand_stored_J3_integ_[is1][is2].resize(num_rmesh_J3_integ);
            integrand_stored_J3_integ_[is1][is2](0) = 0.0; // r = 0
            for (int ir=1; ir<num_rmesh_J3_integ; ir++)
            {
                double expr = std::exp(-rmesh_J3_integ_[is1][is2][ir]);
                integrand_stored_J3_integ_[is1][is2](ir) = 
                    (expr - 1.0)/rmesh_J3_integ_[is1][is2](ir) + expr;
            }
            integrand_stored_J3_integ_[is1][is2]
                = integrand_stored_J3_integ_[is1][is2].array() * coeff[is1][is2].array();
        } // is2
    } // is1
}

void Jastrow::bcast_polynomial_parameters(const bool am_i_mpi_rank0)
{
    // deg_poly_
    MPI_Bcast(&deg_poly_, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // C_damp_poly_
    MPI_Bcast(&C_damp_poly_, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // L_poly_
    assert(L_poly_.size()==2 && L_poly_[0].size()==2 && L_poly_[1].size()==2);
    MPI_Bcast(&L_poly_[0][0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&L_poly_[1][0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // c_poly_input/internal_
    if (am_i_mpi_rank0)
    {
        assert(c_poly_input_.size()==2 && c_poly_input_[0].size()==2 && c_poly_input_[1].size()==2);
        assert(c_poly_internal_.size()==2 && c_poly_internal_[0].size()==2 && c_poly_internal_[1].size()==2);
        for (int is1=0; is1<2; is1++)
        {
            for (int is2=0; is2<2; is2++)
            {
                assert(c_poly_input_[is1][is2].size()==deg_poly_ && c_poly_internal_[is1][is2].size()==deg_poly_+C_damp_poly_);
            }
        }
    }
    else
    {
        c_poly_input_.resize(2);
        c_poly_internal_.resize(2);
        for (int is1=0; is1<2; is1++)
        {
            c_poly_input_[is1].resize(2);
            c_poly_internal_[is1].resize(2);
            for (int is2=0; is2<2; is2++)
            {
                c_poly_input_[is1][is2].resize(deg_poly_);
                c_poly_internal_[is1][is2].resize(deg_poly_+C_damp_poly_);
            }
        }
    }
    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            MPI_Bcast(&c_poly_input_[is1][is2][0], deg_poly_, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&c_poly_internal_[is1][is2][0], deg_poly_+C_damp_poly_, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
}

void Jastrow::print_polynomial_parameters(std::ostream *ost) const
{
    *ost << "  Polynomial Jastrow parameters in atomic unit" << std::endl;
    *ost << "   e-e basis: natural power with deg_poly = " << deg_poly_ << std::endl;
    for (int is1=0; is1<2; is1++)
    {
        std::string is1_name = is1==0 ? "up" : "down";
        for (int is2=0; is2<2; is2++)
        {
            std::string is2_name = is2==0 ? "up" : "down";

            *ost << "   " << is1_name << "-" << is2_name << " channel:" << std::endl;
            *ost << "    L_poly = " << L_poly_[is1][is2] << std::endl;
            *ost << "     for a polynomial cutoff function with C = " << C_damp_poly_ << ": (1-r/L)^" << C_damp_poly_ << "* step_function(L-r)" << std::endl;

            for (int ipoly=0; ipoly<deg_poly_; ipoly++)
            {
                *ost << "    c_poly (r^" << ipoly << " term = c_" << ipoly+1 << ") = " << c_poly_input_[is1][is2][ipoly] << std::endl;
            }
            *ost << "     Note: these coefficients are defined for J = exp[+u], u = sum c_poly * r^l * cutoff function (see above)" << std::endl;
        }
    }
    *ost << std::endl;
}

void Jastrow::print_RPA_polynomial_parameters(std::ostream *ost) const
{
    *ost << "  RPA + Polynomial Jastrow parameters" << std::endl;
    *ost << "   num_rmesh_J3_integ = " << num_rmesh_J3_integ_ << std::endl;
    *ost << std::endl;
}
