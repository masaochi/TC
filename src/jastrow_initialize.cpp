// [class Jastrow]
// jastrow parameters. "Jastrow jastrow" is a member of "Potential potential".

#include "include/header.hpp"

// A_long is at present given in the unit of sqrt(V/(4piN)).
// An appropriate factor is later multiplied in unnormalize_A_long() to make it in Hartree unit.
// F_long_ is also initialized there.
void Jastrow::set_A_long(const std::vector<std::vector<double> > &A_long)
{
    assert(A_long.size()==2 && A_long[0].size()==2 && A_long[1].size()==2);
    assert(std::abs(A_long[0][1] - A_long[1][0])<1e-8); // A_long_up_dn = A_long_dn_up

    A_long_ = A_long;

    is_A_zero_ = true;
    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            if (std::abs(A_long_[is1][is2])>1e-8) { is_A_zero_ = false; }
        }
    }
}

// In "input.in", multiply "sqrt(V/(4piN))" (V:volume, N:no. of valence electrons) with the A parameter
// to make it in Hartree unit. F_long is also determined so as to satisfy the cusp condition.
void Jastrow::unnormalize_A_long(const double &volume,
                                 const double &num_electrons, // incl. spin degeneracy
                                 const Spin &spin)
{
    assert(A_long_.size()==2 && A_long_[0].size()==2 && A_long_[1].size()==2);
    if (spin.num_independent_spins()==1 && !spin.is_spinor()) // no-spin, up & down are equivalent
    {
        if (std::abs(A_long_[0][0] - A_long_[1][1])>1e-8)
        {
            error_messages::stop("This calculation assumes that up & down spins are equivalent. Please specify the same value for A_up_up & A_dn_dn");
        }
    }
    
    if (volume < 1e-8) { error_messages::stop("Non-positive volume of the unit cell was detected in potentials.jastrow.unnormalize_A_long()."); }
    if (num_electrons < 1e-8) { error_messages::stop("Non-positive num. of electrons was detected in potentials.jastrow.unnormalize_A_long()."); }
    
    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            A_long_[is1][is2] *= std::sqrt(volume/(FourPI*num_electrons));
        } // is2
    } // is1
}

void Jastrow::set_F_long(const std::vector<std::vector<double> > &F_long)
{
    assert(F_long.size()==2 && F_long[0].size()==2 && F_long[1].size()==2);
    assert(std::abs(F_long[0][1] - F_long[1][0])<1e-8); // F_long_up_dn = F_long_dn_up

    F_long_ = F_long;
}

void Jastrow::impose_cusp(const bool A_to_F)
{
    assert(A_long_.size()==2 && A_long_[0].size()==2 && A_long_[1].size()==2);
    assert(F_long_.size()==2 && F_long_[0].size()==2 && F_long_[1].size()==2);

    if (A_to_F) // set F using A so that the cusp condition is satisfied
    {
        if (is_A_zero_)
        {
            set_F_long({{1.0, 1.0}, {1.0, 1.0}}); // since we calculate 1/F in set_derived_RPA_parameters()
        }
        else
        {
            set_F_long({{std::sqrt(2*A_long_[0][0]), std::sqrt(A_long_[0][1])},
                    {std::sqrt(A_long_[1][0]), std::sqrt(2*A_long_[1][1])}});
        }
    }
    else // set A using F so that the cusp condition is satisfied
    {
        set_A_long({{F_long_[0][0]*F_long_[0][0]/2.0, F_long_[0][1]*F_long_[0][1]},
                {F_long_[1][0]*F_long_[1][0], F_long_[1][1]*F_long_[1][1]/2.0}});
    }
}

// set 4*pi*A, F*F, 1/(F*F) etc.
void Jastrow::set_derived_RPA_parameters()
{
    assert(A_long_.size()==2 && A_long_[0].size()==2 && A_long_[1].size()==2);
    assert(F_long_.size()==2 && F_long_[0].size()==2 && F_long_[1].size()==2);

    assert(FourPI_A_long_.size()==2 && FourPI_A_long_[0].size()==2 
           && FourPI_A_long_[1].size()==2);
    assert(F2_long_.size()==2 && F2_long_[0].size()==2 && F2_long_[1].size()==2);
    assert(Finv2_long_.size()==2 && Finv2_long_[0].size()==2 && Finv2_long_[1].size()==2);
    assert(PI_AAFinv2_long_.size()==2 && PI_AAFinv2_long_[0].size()==2 
           && PI_AAFinv2_long_[1].size()==2);
    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            assert(std::abs(F_long_[is1][is2])>1e-8); // guaranteed in unnormalize_A_long()

            FourPI_A_long_[is1][is2] = FourPI*A_long_[is1][is2];
            Finv_long_[is1][is2] = 1.0/F_long_[is1][is2];
            F2_long_[is1][is2] = F_long_[is1][is2] * F_long_[is1][is2];
            Finv2_long_[is1][is2] = 1.0/F2_long_[is1][is2];
            PI_AAFinv2_long_[is1][is2] = PI*A_long_[is1][is2]*A_long_[is1][is2]*Finv2_long_[is1][is2];
        }
    }
}

void Jastrow::set_is_A_given_in_input_in(const bool is_A_given_in_input_in)
{
    is_A_given_in_input_in_ = is_A_given_in_input_in;
}

void Jastrow::bcast_A_long()
{
    assert(A_long_.size()==2 && A_long_[0].size()==2 && A_long_[1].size()==2);

    MPI_Bcast(&A_long_[0][0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&A_long_[1][0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&is_A_zero_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&is_A_given_in_input_in_, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
}

void Jastrow::bcast_F_long()
{
    assert(F_long_.size()==2 && F_long_[0].size()==2 && F_long_[1].size()==2);

    MPI_Bcast(&F_long_[0][0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&F_long_[1][0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void Jastrow::print_RPA_parameters(std::ostream *ost) const
{
    *ost << "  RPA-type Jastrow parameters in atomic unit" << std::endl;
    *ost << "   A_up_up = " << A_long_[0][0] << std::endl;
    *ost << "   F_up_up = " << F_long_[0][0] << std::endl;
    *ost << "   A_up_dn = A_dn_up = " << A_long_[0][1] << std::endl;
    assert(std::abs(A_long_[0][1] - A_long_[1][0])<1e-8);
    *ost << "   F_up_dn = F_dn_up = " << F_long_[0][1] << std::endl;
    assert(std::abs(F_long_[0][1] - F_long_[1][0])<1e-8);
    *ost << "   A_dn_dn = " << A_long_[1][1] << std::endl;
    *ost << "   F_dn_dn = " << F_long_[1][1] << std::endl;
    *ost << std::endl;
}
