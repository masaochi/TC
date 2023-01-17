// [class Jastrow]
// jastrow parameters. "Jastrow jastrow" is a member of "Potential potential".

#include "include/header.hpp"

// private

// A_long is at present given in the unit of sqrt(V/(4piN)).
// An appropriate factor is later multiplied in unnormalize_A_long() to make it in Hartree unit.
// F_long_ is also initialized there.
void Jastrow::set_A_long_normalized(const std::vector<std::vector<double> > &A_long)
{
    assert(A_long_.size()==2 && A_long_[0].size()==2 && A_long_[1].size()==2);
    assert(std::abs(A_long[0][1] - A_long[1][0])<1e-8); // A_long_up_dn = A_long_dn_up
    A_long_ = A_long;
}

void Jastrow::print_RPAjastrow(std::ostream *ost)
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
}

// public

void Jastrow::set(const std::vector<std::vector<double> > &A_long)
{
    set_A_long_normalized(A_long);
}

// In "input.in", multiply "sqrt(V/(4piN))" (V:volume, N:no. of valence electrons) with the A parameter
// to make it in Hartree unit. F_long is also determined so as to satisfy the cusp condition.
void Jastrow::unnormalize_A_long(const CrystalStructure &crystal_structure,
                                 const BlochStates &bloch_states,
                                 const Spin &spin,
                                 std::ostream *ost)
{
    assert(A_long_.size()==2 && A_long_[0].size()==2 && A_long_[1].size()==2);
    assert(F_long_.size()==2 && F_long_[0].size()==2 && F_long_[1].size()==2);
    if (spin.num_independent_spins()==1 && !spin.is_spinor()) // no-spin, up & down are equivalent
    {
        if (std::abs(A_long_[0][0] - A_long_[1][1])>1e-8)
        {
            error_messages::stop("This calculation assumes that up & down spins are equivalent, so please specify the same value for A_up_up & A_dn_dn");
        }
    }

    const double volume = crystal_structure.unit_cell_volume();
    assert(volume > 1e-8);
    const double num_electrons = bloch_states.num_electrons(); // including the spin degree of freedom
    assert(num_electrons > 1e-8);

    for (int is1=0; is1<2; is1++)
    {
        for (int is2=0; is2<2; is2++)
        {
            A_long_[is1][is2] *= std::sqrt(volume/(FourPI*num_electrons));

            // F_long_ is set so as to satisfy the cusp condition
            // If A=0 then F=1 (default) is used.
            if (std::abs(A_long_[is1][is2])>1e-8) 
            {
                if (is1==is2) { F_long_[is1][is2] = std::sqrt(2*A_long_[is1][is2]); }
                if (is1!=is2) { F_long_[is1][is2] = std::sqrt(A_long_[is1][is2]); }
            }
        } // is2
    } // is1
    print_RPAjastrow(ost);
}

void Jastrow::bcast()
{
    assert(A_long_.size()==2 && A_long_[0].size()==2 && A_long_[1].size()==2);
    assert(F_long_.size()==2 && F_long_[0].size()==2 && F_long_[1].size()==2);
    MPI_Bcast(&A_long_[0][0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&A_long_[1][0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&F_long_[0][0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&F_long_[1][0], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // 4*pi*A, F*F, 1/(F*F) etc. are set here
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
            assert(std::abs(F_long_[is1][is2])>1e-8); // guaranteed in normalize_A_long()

            FourPI_A_long_[is1][is2] = FourPI*A_long_[is1][is2];
            F2_long_[is1][is2] = F_long_[is1][is2] * F_long_[is1][is2];
            Finv2_long_[is1][is2] = 1.0/F2_long_[is1][is2];
            PI_AAFinv2_long_[is1][is2] = PI*A_long_[is1][is2]*A_long_[is1][is2]*Finv2_long_[is1][is2];
        }
    }
}
