// [class Potentials]
// pseudopotential, Jastrow, Ewald sum

// given in atomic unit (Hartree)

#include "include/header.hpp"

// calculate "derivative_pseudo_local"
void Potentials::differentiate_local_potential(const CrystalStructure &crystal_structure,
                                               PlaneWaveBasis &plane_wave_basis,
                                               std::vector<std::vector<Eigen::VectorXcd> > &derivative_pseudo_local) const
{
    const int npw = plane_wave_basis.size_FFT_grid();
    const int num_atoms = crystal_structure.num_atoms();
    const int num_atomic_species = crystal_structure.num_atomic_species();
    const double unit_cell_volume = crystal_structure.unit_cell_volume();

    // initialize
    derivative_pseudo_local.resize(num_atoms);
    for (int iatom=0; iatom<num_atoms; iatom++)
    {
        derivative_pseudo_local[iatom].resize(3);
        for (int idim=0; idim<3; idim++) 
        {
            derivative_pseudo_local[iatom][idim].resize(npw);
        }
    }

    Eigen::MatrixXd Gvector(npw, 3);
    for (int ipw=0; ipw<npw; ipw++)
    {
        // Gvector[ipw] = i1*b1 + i2*b2 + i3*b3. (ipw = (i1, i2, i3))
        Gvector.row(ipw) =
            ((plane_wave_basis.get_Gvector(ipw)).cast<double>()).transpose()
            * crystal_structure.reciprocal_vectors();
    }

    // set Vpp_local_atom_G (Vpp_local_G for each atom) (G = in G space)
    std::vector<Eigen::VectorXcd> Vpp_local_atom_G(num_atomic_species);
    for (int iatomic_species=0; iatomic_species<num_atomic_species; iatomic_species++)
    {
        set_Vpp_local_atom(crystal_structure, plane_wave_basis, Vpp_local_atom_G[iatomic_species], iatomic_species);
    }

    for (int iatom=0; iatom<num_atoms; iatom++)
    {
        int iatomic_species = crystal_structure.index_of_atoms()[iatom];
        Eigen::VectorXd GdotR = Gvector * crystal_structure.atomic_position_cartesian()[iatom]; // array(npw)

        for (int ipw=0; ipw<npw; ipw++)
        {
            Complex prefactor = (-I*std::cos(GdotR(ipw)) - std::sin(GdotR(ipw))) // (-I)*exp(-I*GdotR)
                * Vpp_local_atom_G[iatomic_species](ipw);

            for (int idim=0; idim<3; idim++)
            {
                derivative_pseudo_local[iatom][idim](ipw) = prefactor * Gvector(ipw, idim);
            }
        } // ipw
        for (int idim=0; idim<3; idim++)
        {
            plane_wave_basis.FFT_backward(derivative_pseudo_local[iatom][idim],
                                          derivative_pseudo_local[iatom][idim]);
        }
    } // iatom
}
