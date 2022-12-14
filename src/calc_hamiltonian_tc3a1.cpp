// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

// *** tc3b1 ***
// calculate -0.5*sum_q1 sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 |j,q1,q2>
void calc_hamiltonian::tc3a1(const Parallelization &parallelization, 
                             const Method &method,
                             const CrystalStructure &crystal_structure,
                             const Potentials &potentials,
                             const Spin &spin, const Kpoints &kpoints,
                             PlaneWaveBasis &plane_wave_basis, 
                             const BlochStates &bloch_states, 
                             const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
                             std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H3phi,
                             std::ostream *ost)
{
    const int num_independent_spins = spin.num_independent_spins();
    const int num_irreducible_kpoints_scf = kpoints.num_irreducible_kpoints_scf();
    const int num_irreducible_kpoints = method.calc_mode()=="SCF" ? 
        kpoints.num_irreducible_kpoints_scf() : kpoints.num_irreducible_kpoints_band();
    const int num_spinor = (spin.is_spinor()==false) ? 1 : 2;

    // (-0.5) * I^2 = +0.5
    const double three_body_factor = 0.5/(kpoints.num_kpoints()*kpoints.num_kpoints()*
                                          crystal_structure.unit_cell_volume()*
                                          crystal_structure.unit_cell_volume());

    std::vector<Eigen::Vector3d> Gvect(plane_wave_basis.size_FFT_grid());
    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
    {
        Gvect[ipw] = crystal_structure.reciprocal_vectors().transpose()
            * plane_wave_basis.get_Gvector(ipw).cast<double>();
    }

    Eigen::VectorXcd orbital(plane_wave_basis.size_FFT_grid());
    std::vector<std::vector<Eigen::VectorXcd> > dnu(num_independent_spins); // spin index of x1
    for (int ispin=0; ispin<num_independent_spins; ispin++) 
    {
        dnu[ispin].resize(3);
        for (int idim=0; idim<3; idim++)
        {
            dnu[ispin][idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
        }
    }
    std::vector<Eigen::VectorXcd> V_3body(num_independent_spins); // spin index of x1
    for (int ispin=0; ispin<num_independent_spins; ispin++) 
    {
        V_3body[ispin] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }

    // set dnu = \sum_q <*,q| \nabla_1 u_12 | *,q>
    for (int ispin=0; ispin<num_independent_spins; ispin++) 
    {
        for (int jspin=0; jspin<2; jspin++)
        {
            orbital = bloch_states.density()[jspin]; // density(R). Name "orbital" is meaningless here
            plane_wave_basis.FFT_forward(orbital, orbital); // -> density(G)

            for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
            {
                double uk = potentials.jastrow.uk(Gvect[ipw], ispin, jspin);
                for (int idim=0; idim<3; idim++)
                {
                    dnu[ispin][idim](ipw) +=
                        uk * Gvect[ipw](idim) * orbital(ipw);
                }
            }
        }
        for (int idim=0; idim<3; idim++)
        {
            plane_wave_basis.FFT_backward(dnu[ispin][idim], dnu[ispin][idim]); // -> dnu(R)
        }
    }

    // set V_3body = \sum_q1 \sum_q2 <*,q1,q2| \nabla_1 u_12 \nabla_1 u_13 | *,q1,q2>
    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int idim=0; idim<3; idim++)
        {
            V_3body[ispin] = V_3body[ispin].array() + dnu[ispin][idim].array().square(); // dnu*dnu
        }
    }

    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            const int num_G_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.num_G_at_k_scf()[ik] : plane_wave_basis.num_G_at_k_band()[ik];
            const Eigen::VectorXi Gindex_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.Gindex_at_k_scf()[ispin][ik][0] :  plane_wave_basis.Gindex_at_k_band()[ispin][ik][0];

            for (int jband=0; jband<H3phi[ispin][ik].size(); jband++)
            {
                if (!parallelization.is_assigned_irreducible_kpoints_all_bands()[ispin][ik][jband]) { continue; }
                for (int jspinor=0; jspinor<num_spinor; jspinor++)
                {
                    // phi -> orbital on the FFT-grid
                    plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, // isym = 0
                                                         false, // time-rersal not used for isym=0
                                                         phi[ispin][ik][jband][jspinor], orbital,
                                                         method.calc_mode());
                    plane_wave_basis.FFT_backward(orbital, orbital); // -> phi_j(R)
                    orbital = orbital.array() * V_3body[ispin].array();
                    plane_wave_basis.FFT_forward(orbital, orbital);
                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                    {
                        H3phi[ispin][ik][jband][jspinor](ipw_at_k) +=
                            three_body_factor * orbital(Gindex_at_k(ipw_at_k));
                    }
                } // jspinor
            } // jband
        } // ik
    } // ispin
}
