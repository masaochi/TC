// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

// *** tc3b1 ***
// calculate -0.5*2*sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_2 u_23 |j,q1,q2>
void calc_hamiltonian::tc3b1(const Parallelization &parallelization, 
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

    // (-0.5*2) * I^2 * (-1) (convolution \nabla_2 u_21 = -\nabla_1 u_12) = -1.0
    const double three_body_factor = -1.0/(kpoints.num_kpoints()*kpoints.num_kpoints()*
                                           crystal_structure.unit_cell_volume()*
                                           crystal_structure.unit_cell_volume());

    std::vector<Eigen::Vector3d> Gvect(plane_wave_basis.size_FFT_grid());
    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
    {
        Gvect[ipw] = crystal_structure.reciprocal_vectors().transpose()
            * plane_wave_basis.get_Gvector(ipw).cast<double>();
    }

    Eigen::VectorXcd orbital(plane_wave_basis.size_FFT_grid());
    std::vector<std::vector<Eigen::VectorXcd> > dnu(num_independent_spins); // spin index of x2
    for (int ispin2=0; ispin2<num_independent_spins; ispin2++) 
    {
        dnu[ispin2].resize(3);
        for (int idim=0; idim<3; idim++)
        {
            dnu[ispin2][idim] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
        }
    }
    std::vector<Eigen::VectorXcd> V_3body(num_independent_spins); // spin index of x1
    for (int ispin1=0; ispin1<num_independent_spins; ispin1++) 
    {
        V_3body[ispin1] = Eigen::VectorXcd::Zero(plane_wave_basis.size_FFT_grid());
    }

    for (int ispin2=0; ispin2<num_independent_spins; ispin2++) // spin index of x2
    {
        for (int ispin3=0; ispin3<2; ispin3++) // spin index of x3
        {
            orbital = bloch_states.density()[ispin3]; // density(R). Name "orbital" is meaningless here
            plane_wave_basis.FFT_forward(orbital, orbital); // -> density(G)

            // set \sum_q2 <*,*,q2| \nabla_2 u_23 |*,*,q2>
            for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
            {
                double uk = potentials.jastrow.uk(Gvect[ipw], ispin2, ispin3);
                for (int idim=0; idim<3; idim++)
                {
                    dnu[ispin2][idim](ipw) +=
                        uk * Gvect[ipw](idim) * orbital(ipw);
                }
            }
        }
        for (int idim=0; idim<3; idim++)
        {
            plane_wave_basis.FFT_backward(dnu[ispin2][idim], dnu[ispin2][idim]); // -> dnu(R)
        }

        // set sum_q1 sum_q2 <*,q1,q2| \nabla_2 u_23 |*,q1,q2>
        orbital = bloch_states.density()[ispin2]; // density(R). Name "orbital" is meaningless here
        for (int idim=0; idim<3; idim++)
        {
            dnu[ispin2][idim] = dnu[ispin2][idim].array() * orbital.array();
            plane_wave_basis.FFT_forward(dnu[ispin2][idim], dnu[ispin2][idim]); // -> (G)
        }
    }

    // set V_3body = sum_q1 \sum_q2 <*,q1,q2| \nabla_2 u_21 \nabla_1 u_23 | *,q1,q2>
    for (int ispin1=0; ispin1<num_independent_spins; ispin1++) // spin index of x1
    {
        for (int ispin2=0; ispin2<2; ispin2++) // spin index of x2
        {
            int ispin2_ref = num_independent_spins==1 ? 0 : ispin2;
            for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++)
            {
                double uk = potentials.jastrow.uk(Gvect[ipw], ispin1, ispin2);
                for (int idim=0; idim<3; idim++)
                {
                    V_3body[ispin1](ipw) +=
                        uk * Gvect[ipw](idim) * dnu[ispin2_ref][idim](ipw);
                }
            }
        }
        plane_wave_basis.FFT_backward(V_3body[ispin1], V_3body[ispin1]); // -> (R)
    }

    for (int ispin=0; ispin<num_independent_spins; ispin++) // spin index of x1
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
                                                         jband, false, // time-rersal not used for isym=0
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
