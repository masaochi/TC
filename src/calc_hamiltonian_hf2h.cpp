// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

// two-body, Hartree term (Hartree-Fock)
void calc_hamiltonian::hf2h(const Parallelization &parallelization, 
                            const Method &method,
                            const CrystalStructure &crystal_structure,
                            const Potentials &potentials,
                            const Spin &spin, const Kpoints &kpoints,
                            PlaneWaveBasis &plane_wave_basis, 
                            const BlochStates &bloch_states, 
                            const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
                            std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H2phi,
                            std::ostream *ost)
{
    const int num_independent_spins = spin.num_independent_spins();
    const int num_irreducible_kpoints = method.calc_mode()=="SCF" ? 
        kpoints.num_irreducible_kpoints_scf() : kpoints.num_irreducible_kpoints_band();
    const int num_spinor = (spin.is_spinor()==false) ? 1 : 2;

    const double two_body_factor = 1.0/(kpoints.num_kpoints()*crystal_structure.unit_cell_volume());

    Eigen::VectorXcd orbital(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd V_orbital(plane_wave_basis.size_FFT_grid());
    Eigen::VectorXcd V_hartree(plane_wave_basis.size_FFT_grid());

    // set Hartree potential
    V_hartree = bloch_states.density()[0] + bloch_states.density()[1]; // density(R)
    plane_wave_basis.FFT_forward(V_hartree, V_hartree); // density(G)
    V_hartree(0) = 0.0; // G=0
    for (int ipw=1; ipw<plane_wave_basis.size_FFT_grid(); ipw++) // Note that ipw=0 is excluded
    {
        Eigen::Vector3d Gvect = 
            crystal_structure.reciprocal_vectors().transpose()
            * plane_wave_basis.get_Gvector(ipw).cast<double>();

        V_hartree(ipw) *= FourPI/Gvect.squaredNorm();
    }
    plane_wave_basis.FFT_backward(V_hartree, V_hartree);

// debug begin
//    Complex tmp_debg = 0.0;
//    for (int ipw=0; ipw<plane_wave_basis.size_FFT_grid(); ipw++) 
//    {
//        tmp_debg += V_hartree(ipw) * (bloch_states.density[0](ipw) + bloch_states.density[1](ipw));
//    }
//    // 2 = two-body term
//    tmp_debg *= two_body_factor/static_cast<double>(2*plane_wave_basis.size_FFT_grid()*kpoints.num_kpoints());
//    std::cout << "Debug Hartree: " << tmp_debg << std::endl;
// debug end

    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            const int num_G_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.num_G_at_k_scf()[ik] : plane_wave_basis.num_G_at_k_band()[ik];
            const Eigen::VectorXi Gindex_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.Gindex_at_k_scf()[ispin][ik][0] :  plane_wave_basis.Gindex_at_k_band()[ispin][ik][0];

            for (int jband=0; jband<H2phi[ispin][ik].size(); jband++)
            {
                if (!parallelization.is_assigned_irreducible_kpoints_all_bands()[ispin][ik][jband]) { continue; }
                for (int jspinor=0; jspinor<num_spinor; jspinor++)
                {
                    // phi -> orbital on the FFT-grid
                    plane_wave_basis.get_orbital_FFTgrid(ispin, ik, 0, // isym = 0
                                                         false, // time-rersal not used for isym=0
                                                         phi[ispin][ik][jband][jspinor], orbital,
                                                         method.calc_mode());
                    plane_wave_basis.FFT_backward(orbital, orbital);
                    V_orbital = V_hartree.array() * orbital.array();
                    plane_wave_basis.FFT_forward(V_orbital, V_orbital);
                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++)
                    {
                        H2phi[ispin][ik][jband][jspinor](ipw_at_k)
                            += two_body_factor * V_orbital(Gindex_at_k(ipw_at_k));
                    }
                } // jspinor
            } // jband
        } // ik
    } // ispin
}
