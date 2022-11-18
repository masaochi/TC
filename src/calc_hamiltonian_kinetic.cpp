// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#include "include/header.hpp"

void calc_hamiltonian::kinetic(const Parallelization &parallelization, 
                               const Method &method,
                               const CrystalStructure &crystal_structure,
                               const Spin &spin, const Kpoints &kpoints,
                               PlaneWaveBasis &plane_wave_basis, 
                               const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
                               std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H1phi,
                               std::ostream *ost)
{
    const int num_independent_spins = spin.num_independent_spins();
    const int num_irreducible_kpoints = method.calc_mode()=="SCF" ? 
        kpoints.num_irreducible_kpoints_scf() : kpoints.num_irreducible_kpoints_band();
    const int num_spinor = (spin.is_spinor()==false) ? 1 : 2;

    Eigen::Vector3d kGvect;

    for (int ispin=0; ispin<num_independent_spins; ispin++)
    {
        for (int ik=0; ik<num_irreducible_kpoints; ik++)
        {
            const Eigen::Vector3d kvector = method.calc_mode()=="SCF" ? 
                kpoints.kvectors_scf()[ik][0] : kpoints.kvectors_band()[ik][0];
            const int num_G_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.num_G_at_k_scf()[ik] : plane_wave_basis.num_G_at_k_band()[ik];
            const Eigen::VectorXi Gindex_at_k = method.calc_mode()=="SCF" ?
                plane_wave_basis.Gindex_at_k_scf()[ispin][ik][0] :  plane_wave_basis.Gindex_at_k_band()[ispin][ik][0];

            for (int jband=0; jband<H1phi[ispin][ik].size(); jband++)
            {
                if (!parallelization.is_assigned_irreducible_kpoints_all_bands()[ispin][ik][jband]) { continue; }
                for (int jspinor=0; jspinor<num_spinor; jspinor++)
                {
                    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++) 
                    {
                        kGvect = crystal_structure.reciprocal_vectors().transpose()
                            *(kvector + plane_wave_basis.get_Gvector(Gindex_at_k(ipw_at_k)).cast<double>());
                        
                        // -0.5*I^2 = +0.5
                        H1phi[ispin][ik][jband][jspinor](ipw_at_k)
                            += 0.5 * kGvect.squaredNorm() * phi[ispin][ik][jband][jspinor](ipw_at_k);
                    }

                } // jspinor
            } // jband
        } // ik
    } // ispin
}
