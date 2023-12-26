// [class TotalEnergy]
// total energy, ewald energy, etc.

#include "include/header.hpp"

// private

double TotalEnergy::calculate_ewald_energy(const CrystalStructure &crystal_structure) const
{
    double return_val = 0.0;

    const int num_atoms = crystal_structure.num_atoms();
    const int num_atomic_species = crystal_structure.num_atomic_species();
    const std::vector<double> &Z_valence = crystal_structure.Z_valence();
    const double unit_cell_volume = crystal_structure.unit_cell_volume();
    const Eigen::Matrix3d lattice_vectors = crystal_structure.lattice_vectors();
    const Eigen::Matrix3d reciprocal_vectors = crystal_structure.reciprocal_vectors();

    const double param_ewald = PI/std::cbrt(unit_cell_volume); // cbrt = cubic root
    const double factor_param_ewald = 1.0/(4*param_ewald*param_ewald);
    const double threshold_ewald = 1e-8;
    double temp_ewald, temp_ewald_sum; // temporary variable for Ewald sum

    // single i-sum terms
    std::vector<double> term3_Gsum(num_atomic_species);
    std::vector<double> term4_Rsum(num_atomic_species);
    for (int iatomic_species=0; iatomic_species<num_atomic_species; iatomic_species++)
    {
        // Term 3 (except the Z^2 factor)
        term3_Gsum[iatomic_species] = 0.0;
        for (int nshell=1; ; nshell++) // Note: G=0 is excluded
        {
            if (nshell%2==1) { temp_ewald = 0.0; } // odd layer
            for (int Gx=-nshell; Gx<=nshell; Gx++)
            {
                for (int Gy=-nshell; Gy<=nshell; Gy++)
                {
                    for (int Gz=-nshell; Gz<=nshell; Gz++)
                    {
                        if (Gx!=nshell && Gy!=nshell && Gz!=nshell
                            && -Gx!=nshell && -Gy!=nshell && -Gz!=nshell) { continue; }
                        Eigen::Vector3i tmpG = {Gx, Gy, Gz};
                        Eigen::Vector3d Gvector = reciprocal_vectors.transpose() * (tmpG.cast<double>());
                        double G2 = Gvector.squaredNorm();
                        temp_ewald += std::exp(-G2*factor_param_ewald)/G2;
                    } // Gz
                } // Gy
            } // Gx
            if (nshell%2==0)
            {
                term3_Gsum[iatomic_species] += temp_ewald;
                if (std::abs(temp_ewald) < threshold_ewald) { break; }
            }
        } // nshell
        term3_Gsum[iatomic_species] *= 2*PI/unit_cell_volume;

        // Term 4 (except the Z^2 factor)
        term4_Rsum[iatomic_species] = 0.0;
        for (int nshell=1; ; nshell++) // Note: R=0 is excluded
        {
            if (nshell%2==1) { temp_ewald = 0.0; } // odd layer
            for (int Rx=-nshell; Rx<=nshell; Rx++)
            {
                for (int Ry=-nshell; Ry<=nshell; Ry++)
                {
                    for (int Rz=-nshell; Rz<=nshell; Rz++)
                    {
                        if (Rx!=nshell && Ry!=nshell && Rz!=nshell
                            && -Rx!=nshell && -Ry!=nshell && -Rz!=nshell) { continue; }
                        Eigen::Vector3i tmpR = {Rx, Ry, Rz};
                        Eigen::Vector3d Rvector = lattice_vectors.transpose() * (tmpR.cast<double>());
                        double R = Rvector.norm();
                        temp_ewald += std::erfc(R*param_ewald)/R; 
                    } // Rz
                } // Ry
            } // Rx
            if (nshell%2==0)
            {
                term4_Rsum[iatomic_species] += temp_ewald;
                if (std::abs(temp_ewald) < threshold_ewald) { break; }
            }
        } // nshell
        term4_Rsum[iatomic_species] /= 2;
    } // iatomic_species

    double Z_total = 0.0; // should = num. of electrons (in BlochStates)
    for (int iatom=0; iatom<num_atoms; iatom++)
    {
        int iatomic_species = crystal_structure.index_of_atoms()[iatom];
        double Zi = Z_valence[iatomic_species];
        Eigen::Vector3d ri = crystal_structure.atomic_position_cartesian()[iatom];

        for (int jatom=0; jatom<iatom; jatom++) // jatom < iatom, thus no (1/2) factor is required
        {
            int jatomic_species = crystal_structure.index_of_atoms()[jatom];
            double Zj = Z_valence[jatomic_species];
            Eigen::Vector3d rirj = ri - crystal_structure.atomic_position_cartesian()[jatom];

            // Term 1 (G-sum term)
            temp_ewald_sum = 0.0;
            for (int nshell=1; ; nshell++) // Note: G=0 is excluded
            {
                if (nshell%2==1) { temp_ewald = 0.0; } // odd layer
                for (int Gx=-nshell; Gx<=nshell; Gx++)
                {
                    for (int Gy=-nshell; Gy<=nshell; Gy++)
                    {
                        for (int Gz=-nshell; Gz<=nshell; Gz++)
                        {
                            if (Gx!=nshell && Gy!=nshell && Gz!=nshell
                                && -Gx!=nshell && -Gy!=nshell && -Gz!=nshell) { continue; }
                            Eigen::Vector3i tmpG = {Gx, Gy, Gz};
                            Eigen::Vector3d Gvector = reciprocal_vectors.transpose() * (tmpG.cast<double>());
                            double G2 = Gvector.squaredNorm();
                            double G_dot_rirj = Gvector.dot(rirj);
                            temp_ewald += std::exp(-G2*factor_param_ewald)/G2*std::cos(G_dot_rirj);
                        } // Gz
                    } // Gy
                } // Gx
                if (nshell%2==0)
                {
                    temp_ewald_sum += temp_ewald;
                    if (std::abs(temp_ewald) < threshold_ewald) { break; }
                }
            } // nshell
            return_val += temp_ewald_sum * FourPI * Zi * Zj / unit_cell_volume;

            // Term 2 (R-sum term)
            temp_ewald_sum = 0.0;
            for (int nshell=0; ; nshell++) // Note: R=0 is included
            {
                if (nshell%2==0) { temp_ewald = 0.0; } // even layer
                for (int Rx=-nshell; Rx<=nshell; Rx++)
                {
                    for (int Ry=-nshell; Ry<=nshell; Ry++)
                    {
                        for (int Rz=-nshell; Rz<=nshell; Rz++)
                        {
                            if (Rx!=nshell && Ry!=nshell && Rz!=nshell
                                && -Rx!=nshell && -Ry!=nshell && -Rz!=nshell) { continue; }
                            Eigen::Vector3i tmpR = {Rx, Ry, Rz};
                            Eigen::Vector3d Rrirj_vector =
                                lattice_vectors.transpose() * (tmpR.cast<double>()) + rirj;
                            double Rrirj = Rrirj_vector.norm();
                            temp_ewald += std::erfc(Rrirj*param_ewald)/Rrirj; 
                        } // Rz
                    } // Ry
                } // Rx
                if (nshell%2==1)
                {
                    temp_ewald_sum += temp_ewald;
                    if (std::abs(temp_ewald) < threshold_ewald) { break; }
                }
            } // nshell
            return_val += temp_ewald_sum * Zi *Zj;
        } // jatom

        // Terms 3 & 4 & 5
        return_val += Zi*Zi*
            (term3_Gsum[iatomic_species] + term4_Rsum[iatomic_species]
             - param_ewald/std::sqrt(PI));

        // calculate Z_total
        Z_total += Zi;
    } // iatom

    // Term 6: -pi*N^2/(2*Omega*a^2)
    return_val -= PI*Z_total*Z_total/(2*unit_cell_volume*param_ewald*param_ewald);

    return return_val;
}


// public

void TotalEnergy::set_ewald_energy(const bool is_heg, const double &ewald_energy)
{
    if (!is_heg) { ewald_energy_ = ewald_energy; } // ewald_energy = 0 for is_heg==true
}

// used in structural optimization. should not be called for is_heg==true
void TotalEnergy::reset_ewald_energy(const CrystalStructure &crystal_structure,
                                     const bool is_heg, const bool am_i_mpi_rank0)
{
    assert(!is_heg);
    if (am_i_mpi_rank0) { ewald_energy_ = calculate_ewald_energy(crystal_structure); }
    MPI_Bcast(&ewald_energy_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
}

std::vector<Eigen::Vector3d> TotalEnergy::differentiate_ewald_energy(const CrystalStructure &crystal_structure) const
{
    const int num_atoms = crystal_structure.num_atoms();
    const std::vector<double> &Z_valence = crystal_structure.Z_valence();
    const double unit_cell_volume = crystal_structure.unit_cell_volume();
    const Eigen::Matrix3d lattice_vectors = crystal_structure.lattice_vectors();
    const Eigen::Matrix3d reciprocal_vectors = crystal_structure.reciprocal_vectors();

    const double param_ewald = PI/std::cbrt(unit_cell_volume); // cbrt = cubic root
    const double param_ewald2 = param_ewald*param_ewald;
    const double factor_param_ewald = 1.0/(4*param_ewald*param_ewald);
    const double prefactor_ewald = 2*param_ewald/std::sqrt(PI);
    const double threshold_ewald = 1e-8;
    Eigen::Vector3d temp_ewald, temp_ewald_sum; // temporary variable for Ewald sum

    std::vector<Eigen::Vector3d> return_vector(num_atoms);

    for (int iatom=0; iatom<num_atoms; iatom++)
    {
        return_vector[iatom] = {0.0, 0.0, 0.0}; // initialize

        int iatomic_species = crystal_structure.index_of_atoms()[iatom];
        double Zi = Z_valence[iatomic_species];
        Eigen::Vector3d ri = crystal_structure.atomic_position_cartesian()[iatom];

        for (int jatom=0; jatom<num_atoms; jatom++)
        {
            if (iatom==jatom) { continue; }
            int jatomic_species = crystal_structure.index_of_atoms()[jatom];
            double Zj = Z_valence[jatomic_species];
            Eigen::Vector3d rirj = ri - crystal_structure.atomic_position_cartesian()[jatom];

            // Term 1 (G-sum term)
            temp_ewald_sum = {0.0, 0.0, 0.0};
            for (int nshell=1; ; nshell++) // Note: G=0 is excluded
            {
                if (nshell%2==1) { temp_ewald = {0.0, 0.0, 0.0}; } // odd layer
                for (int Gx=-nshell; Gx<=nshell; Gx++)
                {
                    for (int Gy=-nshell; Gy<=nshell; Gy++)
                    {
                        for (int Gz=-nshell; Gz<=nshell; Gz++)
                        {
                            if (Gx!=nshell && Gy!=nshell && Gz!=nshell
                                && -Gx!=nshell && -Gy!=nshell && -Gz!=nshell) { continue; }
                            Eigen::Vector3i tmpG = {Gx, Gy, Gz};
                            Eigen::Vector3d Gvector = reciprocal_vectors.transpose() * (tmpG.cast<double>());
                            double G2 = Gvector.squaredNorm();
                            double G_dot_rirj = Gvector.dot(rirj);

                            temp_ewald += std::exp(-G2*factor_param_ewald)/G2*std::sin(G_dot_rirj)*Gvector;
                        } // Gz
                    } // Gy
                } // Gx
                if (nshell%2==0)
                {
                    temp_ewald_sum += temp_ewald;
                    if (temp_ewald.squaredNorm() < threshold_ewald) { break; }
                }
            } // nshell
            return_vector[iatom] -= temp_ewald_sum * FourPI * Zi * Zj / unit_cell_volume;

            // Term 2 (R-sum term)
            temp_ewald_sum = {0.0, 0.0, 0.0};
            for (int nshell=0; ; nshell++) // Note: R=0 is included
            {
                if (nshell%2==0) { temp_ewald = {0.0, 0.0, 0.0}; } // even layer
                for (int Rx=-nshell; Rx<=nshell; Rx++)
                {
                    for (int Ry=-nshell; Ry<=nshell; Ry++)
                    {
                        for (int Rz=-nshell; Rz<=nshell; Rz++)
                        {
                            if (Rx!=nshell && Ry!=nshell && Rz!=nshell
                                && -Rx!=nshell && -Ry!=nshell && -Rz!=nshell) { continue; }
                            Eigen::Vector3i tmpR = {Rx, Ry, Rz};
                            Eigen::Vector3d Rrirj_vector =
                                lattice_vectors.transpose() * (tmpR.cast<double>()) + rirj;
                            double Rrirj2 = Rrirj_vector.squaredNorm();
                            double Rrirj = std::sqrt(Rrirj2);

                            temp_ewald -= Rrirj_vector/(Rrirj*Rrirj)*
                                (std::erfc(Rrirj*param_ewald)/Rrirj
                                 + prefactor_ewald*std::exp(-param_ewald2*Rrirj2)); 
                        } // Rz
                    } // Ry
                } // Rx
                if (nshell%2==1)
                {
                    temp_ewald_sum += temp_ewald;
                    if (temp_ewald.squaredNorm() < threshold_ewald) { break; }
                }
            } // nshell
            return_vector[iatom] += temp_ewald_sum * Zi * Zj;
        } // jatom
    } // iatom

    return return_vector;
}
