// [class TotalEnergy]
// total energy, ewald energy, etc.

#include "include/header.hpp"

void TotalEnergy::set_force(const std::vector<Eigen::Vector3cd> &force)
{
    force_ = force;
}

void TotalEnergy::set_force_old()
{
    force_old_ = force_;
}

void TotalEnergy::set_Hessian_inv_old(const Eigen::MatrixXd &Hessian_inv_old)
{
    Hessian_inv_old_ = Hessian_inv_old;
}

void TotalEnergy::print_force(const bool is_scf_converged, std::ostream *ost) const
{
    const int num_atoms = force_.size();
    const double factor = Ht_in_eV / Bohr_in_ang; // a.u. (i.e., Hartree/Bohr) to eV/ang

    *ost << std::endl;
    if (!is_scf_converged)
    {
        *ost << "   Hellmann-Feynman force (eV/ang) in Cartesian coordinate" << std::endl;
        *ost << "   NOTE: SCF eq. has not been converged. Force might be incorrect." << std::endl;
    }
    else
    {
        *ost << "   Hellmann-Feynman force (eV/ang) in Cartesian coordinate (with converged SCF)" << std::endl;
    }
    *ost << "    real part" << std::endl;
    for (int iatom=0; iatom<num_atoms; iatom++)
    {
        *ost << "    Atom " << iatom 
             << " ( " << factor * force_[iatom](0).real()
             << " , " << factor * force_[iatom](1).real()
             << " , " << factor * force_[iatom](2).real() << " )" << std::endl;;
    } // iatom

    *ost << "    imagniary part" << std::endl;
    for (int iatom=0; iatom<num_atoms; iatom++)
    {
        *ost << "    Atom " << iatom 
             << " ( " << factor * force_[iatom](0).imag()
             << " , " << factor * force_[iatom](1).imag()
             << " , " << factor * force_[iatom](2).imag() << " )" << std::endl;;
    } // iatom
    *ost << std::endl;
}

double TotalEnergy::maximum_force_in_eV_ang() const
{
    const double factor = Ht_in_eV / Bohr_in_ang; // a.u. (i.e., Hartree/Bohr) to eV/ang

    double maximum_force = 0.0;
    for (int iatom=0; iatom<force_.size(); iatom++)
    {
        double temp_force = 0.0;
        for (int idim=0; idim<3; idim++)
        {
            temp_force += force_[iatom](idim).real() * force_[iatom](idim).real();
        }
        temp_force = std::sqrt(temp_force);
        if (maximum_force < temp_force) { maximum_force = temp_force; }
    }
    return factor * maximum_force;
}

