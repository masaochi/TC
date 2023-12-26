// [class Potentials]
// pseudopotential, Jastrow, Ewald sum

// given in atomic unit (Hartree)

#include "include/header.hpp"

// private

// Simpson's rule integration on the pseudo-potential radial mesh. Integrate "function".
double Potentials::simpson(const std::vector<double> &function, const int &iatomic_species) const
{
    int mesh_size = pseudo_weight_rmesh_[iatomic_species].size();
    assert(function.size() == mesh_size);
    assert(iatomic_species < pseudo_weight_rmesh_.size());

    double sum = 0.0;
    for (int imesh=0; imesh<mesh_size-2; imesh=imesh+2)
    {
        sum += function[imesh]*pseudo_weight_rmesh_[iatomic_species][imesh];
        sum += function[imesh+1]*pseudo_weight_rmesh_[iatomic_species][imesh+1]*4;
        sum += function[imesh+2]*pseudo_weight_rmesh_[iatomic_species][imesh+2];
    }
    return sum/3.0;
}

// Rydberg -> Hartree
// (also required for pseudo_nlcc when TC++ supports non-linear core correction...)
void Potentials::set_unit(std::vector<std::vector<double> > &pseudo_local,
                          std::vector<Eigen::MatrixXd> &pseudo_dij) const
{
    const int num_atomic_species = pseudo_local.size();
    assert(pseudo_local.size() == pseudo_dij.size());
    for (int iatomic_species=0; iatomic_species<num_atomic_species; iatomic_species++)
    {
        int mesh_size = pseudo_local[iatomic_species].size();
        int num_projectors = pseudo_dij[iatomic_species].rows();

        assert(pseudo_dij[iatomic_species].rows() == pseudo_dij[iatomic_species].cols());
        for (int imesh=0; imesh<mesh_size; imesh++) { pseudo_local[iatomic_species][imesh] /= 2.0; }
        pseudo_dij[iatomic_species] /= 2.0;
    }    
}

// calculate an "iatomic_species" component
void Potentials::set_Vpp_local_atom(const CrystalStructure &crystal_structure,
                                    PlaneWaveBasis &plane_wave_basis,
                                    Eigen::VectorXcd &Vpp_local_atom_G,
                                    const int &iatomic_species) const
{
    const int npw = plane_wave_basis.size_FFT_grid();
    const int num_atoms = crystal_structure.num_atoms();
    const int num_atomic_species = crystal_structure.num_atomic_species();
    const std::vector<double> &Z_valence = crystal_structure.Z_valence();
    const double unit_cell_volume = crystal_structure.unit_cell_volume();
    const double param_ewald = PI/std::cbrt(unit_cell_volume); // cbrt = cubic root

    assert(iatomic_species >= 0 && iatomic_species < num_atomic_species);
    std::vector<double> local_temp, local_temp2;

    Eigen::MatrixXd Gvector(npw, 3);
    Eigen::VectorXd Gnorm(npw);
    for (int ipw=0; ipw<npw; ipw++)
    {
        // Gvector[ipw] = i1*b1 + i2*b2 + i3*b3. (ipw = (i1, i2, i3))
        Gvector.row(ipw) =
            ((plane_wave_basis.get_Gvector(ipw)).cast<double>()).transpose()
            * crystal_structure.reciprocal_vectors();
        Gnorm(ipw) = std::sqrt(Gvector.row(ipw) * Gvector.row(ipw).transpose()); // |G|
    }

    // set Vpp_local_atom_G (Vpp_local_G for each atom) (G = in G space)
    Vpp_local_atom_G.resize(npw);

    // Ewald sum (usual): Ewald parameter = a, 
    //   \sum_R 1/|r-R| = \sum_R erfc(a|r-R|)/|r-R| + 4pi/V \sum_G exp(-G^2/(4a^2))/(G^2) exp(iGr)
    // Formula used here, which can be obtained from the above equation: (Note! erfc = 1 - erf)
    //   \sum_R erf(a|r-R|)/|r-R| = 4pi/V \sum_G exp(-G^2/(4a^2))/(G^2) exp(iGr)
    // Below we call left-hand-side "term 1" and right-hand-side "term 2", respectively.

    // For evaluating Vpp_local_atom_G,
    // (1) subtract "term 1" from pseudo_local to make a short-ranged function
    //     i.e., make Vpp_short[r] = Vpp[r] - Vpp_long[r]
    // (2) perform Fourier-Transform of Vpp_short (here, summation over lattice-periodicity is performed)
    // (3) add "term 2" to compensate (1). 
    //     At the same time, G=0 of -Z/r is subtracted, which should be canceled with G=0 of el-el + nucl-nucl.
    int mesh_size = pseudo_rmesh_[iatomic_species].size();
    local_temp.resize(mesh_size);
    local_temp2.resize(mesh_size);
    // (1) subtract "term 1" in R-space = subtract (-Z/r)*.... = add (Z/r)*...
    for (int imesh=0; imesh<mesh_size; imesh++) // radial mesh in pseudopot.
    {
        if (pseudo_rmesh_[iatomic_species][imesh] > 1e-8)
        {
            local_temp[imesh]
                = pseudo_local_[iatomic_species][imesh]
                + (Z_valence[iatomic_species]/pseudo_rmesh_[iatomic_species][imesh])
                * std::erf(param_ewald*pseudo_rmesh_[iatomic_species][imesh]);
        }
        else // erf(ax)/x = a*(erf(ax)/(ax)) -> 2a/sqrt(PI) for x->0
        {
            local_temp[imesh]
                = pseudo_local_[iatomic_species][imesh]
                + Z_valence[iatomic_species]*2.0*param_ewald/std::sqrt(PI);
        }
    } // imesh
    
    for (int ipw=0; ipw<npw; ipw++)
    {
        // multiply r^2*sin(|G|r)/(|G|r), see below
        for (int imesh=0; imesh<mesh_size; imesh++) // radial mesh in pseudopot.
        {
            double rmesh = pseudo_rmesh_[iatomic_species][imesh];
            double Gr = Gnorm(ipw) * rmesh;
            if (Gr>0.1)
            {
                local_temp2[imesh] = local_temp[imesh]*std::sin(Gr)/Gr;
            }
            else // For small Gr, use Taylor expansion
            {
                double temp = 1.0 - Gr*Gr/(11*10);
                temp = 1.0 - temp*Gr*Gr/(9*8);
                temp = 1.0 - temp*Gr*Gr/(7*6);
                temp = 1.0 - temp*Gr*Gr/(5*4);
                temp = 1.0 - temp*Gr*Gr/(3*2);
                local_temp2[imesh] = local_temp[imesh]*temp;
            }
            local_temp2[imesh] *= (rmesh*rmesh);
        }
        
        // (2) Fourier-Transform using Simpson's rule integration
        // (1/V) \int pseudo_local(r) exp(-iGr) d^3 r = (4pi/V) \int pseudo_local(r) (sin(|G|r)/|G|r) r^2 dr
        Vpp_local_atom_G(ipw) = simpson(local_temp2, iatomic_species);
        
        // (3) add "term2" in G-space
        if (ipw!=0)
        {
            Vpp_local_atom_G(ipw)
                -= Z_valence[iatomic_species]
                *std::exp(-Gnorm(ipw)*Gnorm(ipw)/(4.0*param_ewald*param_ewald))/(Gnorm(ipw)*Gnorm(ipw));
        }
        else
        {
            // Note! we consider G=0 limit of [exp(-G^2/(4*param^2)) - 1]/(G^2) -> -1/(4*param^2),
            // instead of exp(-G^2/(4*param^2))/(G^2).
            // This is because G=0 of -Z/r should be canceled with G=0 of el-el + nucl-nucl.
            Vpp_local_atom_G(ipw)
                += Z_valence[iatomic_species]/(4.0*param_ewald*param_ewald);
        }
        Vpp_local_atom_G(ipw) *= (FourPI/unit_cell_volume);
    } // ipw
}

void Potentials::set_Vpp_local(const CrystalStructure &crystal_structure,
                               PlaneWaveBasis &plane_wave_basis)
{
    const int npw = plane_wave_basis.size_FFT_grid();
    const int num_atoms = crystal_structure.num_atoms();
    const int num_atomic_species = crystal_structure.num_atomic_species();

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

    // set Vpp_local_G using Vpp_local_atom_G
    // Vpp_local(r) = \sum_{iatom} Vpp_local_atom(r-r_{iatom})
    // where FT[Vpp_local_atom(r-r_{iatom})] (G) = FT[Vpp_local_atom(r)] (G) * exp^{-iGr_{iatom}}
    Eigen::VectorXcd Vpp_local_G = Eigen::VectorXcd::Zero(npw);
    for (int iatom=0; iatom<num_atoms; iatom++)
    {
        int iatomic_species = crystal_structure.index_of_atoms()[iatom];
        Eigen::VectorXd GdotR = Gvector * crystal_structure.atomic_position_cartesian()[iatom]; // array(npw)
        for (int ipw=0; ipw<npw; ipw++)
        {
            Vpp_local_G(ipw) += (std::cos(GdotR(ipw)) - I*std::sin(GdotR(ipw))) 
                * Vpp_local_atom_G[iatomic_species](ipw);
        } // ipw
    } // iatom
    plane_wave_basis.FFT_backward(Vpp_local_G, Vpp_local_G); // G-space to R-space

    Vpp_local_.resize(npw);
    for (int ipw=0; ipw<npw; ipw++) { Vpp_local_[ipw] = Vpp_local_G[ipw].real(); } // take a real part
}

// public

void Potentials::set_upf(const std::vector<std::vector<double> > &pseudo_rmesh,
                         const std::vector<std::vector<double> > &pseudo_weight_rmesh,
                         std::vector<std::vector<double> > &pseudo_local,
                         const std::vector<std::vector<std::vector<double> > > &pseudo_beta,
                         const std::vector<std::vector<int> > &pseudo_lbeta,
                         std::vector<Eigen::MatrixXd> &pseudo_dij,
                         const CrystalStructure &crystal_structure,
                         PlaneWaveBasis &plane_wave_basis)
{   
    set_unit(pseudo_local, pseudo_dij); // Rydberg -> Hartree

    pseudo_rmesh_ = pseudo_rmesh;
    pseudo_weight_rmesh_ = pseudo_weight_rmesh;
    pseudo_local_ = pseudo_local;
    pseudo_beta_ = pseudo_beta;

    pseudo_lbeta_ = pseudo_lbeta;
    pseudo_dij_ = pseudo_dij;

    // set Vpp_local_
    set_Vpp_local(crystal_structure, plane_wave_basis);
}

void Potentials::bcast_upf(const bool am_i_mpi_rank0)
{
    int npw = Vpp_local_.rows();
    MPI_Bcast(&npw, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!am_i_mpi_rank0) { Vpp_local_.resize(npw); }
    MPI_Bcast(Vpp_local_.data(), npw, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int num_atomic_species = pseudo_rmesh_.size();
    MPI_Bcast(&num_atomic_species, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!am_i_mpi_rank0)
    {
        pseudo_rmesh_.resize(num_atomic_species);
        pseudo_weight_rmesh_.resize(num_atomic_species);
        pseudo_beta_.resize(num_atomic_species);

        // pseudo_local_ is not Bcasted (used only in rank0)

        pseudo_lbeta_.resize(num_atomic_species);
        pseudo_dij_.resize(num_atomic_species);
    }
    for (int iatomic_species=0; iatomic_species<num_atomic_species; iatomic_species++)
    {
        int mesh_size = pseudo_rmesh_[iatomic_species].size();
        MPI_Bcast(&mesh_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int num_projectors = pseudo_beta_[iatomic_species].size();
        MPI_Bcast(&num_projectors, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (!am_i_mpi_rank0)
        {
            pseudo_rmesh_[iatomic_species].resize(mesh_size);
            pseudo_weight_rmesh_[iatomic_species].resize(mesh_size);
            pseudo_beta_[iatomic_species].resize(num_projectors);
            for (int iproj=0; iproj<num_projectors; iproj++)
            {
                pseudo_beta_[iatomic_species][iproj].resize(mesh_size);
            }
            pseudo_lbeta_[iatomic_species].resize(num_projectors);
            pseudo_dij_[iatomic_species].resize(num_projectors, num_projectors);
        }
        MPI_Bcast(&pseudo_rmesh_[iatomic_species][0], mesh_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&pseudo_weight_rmesh_[iatomic_species][0], mesh_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (int iproj=0; iproj<num_projectors; iproj++)
        {
            MPI_Bcast(&pseudo_beta_[iatomic_species][iproj][0], mesh_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        MPI_Bcast(&pseudo_lbeta_[iatomic_species][0], num_projectors, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(pseudo_dij_[iatomic_species].data(), num_projectors*num_projectors, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

// return <beta|k+G> except the exp(i(k+G)R_at) factor
// By using the Reyleigh expansion
//   exp(i(k+G)r) = sum_{l,m} 4pi*i^l j_l(|k+G| |r|) Y_lm(e_{K+G})^* Y_lm(e_r)
//   [Note: For k+G=0, j_0(0) = 1, Y_00 = 1/sqrt(4pi) have non-zero contribution. (j_l(0)=0 for l!=0)]
// we evaluate <beta(r) Y_lm(e_r)|k+G> where the atom is place at R_at:
// = \int d^3r beta^*(r-R_at) Y_lm(e_{r-R_at})^* exp(i(k+G)r)/sqrt(V)
// = exp(i(k+G)R_at) \int d^3r beta^*(r) Y_lm(e_r)^* exp(i(k+G)r)/sqrt(V)
// = exp(i(k+G)R_at) i^l Y_lm(e_{k+G})^* [4pi/sqrt(V) \int dr r^2 beta^*(r) j_l(|k+G|r)]
//
// In projector_nonlocal(), exp(i(k+G)R_at) is not multiplied.
// i^l is also eliminated since projector(l) is always multiplied with projector(l)^*.
std::vector<Eigen::VectorXcd> Potentials::projector_nonlocal(const CrystalStructure &crystal_structure,
                                                             const PlaneWaveBasis &plane_wave_basis,
                                                             const Eigen::VectorXi &Gindex_at_k,
                                                             const Eigen::Vector3d &kvector,
                                                             const std::vector<std::vector<Eigen::VectorXcd> > &ylm,
                                                             const int &iatomic_species, const int &iproj) const
{
    const int num_atomic_species = pseudo_rmesh_.size();
    const int num_G_at_k = Gindex_at_k.size();
    const int num_projectors = pseudo_lbeta_[iatomic_species].size();
    const int lbeta = pseudo_lbeta_[iatomic_species][iproj];
    const int mesh_size = pseudo_rmesh_[iatomic_species].size();

    if (lbeta >= ylm.size()) { error_messages::stop("lbeta too large: need to modify the code."); }
    // m = 0, 1, ..., l. For negative m, we just use (Y_lm)^*. (-1)^m does not matter when pseudo_dij is diagonal w.r.t. l.
    assert(ylm[lbeta].size() == lbeta + 1); 

    // m = -l, -l+1, ..., l (2*lbeta +1 indices)
    std::vector<Eigen::VectorXcd> return_vector(2*lbeta + 1); // return this array
    for (int im=0; im<2*lbeta+1; im++) { return_vector[im].resize(num_G_at_k); }

    double beta_jl;
    std::vector<double> temp_rfunc(mesh_size);

    const double factor = FourPI/std::sqrt(crystal_structure.unit_cell_volume());
// i^l is not necessary since projector(l) is always multiplied with projector(l)^*.
//    Complex cfactor = 1.0;
//    for (int il=0; il<lbeta; il++) { cfactor *= I; } // I^l

    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++) 
    {
        Eigen::Vector3d kGvect = crystal_structure.reciprocal_vectors().transpose()
            * (kvector + plane_wave_basis.get_Gvector(Gindex_at_k(ipw_at_k)).cast<double>());
        double kG = kGvect.norm();

        // calculate beta_jl(|k+G|) = [4pi/sqrt(V) \int dr r^2 beta^*(r) j_l(|k+G|r)].
        for (int imesh=0; imesh<mesh_size; imesh++) // radial mesh in pseudopot.
        {
            double rmesh = pseudo_rmesh_[iatomic_species][imesh];
            temp_rfunc[imesh] = rmesh * pseudo_beta_[iatomic_species][iproj][imesh] // beta in UPF = r*beta (Length^{-1/2})
                * boost::math::sph_bessel(lbeta, kG*rmesh);
//            temp_rfunc[imesh] = rmesh*rmesh*pseudo_beta_[iatomic_species][iproj][imesh]
//                * boost::math::sph_bessel(lbeta, kG*rmesh);
        } // imesh
        beta_jl = factor*simpson(temp_rfunc, iatomic_species);

        // calculate Y_lm(e_{k+G})^* beta_jl(|k+G|)
        // Here, cfactor (i^l) is eliminated.
        return_vector[lbeta](ipw_at_k) = ylm[lbeta][0](ipw_at_k) * beta_jl; // im=0 
        for (int im=1; im<lbeta+1; im++)
        {
            return_vector[lbeta - im](ipw_at_k) = std::conj(ylm[lbeta][im](ipw_at_k)) * beta_jl;
            return_vector[lbeta + im](ipw_at_k) = ylm[lbeta][im](ipw_at_k) * beta_jl;
        }
    }
    return return_vector;
}

void Potentials::reset_Vpp_local(const CrystalStructure &crystal_structure,
                                 PlaneWaveBasis &plane_wave_basis,
                                 const bool am_i_mpi_rank0)
{
    const int npw = plane_wave_basis.size_FFT_grid();

    if (am_i_mpi_rank0)
    {
        set_Vpp_local(crystal_structure, plane_wave_basis);
    }
    else
    {
        Vpp_local_.resize(npw); 
    }
    MPI_Bcast(Vpp_local_.data(), npw, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
