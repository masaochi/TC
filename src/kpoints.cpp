// [class Kpoints]
// list of k-point, indicies of symmetry operations for k-points

#include "include/header.hpp"

namespace
{
    
// check whether vec is included in vec_list
bool is_matched(const Eigen::Vector3d &vec, const std::vector<Eigen::Vector3d> &vec_list)
{
    double diff, norm;
    for (auto v : vec_list)
    {
        norm = 0.0;
        for (int idim=0; idim<3; idim++)
        {
            diff = vec(idim) - v(idim);
            while (diff < -0.5) { diff += 1.0; }
            while (diff > 0.5) { diff -= 1.0; }
            norm += diff*diff;
        }
        if (norm < 1e-8) { return true; }
    }
    return false;
}    

} // namespace

// private

// called only in the "SCF" mode. called only by the rank=0 MPI process.
// set kvectors_scf_[ik][is] using kvectors_scf_[ik][0] + symmetry
void Kpoints::set_kvectors_by_symmetry(const Spin &spin, const Symmetry &symmetry, std::ostream *ost)
{
    Eigen::Vector3d vec_tmp(3);
    std::vector<Eigen::Vector3d> vec_list(0);
    for (int ik=0; ik<num_irreducible_kpoints_scf_; ik++)
    {
        if (kweight_scf_[ik] < 1e-8) { continue; } // skip for zero k-weight
        vec_list.push_back(kvectors_scf_[ik][0]);
    }
    
    for (int ik=0; ik<num_irreducible_kpoints_scf_; ik++) 
    {
        if (kweight_scf_[ik] < 1e-8) { continue; } // skip for zero k-weight
        for(int isym=0; isym<symmetry.num_symmetry_operations(); isym++) 
        {
            // k' = U^+ k (Note: rotation matrix is transposed)
            vec_tmp = symmetry.rotation()[isym].cast<double>().transpose() * kvectors_scf_[ik][0]; // 0 = irreducible k-point
            for (int idim=0; idim<3; idim++) { if (std::abs(vec_tmp(idim))<1e-4) { vec_tmp(idim) = 0.0; } }  // almost 0
            
            if (!is_matched(vec_tmp, vec_list)) {
                vec_list.push_back(vec_tmp);
                kvectors_scf_[ik].push_back(vec_tmp);

                // symmetry operation used here
                index_of_rotation_at_k_[ik].push_back(isym);
                is_time_reversal_used_at_k_[ik].push_back(false);
            }
            
            // k' -> -k'
            if (symmetry.is_time_reversal_symmetric() && !is_matched(-vec_tmp, vec_list)) {
                vec_list.push_back(-vec_tmp);
                kvectors_scf_[ik].push_back(-vec_tmp);

                // symmetry operation used here
                index_of_rotation_at_k_[ik].push_back(isym);
                is_time_reversal_used_at_k_[ik].push_back(true);
            }
        } // isym
    } // ik

    num_kpoints_ = 0; // total num. of non-zero-weight k-points
    num_kpoints_all_scf_ = 0; // total num. of k-points including zero-weight ones
    for (int ik=0; ik<num_irreducible_kpoints_scf_; ik++) 
    {
        if (kweight_scf_[ik] > 1e-8) { num_kpoints_ += kvectors_scf_[ik].size(); } // non-zero k-weight
        num_kpoints_all_scf_ += kvectors_scf_[ik].size();
    }

    index_all_kscf_.resize(num_kpoints_all_scf_);
    int ik_isymk = 0;
    for (int ik=0; ik<num_irreducible_kpoints_scf_; ik++) 
    {
        for (int isymk=0; isymk<kvectors_scf_[ik].size(); isymk++)
        {
            index_all_kscf_[ik_isymk] = {ik, isymk}; // i.e. index_all_kscf_[ik_isymk][0 or 1] = ik or isymk
            ik_isymk++;
        }
    }
    assert(ik_isymk == num_kpoints_all_scf_);
    
    // check consistency with k_weight
    double fact = (spin.num_independent_spins()==1 && !spin.is_spinor()) ? 2.0 : 1.0; // 2 for no-spin
    fact /= static_cast<double>(num_kpoints_);
    for (int ik=0; ik<num_irreducible_kpoints_scf_; ik++) 
    {
        if ((kweight_scf_[ik] > 1e-6) &&
            (std::abs(fact*kvectors_scf_[ik].size() - kweight_scf_[ik]) > 1e-6))
        {
            std::cout << "ik: " << ik << "  kweight in QE: " << kweight_scf_[ik] << ", kweight in TC: " << fact*kvectors_scf_[ik].size() << std::endl;
            error_messages::stop("Error in k-weight. (e.g., a k-mesh that breaks a space-group symmetry is not applicable)"); 
        }
    }

    // output
    show_kvectors(symmetry, ost);
}

// called in set_kvectors_by_symmetry() (SCF calc. only)
void Kpoints::show_kvectors(const Symmetry &symmetry, std::ostream *ost) const
{
    *ost << "  Print symmetrically-repiclated kvectors:" << std::endl;

   for (int ik=0; ik<num_irreducible_kpoints_scf_; ik++) 
    {
        *ost << "  Irreducible k-vector (ik = " << ik << ") "
             << kvectors_scf_[ik][0](0) << " " 
             << kvectors_scf_[ik][0](1) << " " 
             << kvectors_scf_[ik][0](2) << std::endl;
        for (int isymk=0; isymk<kvectors_scf_[ik].size(); isymk++) 
        {
            *ost << "   symmetrically-equivalent k-vector (isymk = " << isymk << "): "
                 << kvectors_scf_[ik][isymk](0) << " " 
                 << kvectors_scf_[ik][isymk](1) << " " 
                 << kvectors_scf_[ik][isymk](2) << std::endl;
            *ost << "   transformation operator (in k-space, applied to the irreducible k. Transpose of the symmorphic operator in r-space) " << std::endl; 
            *ost << symmetry.rotation()[index_of_rotation_at_k_[ik][isymk]].transpose() << std::endl;
            *ost << "   time-reversal symmetry = ";
            if (is_time_reversal_used_at_k_[ik][isymk]) { *ost << "true" << std::endl;}
            else { *ost << "false" << std::endl; }
            *ost << "   translation vector (while not applied to the k-vector) = ";
            *ost << symmetry.translation()[index_of_rotation_at_k_[ik][isymk]].transpose() << std::endl;
            *ost << std::endl;
        }
    }
   *ost << std::endl;
}

// public

void Kpoints::set_tc_input(const std::string&smearing_mode,const double &smearing_width)
{
    const std::vector<std::string> list_smearing_mode{"fixed", "gaussian"};
    error_messages::compare_list("smearing_mode", smearing_mode, list_smearing_mode); // stop the program if not included in list_smearing_mode
    smearing_mode_ = smearing_mode;

    if (smearing_width_ > -1e-8) { smearing_width_ = smearing_width; } // negative value will be discarded
}

void Kpoints::bcast_tc_input(const bool am_i_mpi_rank0)
{
    int length = smearing_mode_.size();
    MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::string s(length, 'a'); // aaaaa when length=5
    if (am_i_mpi_rank0) { s = smearing_mode_; }
    MPI_Bcast(const_cast<char*>(s.data()), length, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (!am_i_mpi_rank0) { smearing_mode_ = s; }

    MPI_Bcast(&smearing_width_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void Kpoints::set_qe_input(const Spin &spin, const Symmetry &symmetry,
                           const int num_irreducible_kpoints,
                           const std::vector<double> kweight,
                           const std::vector<Eigen::Vector3d> &kvectors_irred,
                           const std::string &calc_mode,
                           std::ostream *ost)
{
    assert(num_irreducible_kpoints > 0);
    assert(kweight.size() == num_irreducible_kpoints);
    assert(kvectors_irred.size() == num_irreducible_kpoints);
    assert(calc_mode == "SCF" || calc_mode == "BAND");

    if (calc_mode=="SCF")
    {
        num_irreducible_kpoints_scf_ = num_irreducible_kpoints;
        kweight_scf_ = kweight;

        kvectors_scf_.resize(num_irreducible_kpoints);
        // set irreducible k-points
        for (int ik=0; ik<kvectors_irred.size(); ik++) 
        {
            kvectors_scf_[ik].resize(1);
            kvectors_scf_[ik][0] = kvectors_irred[ik];
            for (int idim=0; idim<3; idim++) 
            {
                if (std::abs(kvectors_scf_[ik][0](idim))<1e-6) { kvectors_scf_[ik][0](idim) = 0.0; }  // almost 0
            }
        }

        // In the following, we assume that rotation[0] = identity matrix
        const Eigen::Matrix3i id = Eigen::Matrix3i::Identity();
        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
            {
                assert(symmetry.rotation()[0](i,j)==id(i,j));
            }
        }
        
        // set default values
        index_of_rotation_at_k_.resize(num_irreducible_kpoints);
        is_time_reversal_used_at_k_.resize(num_irreducible_kpoints);
        for (int ik=0; ik<kvectors_irred.size(); ik++) 
        {
            index_of_rotation_at_k_[ik].resize(1);
            is_time_reversal_used_at_k_[ik].resize(1);
            
            index_of_rotation_at_k_[ik][0] = 0; // identity matrix
            is_time_reversal_used_at_k_[ik][0] = false;
        }

        // make symmetrically-equivalent k-points
        // num_kpoints_ & num_kpoints_all_scf_ & index_all_kscf_ are also set here.
        set_kvectors_by_symmetry(spin, symmetry, ost);
    }
    else if (calc_mode=="BAND")
    {
        num_irreducible_kpoints_band_ = num_irreducible_kpoints;

        kvectors_band_.resize(num_irreducible_kpoints);
        for (int ik=0; ik<kvectors_irred.size(); ik++) 
        {
            kvectors_band_[ik].resize(1);
            kvectors_band_[ik][0] = kvectors_irred[ik];
            for (int idim=0; idim<3; idim++) 
            {
                if (std::abs(kvectors_band_[ik][0](idim))<1e-6) { kvectors_band_[ik][0](idim) = 0.0; }  // almost 0
            }
        }
    }
}

void Kpoints::bcast_qe_input(const std::string &calc_mode,
                             const bool am_i_mpi_rank0)
{
    assert(calc_mode == "SCF" || calc_mode == "BAND");
    if (calc_mode=="SCF")
    {
        MPI_Bcast(&num_kpoints_, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&num_kpoints_all_scf_, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if (!am_i_mpi_rank0) { index_all_kscf_.resize(num_kpoints_all_scf_); }
        for (int ik_isymk=0; ik_isymk<num_kpoints_all_scf_; ik_isymk++) 
        {
            if (!am_i_mpi_rank0) { index_all_kscf_[ik_isymk].resize(2); }
            MPI_Bcast(&index_all_kscf_[ik_isymk][0], 2, MPI_INT, 0, MPI_COMM_WORLD);
        }

        MPI_Bcast(&num_irreducible_kpoints_scf_, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!am_i_mpi_rank0) { kvectors_scf_.resize(num_irreducible_kpoints_scf_); }
        for (int ik=0; ik<kvectors_scf_.size(); ik++)
        {
            int i = kvectors_scf_[ik].size();
            MPI_Bcast(&i, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (!am_i_mpi_rank0) { kvectors_scf_[ik].resize(i); }
            
            for (int isym=0; isym<i; isym++)
            {
                MPI_Bcast(kvectors_scf_[ik][isym].data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
        }

        if (!am_i_mpi_rank0) { kweight_scf_.resize(num_irreducible_kpoints_scf_); }
        MPI_Bcast(&kweight_scf_[0], kweight_scf_.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // symmetry operation at k
        if (!am_i_mpi_rank0)
        {
            index_of_rotation_at_k_.resize(num_irreducible_kpoints_scf_);
            is_time_reversal_used_at_k_.resize(num_irreducible_kpoints_scf_);
            for (int ik=0; ik<kvectors_scf_.size(); ik++)
            {
                index_of_rotation_at_k_[ik].resize(kvectors_scf_[ik].size());
                is_time_reversal_used_at_k_[ik].resize(kvectors_scf_[ik].size());
            }
        }
        std::vector<int> ivec;
        for (int ik=0; ik<kvectors_scf_.size(); ik++)
        {
            int nsym = kvectors_scf_[ik].size();
            MPI_Bcast(&index_of_rotation_at_k_[ik][0], nsym, MPI_INT, 0, MPI_COMM_WORLD);
            
            ivec.resize(nsym);
            for (int is=0; is<nsym; is++) {
                ivec[is] = is_time_reversal_used_at_k_[ik][is] ? 1 : 0;
            }
            MPI_Bcast(&ivec[0], nsym, MPI_INT, 0, MPI_COMM_WORLD);
            if (!am_i_mpi_rank0)
            {
                for (int is=0; is<nsym; is++) {
                    is_time_reversal_used_at_k_[ik][is] = ivec[is]==1 ? true : false;
                }
            }
        }
    }
    else if (calc_mode=="BAND")
    {
        MPI_Bcast(&num_irreducible_kpoints_band_, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (!am_i_mpi_rank0) { kvectors_band_.resize(num_irreducible_kpoints_band_); }
        for (int ik=0; ik<kvectors_band_.size(); ik++)
        {
            if (!am_i_mpi_rank0) { kvectors_band_[ik].resize(1); }            
            MPI_Bcast(kvectors_band_[ik][0].data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
}

void Kpoints::set_scfinfo(const std::vector<double> &kweight_scf,
                          const std::vector<std::vector<Eigen::Vector3d> > &kvectors_scf,
                          const std::vector<std::vector<int> > &index_of_rotation_at_k,
                          const std::vector<std::vector<bool> > &is_time_reversal_used_at_k)
{
    const int num_irreducible_kpoints = kweight_scf.size();
    num_irreducible_kpoints_scf_ = num_irreducible_kpoints;

    kweight_scf_ = kweight_scf;
    kvectors_scf_ = kvectors_scf;
    index_of_rotation_at_k_ = index_of_rotation_at_k;
    is_time_reversal_used_at_k_ = is_time_reversal_used_at_k;

    // also set num_kpoints & num_kpoints_all_scf_
    num_kpoints_ = 0; // total num. of non-zero-weight k-points
    num_kpoints_all_scf_ = 0; // total num. of scf k-points including zero-weight ones
    for (int ik=0; ik<num_irreducible_kpoints_scf_; ik++) 
    {
        if (kweight_scf_[ik] > 1e-8) { num_kpoints_ += kvectors_scf_[ik].size(); }
        num_kpoints_all_scf_ += kvectors_scf_[ik].size();
    }

    index_all_kscf_.resize(num_kpoints_all_scf_);
    int ik_isymk = 0;
    for (int ik=0; ik<num_irreducible_kpoints_scf_; ik++) 
    {
        for (int isymk=0; isymk<kvectors_scf_[ik].size(); isymk++)
        {
            index_all_kscf_[ik_isymk] = {ik, isymk}; // i.e. index_all_kscf_[ik_isymk][0 or 1] = ik or isymk
            ik_isymk++;
        }
    }
    assert(ik_isymk == num_kpoints_all_scf_);
}

// return band filling including the spin factor
double Kpoints::return_band_filling(const Spin &spin, const double &orbital_energy, const double &fermi_energy,
                                    const int &num_electrons, const int &iband) const
{
    const double spin_factor = (spin.num_independent_spins()==1 && !spin.is_spinor()) ? 2.0 : 1.0; // 2 for no-spin
    if (smearing_mode_=="fixed") // fermi_energy is not used
    {
        int num_occupied_bands;
        if (spin.is_spinor())
        {
            num_occupied_bands = std::round(num_electrons);
            if (std::abs(num_occupied_bands - num_electrons)>1e-5)
            {
                error_messages::stop("fixed smearing cannot be used for this num_electrons in return_band_filling()");
            }
        }
        else
        {
            num_occupied_bands = std::round(num_electrons/2.0);
            if (std::abs(num_occupied_bands - num_electrons/2.0)>1e-5)
            {
                error_messages::stop("fixed smearing cannot be used for this num_electrons in return_band_filling()");
            }
        }
        
        if (iband < num_occupied_bands)
        {
            return spin_factor;
        }
        else
        {
            return 0.0;
        }
    }
    else if (smearing_mode_=="gaussian") // num_electrons and iband are not used
    {
        double x = (orbital_energy - fermi_energy)/smearing_width_;
        return spin_factor*0.5*(1.0 - std::erf(x));
    }
}
