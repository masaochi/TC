// [class Diagonalization]
// scf_loop, convergence check, control paramters, etc.

#include "include/header.hpp"

// make target_vector orthonormalized to vectors[i] (0 <= i < target_band)
bool Diagonalization::Gram_Schmidt(const std::vector<Eigen::VectorXcd> &vectors, const int &target_band,
                                   Eigen::VectorXcd &target_vector) const
{
    assert(target_band <= vectors.size());
    Complex coeff;

    for (int iband=0; iband<target_band; iband++) 
    {
        coeff = vectors[iband].dot(target_vector); // <iband | target_vector>
        target_vector -= coeff * vectors[iband]; // |target_vector> -= <iband | target_vector> * |iband>
    }
    double norm = target_vector.norm();

    if (norm < 1e-8) 
    {
        return false;
    }
    else
    {
        target_vector /= norm;
        return true;
    }
}

bool Diagonalization::Gram_Schmidt_biortho(const std::vector<Eigen::VectorXcd> &vectors_right, 
                                           const std::vector<Eigen::VectorXcd> &vectors_left,
                                           const int &target_band,
                                           Eigen::VectorXcd &target_vector_right,
                                           Eigen::VectorXcd &target_vector_left) const
{
    assert(target_band < vectors_right.size() && target_band < vectors_left.size());
    Complex coeff;
    double norm;

    // Impose <jband | target_band> = 0 on |target_band> (i.e., change |target_band>)
    for (int jband=0; jband<target_band; jband++)
    {
        coeff = vectors_left[jband].dot(target_vector_right); // <jband | target_band>
        target_vector_right -= coeff * vectors_right[jband]; // |target_band> -= <jband | target_band> * |jband>
    }
    
    // Normalize |target_band>
    norm = target_vector_right.norm();
    if (norm < 1e-10) 
    {
        return false;
    }
    else
    {
        target_vector_right /= norm;
    }

    // Impose <target_band | jband> = 0 on <target_band| (i.e., change <target_band|)
    for (int jband=0; jband<target_band; jband++)
    {
        coeff = vectors_right[jband].dot(target_vector_left); // (<target_band | jband>)^*
        target_vector_left -= coeff * vectors_left[jband]; //  <target_band| -= <target_band | jband> * <jband|
    }

    // Impose <target_band | target_band> = 1 on <target_band| (i.e., change <target_band|)
    coeff = target_vector_left.dot(target_vector_right);
    if (std::abs(coeff) < 1e-10)
    {
        return false;
    }
    else
    {
        target_vector_left /= std::conj(coeff);
        return true;
    }
}

// sort eigenvalues as Re(e0) <= Re(e1) <= Re(e2) ...
//   eigenvalues(i) (after sort) =  eigenvalues(sort_index[i]) (before sort)
//   eigenvectors(i,j) (after sort) = eigevectors(i,sort_index[j]) (before sort)
void Diagonalization::sort_eigen(Eigen::VectorXcd &eigenvalues, Eigen::MatrixXcd &eigenvectors) const
{
    const int ndim = eigenvalues.size();
    assert(ndim == eigenvectors.cols() && ndim == eigenvectors.rows());

    // sort eigenvalues & eigenvectors using the real part of eigenvalues
    Complex temp_val;
    for (int idim=0; idim<ndim; idim++)
    {
        for (int jdim=idim+1; jdim<ndim; jdim++)
        {
            // compare eigenvalues(idim) and eigenvalues(jdim)
            if (eigenvalues(idim).real() > eigenvalues(jdim).real()) // swap
            {
                temp_val = eigenvalues(idim);
                eigenvalues(idim) = eigenvalues(jdim);
                eigenvalues(jdim) = temp_val;

                for (int kdim=0; kdim<ndim; kdim++)
                {
                    temp_val = eigenvectors(kdim, idim);
                    eigenvectors(kdim, idim) = eigenvectors(kdim, jdim);
                    eigenvectors(kdim, jdim) = temp_val;
                }
            } // if (swap)
        } // jdim
    } // idim

    // check whether sorting is sucessfully done
    for (int idim=0; idim<ndim-1; idim++)
    {
        if (eigenvalues(idim).real() > eigenvalues(idim+1).real() + 1e-8) 
        {
            error_messages::stop("sorting failed2"); 
        }
    }
}

// sort eigenvalues so that a new eigenvector has the largest overlap with the previous one
void Diagonalization::sort_eigen_overlap(Eigen::VectorXcd &eigenvalues, Eigen::MatrixXcd &eigenvectors) const
{
    const int ndim = eigenvalues.size();
    assert(ndim == eigenvectors.cols() && ndim == eigenvectors.rows());

    Eigen::VectorXcd eigenvalues_temp(ndim);
    Eigen::MatrixXcd eigenvectors_temp(ndim, ndim);

    std::vector<bool> is_jdim_already_chosen(ndim);
    for (int jdim=0; jdim<ndim; jdim++) { is_jdim_already_chosen[jdim] = false; }

    // sort eigenvalues & eigenvectors
    for (int idim=0; idim<ndim; idim++) // determine eigenvalues(idim)
    {
        Complex overlap_save = 0.0;
        int jdim_save = -1; // jdim-eigenvector that has the largest idim-component
        for (int jdim=0; jdim<ndim; jdim++)
        {
            if (is_jdim_already_chosen[jdim]) { continue; } // skipped
            if (std::norm(overlap_save) < std::norm(eigenvectors(idim,jdim))) // idim-component for eigenvalues(jdim)
            {
                overlap_save = eigenvectors(idim,jdim);
                jdim_save = jdim;
            }
        }
        assert(jdim_save >= 0 && jdim_save < ndim);

        // chose jdim_save-eigenvector as a new idim-eigenvector
        is_jdim_already_chosen[jdim_save] = true;
        eigenvalues_temp(idim) = eigenvalues(jdim_save);
        for (int kdim=0; kdim<ndim; kdim++)
        {
            eigenvectors_temp(kdim, idim) = eigenvectors(kdim, jdim_save);
        }
    } // idim
    eigenvalues = eigenvalues_temp;
    eigenvectors = eigenvectors_temp;

    // check whether sorting is sucessfully done
    for (int jdim=0; jdim<ndim; jdim++)
    {
        if (!is_jdim_already_chosen[jdim])
        {
            error_messages::stop("sorting_overlap failed"); 
        }
    }
}

// preconditioning (M. C. Payne, Rev. Mod. Phys.)
void Diagonalization::make_a_new_trial_vector(const Eigen::VectorXd &kGvect2,
                                              Eigen::VectorXcd &phi,
                                              const Eigen::VectorXcd &H1phi,
                                              const Eigen::VectorXcd &H2phi,
                                              const Eigen::VectorXcd &H3phi,
                                              const Complex &eigenvalue,
                                              const bool uses_3body) const
{
    const int num_G_at_k = kGvect2.size();

//    const Complex expectation_value = phi.dot(H1phi + H2phi + H3phi); // <phi|H|phi>

    Eigen::VectorXd phi_norm = phi.array().abs2();
    double kinetic_tot = kGvect2.dot(phi_norm); // <phi|-\nabla^2|phi>
    kinetic_tot /= phi.squaredNorm(); // <phi|-\nabla^2|phi>/<phi|phi> // required for handling left orbitals

    Eigen::VectorXcd coeff_precond(num_G_at_k);
    for (int ipw_at_k=0; ipw_at_k<num_G_at_k; ipw_at_k++) 
    {        
        double k_ratio = kGvect2(ipw_at_k) / kinetic_tot;
        coeff_precond(ipw_at_k) = (27 + k_ratio*(18 + k_ratio*(12 + k_ratio*8)))
            /(27 + k_ratio*(18 + k_ratio*(12 + k_ratio*(8 + k_ratio*16))));
    }
    if (uses_3body)
    {
        phi = - coeff_precond.array() * (H1phi + H2phi + H3phi - eigenvalue * phi).array();
    }
    else
    {
        phi = - coeff_precond.array() * (H1phi + H2phi - eigenvalue * phi).array();
    }

    double norm = phi.norm();
    if (norm < 1e-8)
    {
        phi = Eigen::VectorXcd::Zero(num_G_at_k);
    }
    else
    {
        phi /= norm;
    }
}
