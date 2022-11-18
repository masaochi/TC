// [class Diagonalization]
// scf_loop, convergence check, control paramters, etc.

#include "include/header.hpp"

// make target_vector orthonormalized to vectors[i] (0 <= i < target_band)
bool Diagonalization::Gram_Schmidt(std::vector<Eigen::VectorXcd> &vectors, const int &target_band,
                                   Eigen::VectorXcd &target_vector)
{
    assert(target_band < vectors.size());
    Eigen::VectorXcd coeff(target_band);

    for (int iband=0; iband<target_band; iband++) 
    {
        coeff(iband) = vectors[iband].dot(target_vector); // <iband | target_vector>
        // If you comment-out the below two lines, unmodified Gram Schmidt is performed.
        //  }
        //  for (int iband=0; iband<target_band; iband++) {
        target_vector -= coeff(iband) * vectors[iband]; // |target_vector> -= <iband | target_vector> * |iband>
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

// sort eigenvalues as Re(e0) <= Re(e1) <= Re(e2) ...
//   eigenvalues(i) (after sort) =  eigenvalues(sort_index[i]) (before sort)
//   eigenvectors(i,j) (after sort) = eigevectors(i,sort_index[j]) (before sort)
void Diagonalization::sort_eigen(Eigen::VectorXcd &eigenvalues, Eigen::MatrixXcd &eigenvectors)
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

// preconditioning (M. C. Payne, Rev. Mod. Phys.)
void Diagonalization::make_a_new_trial_vector(const Eigen::VectorXd &kGvect2,
                                              Eigen::VectorXcd &phi,
                                              const Eigen::VectorXcd &H1phi,
                                              const Eigen::VectorXcd &H2phi,
                                              const Eigen::VectorXcd &H3phi,
                                              const Complex &eigenvalue,
                                              const bool uses_3body)
{
    const int num_G_at_k = kGvect2.size();

//    const Complex expectation_value = phi.dot(H1phi + H2phi + H3phi); // <phi|H|phi>

    Eigen::VectorXd phi_norm = phi.array().abs2();
    double kinetic_tot = kGvect2.dot(phi_norm); // <phi|-\nabla^2|phi>

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
