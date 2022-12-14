Tips & FAQ
==========


Tips
----

- How the computational cost scales?

  Computational time :math:`\propto N_k^2 N_b^2 N_{pw}\log N_{pw}`, memory requirement :math:`\propto N_k N_b N_{pw}`,
  where :math:`N_k, N_b, N_{pw}` are the numbers of k-points, bands, and plane waves, respectively.
  This is the same irrespective of the methods (HF, TC, BITC).

- Convergence of (BI)TC calculation is not good.

  Possible ways to resolve the convergence problem are as follows: (i) increase the number of k-points, (ii) increase ``max_num_blocks_david`` (e.g., to 5),
  (iii) increase the number of bands, and (iv) reduce ``mixing_beta`` with ``mixes_density_matrix`` = true. In general, (i) is the most effective in many cases because of the divergence correction in the reciprocal space.
  Since ``A_up_up = A_up_dn = A_dn_dn = 0.0`` in (BI)TC is equivalent to HF, performing (BI)TC calcultaion with a small value of the A parameter,
  and then restarting calculation with increasing the A values, might be another way to get convergence if HF has no problem for convergence.
  
- Band dispersion is not smooth.

  Due to the difficulty in correcting  the :math:`1/(k-k')^2` divergence of the electron-electron interactions in the reciprocal space,
  eigenvalues of the band k-points (:math:`k'`) that are too close to the SCF k-points (:math:`k`) often have a large error.
  We recommend that users should not take such a k-point in BAND calculation.
  If many band k-points have a large error, please increase the number of k-points in SCF calculation.
  A larger number of k-points compared with GGA might be necessary in TC calculation.

- Memory requirement is too demanding.

  Increasing the number of MPI processes can alleviate this issue by distributing large arrays to many MPI processes.
  In addition, using non-unity **OMP_NUM_THREADS** can also be helpful. While this is not effective for acceleration since TC++ does not explicitly use
  OpenMP parallelization, memory requirement **per node** will be reduced in massively-parallelized calculation.


FAQ
---

will be added when we give questions from users!
