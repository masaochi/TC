Introduction
============

About this code
---------------
TC++ is a first-principles calculation code using the Hartree-Fock (HF) and the transcorrelated (TC) methods for periodic systems.

Supported Functionalities
-------------------------

- Free-electron mode (FREE), HF, TC, BITC (biorthogonal TC)
- SCF and band calculations
- (from ver.1.3) Cell-fixed structural optimization for HF and BITC. Hellmann--Feynman theorem holds for BITC but not for TC.
- Solid-state calculation under the periodic boundary condition. Homogeneous-electron-gas calculation using a periodic cell is also possible.
- Plane-wave basis set
- Norm-conserving pseudopotentials without partial core correction (available, e.g., in `Pseudopotential Library <https://pseudopotentiallibrary.org/>`_)
- For spin-polarized calculation, only spin-collinear calculation without spin-orbit coupling is available.
- Monkhorst-Pack k-grid with/without a shift. A k-grid should not break any crystal symmetry.
- RPA-type Jastrow factor
  :math:`u_{\sigma, \sigma'}({\bf r}, {\bf r'}) = \frac{A_{\sigma, \sigma'}}{|{\bf r}-{\bf r'}|}(1-e^{-|{\bf r}-{\bf r'}|/C_{\sigma,\sigma'}})`
  
Terms of use
------------
TC++ is a free/libre open-source software, released under the MIT License. See the file ``LICENSE`` included in the distribution.
Please cite the following paper in all publications resulting from your use of TC++.

- M. Ochi: "TC++: First-principles calculation code for solids using the transcorrelated method", `Comput. Phys. Commun. 287, 108687 (2023) <https://doi.org/10.1016/j.cpc.2023.108687>`_ .  [`arXiv <https://arxiv.org/abs/2302.07420>`_]

Author & Contact
----------------
`Masayuki Ochi <http://ann.phys.sci.osaka-u.ac.jp/ochi/ochi_en.html>`_ (Osaka University, Japan)  ochi@presto.phys.sci.osaka-u.ac.jp

History
-------
- 2023/12/26 ver.1.3

  + Cell-fixed structural optimization is implemented for HF and BITC. Hellmann--Feynman theorem holds for BITC but not for TC.
  + Gamma-only calculation is supported.
  + Dump ``jastrow.plt`` to display the Jastrow function (can be used by ``load "jastrow.plt"`` with **gnuplot**)
  + MPI+OpenMP parallelization. OpenMP is mainly for memory saving, not so efficient with respect to computational time.
  + Test implementation for polynomial Jastrow functions (Not shown in the manual, will be described in a future release)
  + Default values of the energy and charge convergence criteria are lowered by a factor of 10.
  + Some cpp file names are changed (read_qe.cpp -> io_qe_files_read.cpp etc.)
  + Bug fix for spin-polarized calculation (fixed in ver.1.2.2)
  + Bug fix for divergence correction for TC and BITC (fixed in ver.1.2.3)
    
- 2023/1/16 ver.1.2

  + cmake and test are implemented.
  + A bug in reading upf files is fixed. Some upf files were not readable in older versions of TC++.
  + Some include files were not properly updated in ver.1.1.

- 2022/12/14 ver.1.1

  + ``mixes_density_matrix`` is implemented, which can improve convergence of calculation (See :doc:`input_in`).
  + A format of the total energies in ``output.out`` is slightly changed in order to track the total-energy convergence easily (e.g., by ``p `< grep "Total energy =" output.out' u 5`` with **gnuplot**).
  + Some files that are not used are removed (e.g., ``calc_hamiltonian_tc3a2.cpp``).

- 2022/11/18 ver.1.0
