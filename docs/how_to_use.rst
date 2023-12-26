How to use
==========

TC++ calculation is performed after precalculation using Quantum Espresso (QE).
If you prefer Quick Start of TC++, skip this section and go to :doc:`tutorial`.

QE precalculation
-----------------

Before TC++ calculation, users should perform QE calculation using DFT, HF, etc.
Any exchange-correlation energy functionals can be used as long as one can get one-electron orbitals.
Please note the following issues in QE calculation:

- It is recommended to perform QE calculation in the same environment for Fortran as TC++ because TC++ reads binary files containing wave-function data dumped by QE. In TC++, Fortran is used only for this purpose.

- Set ``wfcdir`` to be the same as ``outdir`` to ensure that the xml and wfc files (see below) are dumped in the same directory.
  Since this is the default setting for ``wfcdir`` in QE, there is no need to specify ``wfcdir`` in the input file, at present.

- Please perform BAND calculation in a different directory from SCF calculation. To say, users should first perform SCF calculation,
  copy that directory, and perform BAND calculation there. See :doc:`tutorial` for details.

- Please specify ``noinv = .true.`` and ``no_t_rev = .true.`` for spin-polarized calculation since TC++ at present does not support these symmetry operations.
  
- There are some restrictions on pseudopotentials. See below.


Input files of TC++
-------------------

After QE calculation, users should prepare the following input files for TC++.
Users just need to make ``input.in`` and specify where the first three files (pseudopotentials, xml, and wfc files) are placed.

- Pseudopotential files (upf format)

  Norm-conserving pseudopotentials without partial core correction are acceptable: available, e.g., in `Pseudopotential Library <https://pseudopotentiallibrary.org/>`_.
  Both of the old and new formats of upf are acceptable.
  The pseudopotential files for TC++ should be the same with those used in QE calculation.

- ``data-file-schema.xml``

  This file is dumped by QE, saved in ``outdir``/ ``prefix`` *.save* directory as specified in the QE input file.
  TC++ reads several information from this file.

- ``wfc.dat`` (e.g., ``wfc1.dat``, ``wfc2.dat``, etc.)

  These file are dumped by QE, saved in ``outdir``/ ``prefix`` *.save* directory.
  TC++ reads one-electron orbitals obtained by QE from this file.

- ``input.in``

  This file should be made by users and read by TC++. See :doc:`input_in` for details.

For restarting SCF calculation or performing band calculation after SCF, TC++ requires some other input files dumped by TC++. Please see output_files_.
  
How to run TC++
---------------

Perform TC++ calculation with a command, e.g.,
::

   $ mpirun -np 4 $HOME/TC++/ver.1.0/src/tc++

One can also use MPI+OpenMP parallelization by setting **OMP_NUM_THREADS** to a non-unity value, while MPI is more efficient than OpenMP with respect to computational time in TC++.

.. note::

   If you find that the memory requirement is too demanding (e.g., in massively parallel computation), MPI+OpenMP calculation with a large value of **OMP_NUM_THREADS**
   can reduce the memory requirement **per node**.

.. _output_files:


Output files of TC++
--------------------

The following output files are obtained by TC++ calculation:

- Standard output
  
  **Error messages are shown in the standard output. Please check it when calculation unexpectedly stops.**

- ``output.out``

  Many information obtained by calculation are shown. For example, you can track the total-energy convergence by
  ``p `< grep "Total energy =" output.out' u 5`` with **gnuplot**.

- ``tc_bandplot.dat``

  Band eigenvalues are shown when ``calc_mode = BAND``. For example,
  ``p 'tc_bandplot.dat' u 4:5 w l``
  with **gnuplot** will show the band structure. The Fermi energy obtained by SCF calculation is also shown in this file.
  These information are also shown in ``output.out``.
  For spin-polarized calculation, ``tc_bandplot_up.dat`` and ``tc_bandplot_dn.dat`` are dumped instead.

- ``jastrow.plt`` (from ver.1.3)

  A Jastrow function used in calculation is shown when ``calc_method = TC`` or ``calc_method = BITC``.
  ``load 'jastrow.plt'``
  with **gnuplot** will show the spin-parallel and spin-antiparallel Jastrow functions (in atomic unit).

- ``tc_crystal_structure.dat`` (from ver.1.3)

  An optimized crystal structure is dumped after structural optimization.
  This file can be read for further (subsequent) structural optimization calculation by setting ``reads_crystal_structure = true`` (i.e., users can restart the structural optimization from the new structure).
  Note that the structural optimization is switched on only when ``calc_mode = SCF`` and ``is_heg = false`` (default) and  ``calc_method = HF or BITC`` and ``max_num_ionic_steps > 0``.
  Thus, if you would like to perform SCF and BAND calculations after the structural optimization, please rerun QE with the optimized crystal structure because
  SCF nor BAND calculation with ``max_num_ionic_steps = 0`` (default) cannot read ``tc_crystal_structure.dat``. This is because rerunning QE is safer in accuracy considering that an initial guess of the orbitals is given by QE.
  
The following binary files are dumped in TC++ calculation.
Users will not read them but some subsequent TC++ calculations need them.

- ``tc_energy_scf.dat``

  SCF energy eigenvalues that are used for restarting SCF calculation or performing subsequent BAND calculation. Dumped in SCF calculation.

- ``tc_energy_band.dat``

  BAND energy eigenvalues that are used for restarting BAND calculation. Dumped in BAND calculation.

- ``tc_wfc_scf.dat``

  SCF wave functions that are used for restarting SCF calculation or performing subsequent BAND calculation. Dumped in SCF calculation.

- ``tc_wfc_band.dat``

  BAND wave functions that are used for restarting BAND calculation. Dumped in BAND calculation.


- ``tc_scfinfo.dat``

  Several information of SCF calculation that are used for subsequent BAND calculation.

Here, ``tc_energy_scf/band.dat`` and ``tc_wfc_scf/band.dat`` are dumped in each self-consistent loop so that users can restart calculation when calculation stops.



