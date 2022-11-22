Tutorial
========

Bulk Silicon
------------

First, perform SCF calculation by QE using the following input file:

::

   &control
     prefix = 'prefix'
     calculation = 'scf'
     pseudo_dir = '/home/user/QE/pseudo_potential/'
     outdir = './'
     verbosity = 'high'
     disk_io = 'low'
   /
   &system
     ibrav = 2
     celldm(1) = 10.26
     nat = 2
     ntyp = 1
     nbnd = 10
     ecutwfc = 20.0
     occupations = 'fixed'
   /
   &electrons
     conv_thr = 1.0d-8
   /
   ATOMIC_SPECIES
     Si 1.0 Si.upf
   ATOMIC_POSITIONS {alat}
     Si 0.00 0.00 0.00
     Si 0.25 0.25 0.25
   K_POINTS {automatic}
      8 8 8 0 0 0

Here, the Ne-core pseudopotentials of silicon [1]_ taken from `Pseudopotential Library <https://pseudopotentiallibrary.org/>`_ is used.
Any calculation method in QE, such as DFT and HF, is acceptable as long as one can get one-elecron orbitals.
To obtain a band structure, we also need to perform band calculation using QE.
For this purpose, copy the directory where SCF calculation was performed and perform band calculation there.
Namely, SCF and band calculations should be performed in different directories.
Band calculation in QE is performed with the following input file:

::

   &control
     prefix = 'prefix'
     calculation = 'bands'
     pseudo_dir = '/home/user/QE/pseudo_potential/'
     outdir = './'
     verbosity = 'high'
     disk_io = 'low'
   /
   &system
     ibrav = 2
     celldm(1) = 10.26
     nat = 2
     ntyp = 1
     nbnd = 10
     ecutwfc = 20.0
     occupations = 'fixed'
   /
   &electrons
     conv_thr = 1.0d-8
   /
   ATOMIC_SPECIES
     Si 1.0 Si.upf
   ATOMIC_POSITIONS {alat}
     Si 0.00 0.00 0.00
     Si 0.25 0.25 0.25
   K_POINTS {crystal_b}
   3
    0.5 0.5 0.0 20
    0.0 0.0 0.0 20
    0.5 0.0 0.0 0

Next, perform SCF calculation using the following input file, ``input.in``,

::

   calc_method   HF  # change here for different methods (TC, BITC)
   calc_mode     SCF
   pseudo_dir    /home/user/QE/pseudo_potential
   qe_save_dir   /home/user/where_QE_SCFcalc_was_performed/prefix.save
   smearing_mode fixed

where ``pseudo_dir`` and ``qe_save_dir`` shoud be appropriately specified.
After SCF calculation, we should check whether **convergence is achieved!** is shown in ``output.out``.
If the convergence is not achieved, we can restart calculation using ``input.in`` where the following line is added:

::

   restarts  true

However,

.. [1] M. Chandler Bennett *et al.*, J. Chem. Phys. **149**, 104108 (2018).

Homogeneous Electron Gas
------------------------



