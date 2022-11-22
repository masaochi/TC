Introduction
============

About this code
---------------
TC++ is a first-principles calculation code using the Hartree-Fock (HF) and the transcorrelated (TC) methods for periodic systems.

Supported Functionalities
-------------------------

- Free-electron mode (FREE), HF, TC, BITC (biorthogonal TC)
- SCF and band calculations
- Solid-state calculation under the periodic boundary condition. Homogeneous-electron-gas calculation using a periodic cell is also possible.
- Plane-wave basis set
- Norm-conserving pseudopotentials without partial core correction (available, e.g., in `Pseudopotential Library <https://pseudopotentiallibrary.org/>`_)
- For spin-polarized calculation, only spin-collinear calculation without spin-orbit coupling is available.
- Monkhorst-Pack k-grid with/without a shift. A k-grid should not break any crystal symmetry. Gamma-only calculation is at present not supported.
- RPA-type Jastrow factor

Terms of use
------------
TC++ is a free/libre open-source software, released under the MIT License. See the file ``LICENSE`` included in the distribution.

Author & Contact
----------------
`Masayuki Ochi <http://ann.phys.sci.osaka-u.ac.jp/ochi/ochi_en.html>`_ (Osaka University, Japan)  ochi@presto.phys.sci.osaka-u.ac.jp
