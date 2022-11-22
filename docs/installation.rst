Installation
============

Download
--------
The computational code is available in github
https://github.com/masaochi/TC

.. _label_prerequisites:

Prerequisites
-------------

- C++ compiler (C++11 or newer)
- Fortran compiler (Fortran90 or newer)
- MPI library
- `Boost C++ library <https://www.boost.org/>`_
- `FFTW3 library <https://www.fftw.org/>`_ (Note: if you compiled FFTW with the Intel compiler, please compile TC++ with the Intel compiler to avoid some errors.)
- `Eigen3 library <https://eigen.tuxfamily.org/>`_
- `Quantum ESPRESSO <https://www.quantum-espresso.org/>`_ (ver.6.2 or newer) is used in precalculation to get several information such as an initial guess of one-electron orbitals. TC++ requires xml and wfc files dumped by QE.

Building
--------

Unzip the downloaded source files, type

  cd src

and edit ``Makefile`` to specify the compilers and the libraries listed in :ref:`label_prerequisites` (except Quantum ESPRESSO). Then, typing

  make

will create an execution file named ``tc++`` in ``src``.

