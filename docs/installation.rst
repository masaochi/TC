Installation
============

Download
--------
The computational code is available in github
https://github.com/masaochi/TC

Prerequisites
-------------

- C++ compiler (C++11 or newer)
- Fortran compiler (Fortran90 or newer)
- MPI library
- Boost C++ library
- FFTW library (Note: if you compiled FFTW with the Intel compiler, please compile TC++ with the Intel compiler to avoid some errors.)
- Eigen3 library
- Quantum ESPRESSO (ver.6.2 or newer) is used in precalculation to get several information such as an initial guess of one-electron orbitals. TC++ requires xml and wfc files dumped by QE.




