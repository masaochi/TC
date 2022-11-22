# TC++
Version 1.0 (2022/11/18)

TC++ is a first-principles calculation code using the Hartree-Fock (HF) and the transcorrelated (TC) methods for periodic systems.

## Supported functionalities
- Free-electron mode (FREE), HF, TC, BITC (biorthogonal TC)
- SCF and band calculations
- Solid-state calculation under the periodic boundary condition. Homogeneous-electron-gas calculation using a periodic cell is also possible.
- Plane-wave basis set
- Norm-conserving pseudopotentials without partial core correction (available, e.g., in [Pseudopotential Library](https://pseudopotentiallibrary.org/))
- For spin-polarized calculation, only spin-collinear calculation without spin-orbit coupling is available.
- Monkhorst-Pack k-grid with/without a shift. A k-grid should not break any crystal symmetry. Gamma-only calculation is at present not supported.
- RPA-type Jastrow factor

## Installation
Download the source files and unzip it. Then, `cd src` and edit [`Makefile`](./src/Makefile) to specify the following compilers and libraries except Quantum ESPRESSO. Finally, typing `make` will create an execution file named `tc++` in `src`.

### Prerequisites
- C++ compiler (C++11 or newer)
- Fortran compiler (Fortran90 or newer)
- MPI library
- [Boost C++ library](https://www.boost.org/)
- [FFTW library](https://www.fftw.org/) (Note: If you compiled FFTW with the Intel compiler, please compile TC++ with the Intel compiler to avoid some errors.)
- [Eigen3 library](https://eigen.tuxfamily.org/)
- [Quantum ESPRESSO](https://www.quantum-espresso.org/) (ver.6.2 or newer) is used in precalculation to get several information such as an initial guess of one-electron orbitals. TC++ requires xml and wfc files dumped by QE.

## Documentation
Under construction... (available soon!)

## License
Copyright (c) 2022 Masayuki Ochi

Released under the MIT license. See [`LICENSE`](./LICENSE).

## Author & Contact
[Masayuki Ochi](http://ann.phys.sci.osaka-u.ac.jp/ochi/ochi_en.html) (Osaka University, Japan)
ochi@presto.phys.sci.osaka-u.ac.jp



