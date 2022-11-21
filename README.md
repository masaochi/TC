# TC++
Version 1.0 (2022/11/18)

TC++ is a first-principles calculation code using the Hartree-Fock (HF) and the transcorrelated (TC) methods for periodic systems.

## Supported functionalities
- Free-electron mode (FREE), HF, TC, BITC (biorthogonal TC)
- SCF and band calculations
- Solid-state calculation under the periodic boundary condition. Homogeneous-electron-gas calculation using a periodic cell is also possible.
- Plane-wave basis set
- Norm-conserving pseudopotentials without partial core correction (available, e.g., in [Pseudopotential Library](https://pseudopotentiallibrary.org/))
- Non-spin-polarized calculation or spin-polarized calculation. For the latter one, only spin-collinear calculation without spin-orbit coupling is available.
- Monkhorst-Pack k-grid with/without a shift. k-grid should not break any crystal symmetry. Gamma-only calculation is at present not supported.
- RPA-type Jastrow factor

## Installation
Download the source files and unzip it. Then, `cd src` and edit Makefile to specify the following compilers and libraries. Finally, typing `make` will create an execution file named `tc++` in `src`.

### Prerequisite
- C++ compiler (C++11 or newer)
- Fortran compiler (Fortran90 or newer)
- MPI library
- Boost C++ library
- FFTW library
- Eigen3 library
- Quantum ESPRESSO (ver. 6.2 or newer) is used in precalculation.

## Documentation
Under construction... (available soon!)

## License
Copyright (c) 2022 Masayuki Ochi

Released under the MIT license. See LICENSE file.

## Author
Masayuki Ochi (Osaka University)

## Papers
Under construction...


