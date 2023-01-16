# TC++
Version 1.1 (2022/12/14)

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
Download the compressed source file and unzip it (see Releases or Tags). Then, there are two ways for installation:

(1) `cd src` and edit [`Makefile`](./src/Makefile) to specify the following compilers and libraries except Quantum ESPRESSO. Finally, typing `make` will create an execution file named `tc++` in `src`.

(2) `cmake` is also available for installation. Type `mkdir build && cd build`, `cmake ..`, `make`, and `make install` to create an execution file named `tc++`. For several options for `cmake`, please see [User's Guide](https://TCplusplus.readthedocs.io/).

### Prerequisites
- C++ compiler (C++11 or newer)
- Fortran compiler (Fortran90 or newer)
- MPI library
- [Boost C++ library](https://www.boost.org/)
- [FFTW3 library](https://www.fftw.org/) (Note: If you compiled FFTW with the Intel compiler, please compile TC++ with the Intel compiler to avoid some errors.)
- [Eigen3 library](https://eigen.tuxfamily.org/)
- [Quantum ESPRESSO](https://www.quantum-espresso.org/) (ver.6.2 or newer) is used in precalculation to get several information such as an initial guess of one-electron orbitals. TC++ requires xml and wfc files dumped by QE.

To verify that your installation is successful, a test suite is provided.
Type `cd test` and copy `tc++` into `test` directory. Then, you can perform a test calculation by typing `python3 test.py`.

## Documentation
[User's Guide](https://TCplusplus.readthedocs.io/)

## License
Copyright (c) 2022--2023 Masayuki Ochi

Released under the MIT license. See [`LICENSE`](./LICENSE).

## Author & Contact
[Masayuki Ochi](http://ann.phys.sci.osaka-u.ac.jp/ochi/ochi_en.html) (Osaka University, Japan)
ochi@presto.phys.sci.osaka-u.ac.jp



