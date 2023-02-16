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
- `Quantum ESPRESSO <https://www.quantum-espresso.org/>`_ (QE) (ver.6.2 or newer) is used in precalculation to get several information such as an initial guess of one-electron orbitals. TC++ requires xml and wfc files dumped by QE.

Building 1 (directly editing Makefile)
--------------------------------------

There are two ways for installation. One is directly editing ``Makefile``.
Unzip the downloaded compressed source file, type

::

   $ cd src

and edit ``Makefile`` to specify the compilers and the libraries listed in :ref:`label_prerequisites` (except QE). Then, typing

::

   $ cd make

will create an execution file named ``tc++`` in ``src``.

Building 2 (cmake)
------------------

The other way for installation is using cmake. Unzip the downloaded compressed source file, and type

::

   $ mkdir build && cd build
   $ cmake ..
   $ make
   $ make install

will create an execution file named ``tc++``. Make sure that ``CMakeLists.txt`` exists in ``../`` when you type ``cmake ..``.
Here, you might need to specify some options for ``cmake`` when some compilers or libraries are not properly set by ``cmake``. Available options are shown below:

- CMAKE_CXX_COMPILER= *(C++ compiler name)*
- CMAKE_Fortran_COMPILER= *(Fortran90 comipler name)*
- EIGEN3_INCLUDE= *(Eigen3 directory)*
- FFTW_INCLUDE= *(FFTW include directory)*
- FFTW_LIB = *(FFTW library directory)*
- BOOST_INCLUDE= *(BOOST include directory)*
- BOOST_LIB = *(BOOST library directory)*
- CMAKE_INSTALL_PREFIX = *(where to install tc++)*
  
For example:

::

   $ cmake -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DEIGEN3_INCLUDE=/usr/local/eigen-3.4.0 -DFFTW_INCLUDE=/usr/local/fftw-3.3.10/include -DFFTW_LIB=/usr/local/fftw-3.3.10/lib -DCMAKE_INSTALL_PREFIX=/home/ochi/TC++ ..
  
will create ``tc++`` in ``/home/ochi/TC++/bin`` (NOTE: ``bin`` is added to CMAKE_INSTALL_PREFIX).
Note that MPI is automatically searched by ``cmake``, so please do not specify, e.g., mpiicpc, for CMAKE_CXX_COMPILERS.

Test
----

To verify that your installation was successfully done, a test suite is provided in the ``test`` directory. Type

::

   $ cd test

and copy ``tc++`` compiled above to the ``test`` directory. Then, test calculation will start by typing

::

   $ python3 test.py

, which will take one or few minutes.

.. note::

   python2 is not supported.

.. note::

   DO NOT use the input files (including pseudopot.) provided here for your research. They were made only for test calculation.

.. note::

   ``test_full.py`` is also provided to cover a wider range of functionalities. This test is rather for developers, which will take 10--20 minutes.

.. note::

   ``reads_binary false`` is specified in ``input.in`` in these test suites.
   This option is used only for test calculation, where environment-dependent binary files are not appropriate.
   You will not use this option in your research.




