# paladin
A C++ code for running many _serial_ LAPACK eigendecompositions distributed across MPI ranks on multi-core CPU systems. paladin depends on C++11, CMake, MPI, BLAS, and LAPACK. No additional external libraries are necessary.

This provides a portable, lightweight means of computing a large number (say, tens of thousands) of decompositions of small matrices (say, thousands of rows) using LAPACK in serial on a multi-core system. Up to 128 cores have been used with paladin. The benefit obtained through paladin is efficiency on typical clusters, where a huge number of 1-core jobs would progress through the queue _far_ slower than a small number of larger jobs. For instance, to compute 12,800 decompositions one could request 12,800 1-core jobs or a single 128-core job with 100 matrices per core. The latter solution is likely to complete far faster.

For problems with variety in matrix structure, static load balancing across MPI ranks is performed based on matrix sizes and sparsity.
Several load measures are provided to enable optimal performance across the range of sparser, larger systems (eigendecomposition dominates cost) and denser, smaller systems (I/O dominates cost).

Matrices must be provided in NIST's matrix market exchange format (see http://math.nist.gov/MatrixMarket/formats.html).
At the moment only general (dense), real, general matrices written in coordinate form are supported.

paladin has been used on the following systems and compilers:

- Mac OS X, Apple LLVM version 7.0.2 (clang-700.1.81), MPI configured with gcc 4.2.1
- gcc 4.7.2 and openmpi 1.6.4 on a Red Hat 6 workstation
- Intel 14.0 and openmpi 1.6.4 on the NNSA Tri-Lab Linux capacity clusters (note: the `std=c++11` flag must be explicitly added here...).

paladin has been built with 3.7.1 on Mac OS X and 3.4.3 on Linux. The C++11 requirement means at least CMake 3.1 is required.

# Usage
## Installation

To use paladin first ensure that you have the appropriate dependencies installed, namely CMake 3.1, a C++11-compliant compiler, and BLAS and LAPACK. It is very likely that CMake is the only code you'll have to install or upgrade manually.

Then clone the repository, 

`git clone https://github.com/michael-a-hansen/paladin.git`,

enter the directory, make and enter a build directory,

`cd paladin; mkdir build; cd build`,

### with default configuration
To install paladin with default settings and allow system discovery of cmake, compilers, etc. simply run `cmake`:

`cmake ..`,

and then

`make -j install`

### specifying compilers and flags
To configure paladin specially for your system, which may be necessary if cmake discovers a compiler executable you don't want to use, or if you want a particular cmake executable, make an executable configuration file in your build directory.

`touch configme; chmod +x configme`

A configuration file might look like the following. `[...]` indicates a prefix to the rest of the path. For instance a compiler executable's path might be `/usr/bin/mpic++`. Below we write `[...]/bin/mpic++` so `[...]` would be `/usr`. The purpose of the `EXTRA_ARGS` statement is so that configme can be configured on the fly when it is executed, as below.

```
EXTRA_ARGS=$@
 
/[...]/bin/cmake-3.4.3/bin/cmake \
    -D CMAKE_CXX_COMPILER=/[...]/bin/mpic++ \
    -D CMAKE_C_COMPILER=/[...]/bin/mpic \
    -D CMAKE_CXX_FLAGS=-std=c+=11 \
    $EXTRA_ARGS \
    ..
```

Instead of running `cmake ..` as above, we execute our config file:

`./configme`

We might run this as `./configme -D CMAKE_BUILD_TYPE=DEBUG`, in which case the extra argument is passed directly into the config file's execution of `cmake`. Observe that the configuration file is really just storing a big call to the same `cmake` executable we called earlier without arguments.

After the configuration is successful, proceed in the same way:

`make -j install`.


## Running the tests
After the code has been installed (a successful call to `make install`), simply navigate to the build directory and run

`ctest`

to compute the eigenvalues of several simple test matrices and compare to exact results or gold standards.
The tests cover calculations performed in serial, and in parallel with two and four cores.
Small, simple matrices are compared against exact eigenvalues while several Jacobian matrices from combustion simulations are compared against gold standards.

## Running a calculation


## Command line options

# The name

paladin = **PA**rallel **LA**pack **DI**stributor + **n** :)