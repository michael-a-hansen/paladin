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
    -D CMAKE_INSTALL_PREFIX=./install/share \
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

to compute the eigenvalues of several simple test matrices and compare to exact eigenvalues.
The tests cover calculations performed in serial, and in parallel with two and four cores.

## Running a calculation
To run a calculation you need three things:
1. a successful install build
2. matrix files
3. a matrix listing

### example matrix file
The triangular matrix below has the eigenvalues `{1+0j, 3+0j, 9+0j}`
```
A = [1, 2, 7;
     0, 3, 4;
     0, 0, 9]
```
The file in matrix market form appropriate for paladin is below.
The first line is necessary in that exact form.
The following lines starting with `%` are comments that have no impact on paladin.
The first line following the comments, `3 3 5` is the number of rows, number of columns, and number of nonzero elements.
For eigendecompositions, the matrix must be square and so the numbers of rows and columns must be equivalent.

Following this are the matrix elements, written as `[row index] [column index] [element value]` where the indices are 1-based. Compare the matrix written above and below until it is clear.

```
%%MatrixMarket matrix coordinate real general
%
% comments...
% example triangular matrix
%
3 3 6
1 1 1
1 2 2
1 3 7
2 2 3
2 3 4
3 3 9
```

### example matrix listing
With the matrix file above, we move to the matrix listing.
Suppose we called the matrix file `tri-example1.matrix`.
Our listing file would simply include this file name on its own line, followed by any other matrix files to be included in the calculation.

```
tri-example1.matrix
tri-example2.matrix
tri-example3.matrix
...
```


### running paladin
#### inside the matrix directory
Suppose the listing has the path, `/Users/mike/paladin-example/listing` while the matrix files are `/Users/mike/paladin-example/tri-example1.matrix`, `/Users/mike/paladin-example/tri-example2.matrix`, `/Users/mike/paladin-example/tri-example3.matrix`, etc.

From inside the `/Users/mike/paladin-example` directory we can run paladin in serial with

`[...]/paladin/build/src/exec-paladin --listing=listing`,

and in parallel on four cores with

`mpirun -np 4 [...]/paladin/build/src/exec-paladin --listing=listing`.

#### outside the matrix directory
Suppose we want to run only the `tri-example1.matrix` and `tri-example2.matrix` matrices.
We would make a new listing file,

```
tri-example1.matrix
tri-example2.matrix
```

and perhaps we saved this new listing file in a new directory `/Users/mike/other-listing/listing2`. There is a problem now. If we run paladin from this directory it will look for `tri-example1.matrix` in the _current_ directory, where it is _not_.

To solve this issue, we could move the listing to the old directory and run paladin from there... or we could run paladin from any directory by using the `rootdir` option and specifying a full path on the listing, as below. The `rootdir` is the location at which the matrices exist. When `rootdir` is specified, paladin looks for a matrix at `[rootdir]/tri-example1.matrix`.

```
[...]/paladin/build/src/exec-paladin \
     --listing=/Users/mike/other-listing/listing2 \
     --rootdir=/Users/mike/paladin-example
```


## Summary of command line options
paladin command line options must be given as `--key=value`, as seen above in the examples. Below all of the options are summarized:

1. `--listing=[value]`: the path of the matrix listing
2. `--rootdir=[value]`: prefix to the path of the matrix files in the listing
3. `--showdist`: add this to the command line to show the distribution of matrix files to each MPI rank
4. `--repeats=[value]`: the number of times the decomposition is repeated. Use this if you are assessing performance and want to average the cost of the eigendecomposition over a number of runs.
5. `--load-measure=[value]`: the measure of the matrix used for load balancing in parallel calculations. Options are below. `nnz` does a good job in most cases. If you have a large variety in matrix sizes, `dcb` may do well. For very dense matrices, I/O time is nontrivial and `nnz` should be used.
 - `nnz`: number of nonzeros in the matrix (default load measure)
 - `dim`: dimension of the matrix (number of rows)
 - `dcb`: dimension cubed
 - `zds`: number of nonzeros multiplied by the dimension squared
 - `sps`: sparsity - number of nonzeros over dimension squared
 - `spc`: sparsity multiplied by dimension cubed - nnz * dim


# The name

paladin = **PA**rallel **LA**pack **DI**stributor + **n** :)