# paladin
A C++ code for running many _serial_ LAPACK eigendecompositions distributed across MPI ranks on multi-core CPU systems. paladin depends on C++11, CMake, BLAS, and LAPACK. No additional external libraries are necessary.

This provides a portable, lightweight means of computing a large number (say, tens of thousands) of decompositions of small matrices (say, thousands of rows) using LAPACK in serial on a multi-core system. Up to 128 cores have been used with paladin. The benefit obtained through paladin is efficiency on typical clusters, where a huge number of 1-core jobs would progress through the queue _far_ slower than a small number of larger jobs. For instance, to compute 12,800 decompositions one could request 12,800 1-core jobs or a single 128-core job with 100 matrices per core. The latter solution is likely to complete far faster.

For problems with variety in matrix structure, static load balancing across MPI ranks is performed based on matrix sizes and sparsity.
Several load measures are provided to enable optimal performance across the range of sparser, larger systems (eigendecomposition dominates cost) and denser, smaller systems (I/O dominates cost).

Matrices must be provided in NIST's matrix market exchange format (see http://math.nist.gov/MatrixMarket/formats.html).
At the moment only general (dense), real, general matrices written in coordinate form are supported.

paladin has been used on the following systems and compilers:

- Mac OS X, Apple LLVM version 7.0.2 (clang-700.1.81), MPI configured with gcc 4.2.1
- gcc 4.7.2 and openmpi 1.6.4 on a Red Hat 6 workstation
- Intel 14.0 and openmpi 1.6.4 on the NNSA Tri-Lab Linux capacity clusters (note: the `std=c++11` flag must be explicitly added here...).

C++11 demands CMake 3.1, and paladin has been built with 3.7.1 on Mac OS X and 3.4.3 on Linux.

As for the name... paladin = **PA**rallel **LA**pack **DI**stributor + **n**