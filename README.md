# paladin
A code for running many serial LAPACK eigenvalue solves distributed across MPI ranks on multi-core CPU systems. It is written in C++11 and uses CMake. Aside from BLAS and LAPACK, no external libraries are necessary.

paladin is meant to provide a means to compute a large number (say, tens of thousands) of eigenvalue decompositions of small matrices (say, thousands of rows) using LAPACK in serial on a multi-core system.
Load balancing across MPI ranks is performed based on matrix sizes and sparsity.
Several load measures are provided to enable optimal performance across the range of sparser, larger systems (eigendecomposition dominates cost) and denser, smaller systems (matrix I/O dominates cost).

Matrices must be provided in the matrix market exchange format of NIST (see http://math.nist.gov/MatrixMarket/formats.html).
At the moment only general (dense), real, nonsymmetric (nothing special is done for symmetric) matrices are supported, although the extension to symmetric systems, for instance, would be straightforward.

paladin = PArallel LApack DIstributor