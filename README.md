# paladin
A code for running many serial LAPACK eigenvalue solves distributed across MPI ranks on multi-core CPU systems with automatic load balancing


paladin is meant to provide a means to compute many (say, tens of thousands) of eigenvalue decompositions of small (say, tens of thousands of rows) matrices using LAPACK in serial on a multi-core system.
Load balancing across MPI ranks is performed based on matrix sizes and sparsity.
Matrices must be provided in the matrix market exchange format of NIST (see http://math.nist.gov/MatrixMarket/formats.html).
At the moment only general (dense), real, nonsymmetric (nothing special is done for symmetric) matrices are supported, although the extension to symmetric systems, for instance, would be straightforward.
Several load measures are provided to provide optimal performance across the range of sparse, large systems (eigendecomposition dominates cost) and dense, small systems (matrix I/O dominates cost).
