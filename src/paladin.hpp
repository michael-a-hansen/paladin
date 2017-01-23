/*
 * The MIT License
 *
 * Copyright (c) 2016, 2017 Michael A. Hansen
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef SRC_PALADIN_HPP_
#define SRC_PALADIN_HPP_

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <command-line-parser.hpp>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <header.hpp>
#include <iomanip>
#include <iostream>
#include <lapack-wrapper.hpp>
#include <mpi-util.hpp>
#include <spectrum-util.hpp>
#include <square-matrix.hpp>
#include <sstream>
#include <string>
#include <vector>

namespace paladin {

/**
 * @class Paladin
 * @brief PArallel LApack DIstributor for many eigenvalue calculations
 */
class Paladin {
 protected:
  std::vector<std::string> allMatrixPaths, pod;
  std::vector<int> podColors;
  std::vector<SpectralDecomposition> spectra_;
  double localRunTime_, localReadRunTime_, localEigRunTime_;
  std::string matrixListingPath_, matrixRootPath_;
  int numRuns_, numMatrices_;
  bool showPods_;
  LoadPredictor type_;
  const MpiComm& comm_;
  bool doRight_ = false;
  bool doLeft_ = false;

 public:
  Paladin( int& argc, char* argv[], const MpiComm& comm )
      : localRunTime_( 0.0 ),
        localReadRunTime_( 0.0 ),
        localEigRunTime_( 0.0 ),
        comm_( comm ) {
    if ( comm.amRoot ) {
      CommandLineParser clp( argc, argv );
      matrixListingPath_ = clp.getValue( "listing", "no matrices given!" );
      matrixRootPath_ = clp.getValue( "rootdir", "." );
      numRuns_ = std::stoi( clp.getValue( "repeats", "1" ) );
      showPods_ = clp.checkExists( "showdist" );
      doRight_ = clp.checkExists( "right" );
      doLeft_ = clp.checkExists( "left" );
      type_ = string_to_measure_type(
          std::string( clp.getValue( "measure", "nnz" ) ) );

      paladin::print_header();
      bool parsingSuccess = true;
      if ( !parsingSuccess ) {
        if ( matrixListingPath_ == "no matrices given!" ) {
          std::cout << "ERROR! No matrix listing was provided!\n";
          parsingSuccess = false;
        }
      }
      bool badkeyfound = false;
      for ( const auto& k : clp.get_keys() ) {
        if ( k != "listing" && k != "repeats" && k != "showdist" &&
             k != "measure" && k != "rootdir" && k != "left" && k != "right" ) {
          std::cout << "\n-- POTENTIAL ERROR: key " << k
                    << " is not a recognized option!";
          badkeyfound = true;
        }
      }
      if ( badkeyfound ) {
        std::cout
            << "\n-- Allowable keys: "
            << "listing, repeats, showdist, measure, rootdir, left, right\n\n";
      }

      numMatrices_ = count_nonempty_lines( matrixListingPath_ );
      std::cout << "\nCommand line parsed!" << '\n';
      std::cout << "MPI ranks          : " << comm.numRanks << '\n';
      std::cout << "Matrix listing     : " << matrixListingPath_ << '\n';
      std::cout << "Number of matrices : " << numMatrices_ << '\n';
      std::cout << "Load measure type  : " << measure_type_description( type_ )
                << '\n';
      std::cout << "Timing repeats     : " << numRuns_ << '\n' << '\n';

      // populate and sort the matrix listing
      std::ifstream listingFile( matrixListingPath_ );
      if ( !listingFile ) {
        std::cerr << "Couldn't open the matrix listing file, "
                  << matrixListingPath_ << ", exiting!" << '\n';
        exit( 1 );
      }
      std::vector<NameMeasurePairT> matrixListing;
      while ( listingFile ) {
        std::string matrixPath;
        std::getline( listingFile, matrixPath );
        if ( matrixPath.empty() ) {
          break;
        }
        matrixPath = matrixRootPath_ + "/" + matrixPath;
        double measure = obtain_matrix_measure( matrixPath, type_ );
        matrixListing.push_back( NameMeasurePairT( matrixPath, measure ) );
        allMatrixPaths.push_back( matrixPath );
      }
      sort_name_measure_pairs( matrixListing );
      allMatrixPaths.resize( numMatrices_ );

      for ( int i = 0; i < numMatrices_; ++i ) {
        allMatrixPaths[i] = matrixListing[i].first;
      }

      // color matrices into pods
      std::vector<double> podMeasures( comm.numRanks, 0 );
      for ( const auto& f : matrixListing ) {
        int minPodIdx = std::distance(
            podMeasures.begin(),
            std::min_element( podMeasures.begin(), podMeasures.end() ) );
        podMeasures[minPodIdx] += f.second;
        podColors.push_back( minPodIdx );
      }
    }

    // distribute matrix listing size, resize storage vectors on leaf ranks
    MPI_Bcast( &numMatrices_, 1, MPI_INT, 0, MPI_COMM_WORLD );
    if ( !comm.amRoot ) {
      allMatrixPaths.resize( numMatrices_ );
      podColors.resize( numMatrices_ );
    }

    // distribute coloring vector
    MPI_Bcast( &podColors[0], numMatrices_, MPI_INT, 0, MPI_COMM_WORLD );

    // distribute the matrix listing vector, one string at a time...
    for ( auto& f : allMatrixPaths ) {
      std::string str = f;
      broadcast_string( str, comm );
      f = str;
    }

    // assign pods by color
    for ( int i = 0; i < numMatrices_; ++i ) {
      if ( podColors[i] == comm.myRank ) {
        pod.push_back( allMatrixPaths[i] );
      }
    }

    // distribute show pods condition
    MPI_Bcast( &showPods_, 1, MPI_INT, 0, MPI_COMM_WORLD );
    // distribute measure type in case we show that
    MPI_Bcast( &type_, 1, MPI_INT, 0, MPI_COMM_WORLD );

    // print matrix counts per pod and maybe the pods themselves
    for ( int rank = 0; rank < comm.numRanks; ++rank ) {
      if ( comm.myRank == rank ) {
        std::cout << "rank " << comm.myRank << ": " << pod.size() << " matrices"
                  << '\n';
        if ( showPods_ ) {
          for ( const auto& p : pod ) {
            double measure = obtain_matrix_measure( p, type_ );
            printf( "%s%0.1e%s%s\n", "    ", measure, "   ", p.c_str() );
          }
        }
      }
      comm.barrier();  // to get a contiguous list per pod
    }

    // distribute do left eigenvectors
    MPI_Bcast( &doLeft_, 1, MPI_INT, 0, MPI_COMM_WORLD );
    // distribute do right eigenvectors
    MPI_Bcast( &doRight_, 1, MPI_INT, 0, MPI_COMM_WORLD );

    // distribute number of timing repeats
    MPI_Bcast( &numRuns_, 1, MPI_INT, 0, MPI_COMM_WORLD );
  }

  void compute() {
    std::chrono::time_point<std::chrono::system_clock> start, end, startRead,
        endRead;
    localReadRunTime_ = 0.0;
    start = std::chrono::system_clock::now();

    for ( const auto& p : pod ) {
      bool firstRunOfPod = true;
      for ( int num = 0; num < numRuns_; ++num ) {
        startRead = std::chrono::system_clock::now();
        SquareMatrix mat = read_matrix_from_mm_file( p );
        endRead = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = endRead - startRead;
        localReadRunTime_ += elapsed_seconds.count();

        std::vector<double> A = mat.mat_;
        int N = mat.nrows_;
        SpectralDecomposition spectrum( N );
        std::vector<double> left, right, real( N ), imag( N ), work;
        int lwork, info;
        double wkopt;

        char leftchar = 'N';
        char rightchar = 'N';
        if ( doLeft_ ) {
          leftchar = 'V';
          left.reserve( N * N );
        }
        if ( doRight_ ) {
          rightchar = 'V';
          right.reserve( N * N );
        }

        // first query dgeev for the optimal workspace size
        lwork = -1;
        dgeev_( &leftchar, &rightchar, &N, A.data(), &N, real.data(),
                imag.data(), left.data(), &N, right.data(), &N, &wkopt, &lwork,
                &info );
        // allocate the workspace and compute the decomposition
        lwork = static_cast<std::size_t>( wkopt );
        work = std::vector<double>( lwork );
        dgeev_( &leftchar, &rightchar, &N, A.data(), &N, real.data(),
                imag.data(), left.data(), &N, right.data(), &N, work.data(),
                &lwork, &info );

        if ( firstRunOfPod ) {
          bool skipNext = false;
          for ( int i = 0; i < N - 1; ++i ) {
            spectrum.append_eigenvalue(
                std::complex<double>( real[i], imag[i] ) );
            if ( !skipNext ) {
              if ( imag[i] == -imag[i + 1] &&
                   imag[i] != 0 ) {  // if complex conjugate pair
                std::vector<std::complex<double> > leftComplexI( N ),
                    rightComplexI( N ), leftComplexIp1( N ),
                    rightComplexIp1( N );
                if ( doLeft_ ) {
                  for ( int j = i * N; j < i * N + N; ++j ) {
                    leftComplexI[j - i * N] =
                        std::complex<double>( left[j], -left[j + N] );
                    leftComplexIp1[j - i * N] =
                        std::complex<double>( left[j], left[j + N] );
                  }
                  spectrum.append_left_eigenvector( leftComplexI );
                  spectrum.append_left_eigenvector( leftComplexIp1 );
                }
                if ( doRight_ ) {
                  for ( int j = i * N; j < i * N + N; ++j ) {
                    rightComplexI[j - i * N] =
                        std::complex<double>( right[j], -right[j + N] );
                    rightComplexIp1[j - i * N] =
                        std::complex<double>( right[j], right[j + N] );
                  }
                  spectrum.append_right_eigenvector( rightComplexI );
                  spectrum.append_right_eigenvector( rightComplexIp1 );
                }
                skipNext = true;
              } else {  // if not complex conjugate pair
                std::vector<std::complex<double> > leftComplex( N ),
                    rightComplex( N );

                if ( doLeft_ ) {
                  for ( int j = i * N; j < i * N + N; ++j ) {
                    leftComplex[j - i * N] = std::complex<double>( left[j], 0 );
                  }
                  spectrum.append_left_eigenvector( leftComplex );
                }
                if ( doRight_ ) {
                  for ( int j = i * N; j < i * N + N; ++j ) {
                    rightComplex[j - i * N] =
                        std::complex<double>( right[j], 0 );
                  }
                  spectrum.append_right_eigenvector( rightComplex );
                }
                skipNext = false;
              }
            } else {
              skipNext = false;
            }
          }
          int i = N - 1;
          spectrum.append_eigenvalue(
              std::complex<double>( real[i], imag[i] ) );
          if ( !skipNext ) {
            std::vector<std::complex<double> > leftComplex( N ),
                rightComplex( N );
            if ( doLeft_ ) {
              for ( int j = i * N; j < i * N + N; ++j ) {
                leftComplex[j - i * N] = std::complex<double>( left[j], 0 );
              }
              spectrum.append_left_eigenvector( leftComplex );
            }
            if ( doRight_ ) {
              for ( int j = i * N; j < i * N + N; ++j ) {
                rightComplex[j - i * N] = std::complex<double>( right[j], 0 );
              }
              spectrum.append_right_eigenvector( rightComplex );
            }
          }
          spectrum.sort( "LM", doLeft_, doRight_ );
          spectra_.push_back( spectrum );
        }
        firstRunOfPod = false;
      }
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    localRunTime_ = elapsed_seconds.count() / static_cast<double>( numRuns_ );
    localReadRunTime_ = localReadRunTime_ / static_cast<double>( numRuns_ );
    localEigRunTime_ = localRunTime_ - localReadRunTime_;
  }

  void write_decomposition() {
    int podIdx = 0;
    for ( const auto& spectrum : spectra_ ) {
      spectrum.write_eigenvalues( std::ofstream( pod[podIdx] + ".spectrum" ) );
      if ( doLeft_ || doRight_ ) {
        spectrum.write_eigenvectors(
            std::ofstream( pod[podIdx] + ".leftvecs" ),
            std::ofstream( pod[podIdx] + ".rightvecs" ), doLeft_, doRight_ );
      }
      ++podIdx;
    }
  }

  void display_timings() {
    if ( comm_.amRoot ) {
      std::cout << '\n'
                << "Timings averaged over " << numRuns_ << " runs:" << '\n';
      std::cout << "                      runtimes (s)" << '\t' << '\t' << '\t'
                << "  fractions" << '\n';
      std::cout << "        -----------------------------------------" << '\t'
                << "-------------" << '\n';
      std::cout << "rank" << '\t' << "read" << '\t' << '\t' << "eig." << '\t'
                << '\t' << "total" << '\t' << '\t' << "read" << '\t' << "eig."
                << '\n';
    }
    comm_.barrier();

    const double readFrac = localReadRunTime_ / localRunTime_;
    const double eigFrac = localEigRunTime_ / localRunTime_;

    for ( int rank = 0; rank < comm_.numRanks; ++rank ) {
      if ( comm_.myRank == rank ) {
        printf( "%i\t%0.3e\t%0.3e\t%0.3e\t%0.2f\t%0.2f\n", comm_.myRank,
                localReadRunTime_, localEigRunTime_, localRunTime_, readFrac,
                eigFrac );
      }
      comm_.barrier();
    }

    double minRunTime, maxRunTime, totalRunTime;
    double minReadRunTime, maxReadRunTime, totalReadRunTime;
    double minEigRunTime, maxEigRunTime, totalEigRunTime;
    MPI_Reduce( &localReadRunTime_, &minReadRunTime, 1, MPI_DOUBLE, MPI_MIN, 0,
                MPI_COMM_WORLD );
    MPI_Reduce( &localReadRunTime_, &maxReadRunTime, 1, MPI_DOUBLE, MPI_MAX, 0,
                MPI_COMM_WORLD );
    MPI_Reduce( &localReadRunTime_, &totalReadRunTime, 1, MPI_DOUBLE, MPI_SUM,
                0, MPI_COMM_WORLD );
    MPI_Reduce( &localEigRunTime_, &minEigRunTime, 1, MPI_DOUBLE, MPI_MIN, 0,
                MPI_COMM_WORLD );
    MPI_Reduce( &localEigRunTime_, &maxEigRunTime, 1, MPI_DOUBLE, MPI_MAX, 0,
                MPI_COMM_WORLD );
    MPI_Reduce( &localEigRunTime_, &totalEigRunTime, 1, MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD );
    MPI_Reduce( &localRunTime_, &minRunTime, 1, MPI_DOUBLE, MPI_MIN, 0,
                MPI_COMM_WORLD );
    MPI_Reduce( &localRunTime_, &maxRunTime, 1, MPI_DOUBLE, MPI_MAX, 0,
                MPI_COMM_WORLD );
    MPI_Reduce( &localRunTime_, &totalRunTime, 1, MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD );

    if ( comm_.amRoot ) {
      const double numRanks = static_cast<double>( comm_.numRanks );
      const double avgReadRunTime = totalReadRunTime / numRanks;
      const double avgEigRunTime = totalEigRunTime / numRanks;
      const double avgRunTime = totalRunTime / numRanks;
      const double avgReadFrac = avgReadRunTime / avgRunTime;
      const double avgEigFrac = avgEigRunTime / avgRunTime;

      std::cout << "-----------------------------------------------------------"
                   "----------"
                << '\n';
      printf( "%s\t%0.3e\t%0.3e\t%0.3e\n", "min", minReadRunTime, minEigRunTime,
              minRunTime );
      printf( "%s\t%0.3e\t%0.3e\t%0.3e\t%0.2f\t%0.2f\n", "avg", avgReadRunTime,
              avgEigRunTime, avgRunTime, avgReadFrac, avgEigFrac );
      printf( "%s\t%0.3e\t%0.3e\t%0.3e\n", "max", maxReadRunTime, maxEigRunTime,
              maxRunTime );

      const double li = maxRunTime / avgRunTime - 1.0;
      std::cout << '\n' << "load imbalance = " << 100 * li << "%" << '\n';
    }
  }
};
}

#endif /* SRC_PALADIN_HPP_ */
