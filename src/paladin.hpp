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
#include <types.hpp>
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
  std::vector<SpectrumT> spectra_;
  double localRunTime_, localReadRunTime_, localEigRunTime_;
  std::string matrixListingPath_, matrixRootPath_;
  int numRuns_, numMatrices_;
  bool showPods_;
  LoadPredictor type_;
  const MpiComm& comm_;

 public:
  Paladin( int& argc, char* argv[], const MpiComm& comm )
      : localRunTime_( 0.0 ),
        localReadRunTime_( 0.0 ),
        localEigRunTime_( 0.0 ),
        comm_( comm ) {
    if ( comm.amRoot ) {
      CommandLineParser clp( argc, argv );
      matrixListingPath_ = clp.getValue( "matrices", "no matrices given!" );
      matrixRootPath_ = clp.getValue( "rootdir", "." );
      numRuns_ = std::stoi( clp.getValue( "repeats", "1" ) );
      showPods_ = clp.checkExists( "showdist" );
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
        if ( k != "matrices" && k != "repeats" && k != "showdist" &&
             k != "measure" && k != "rootdir" ) {
          std::cout << "\n-- POTENTIAL ERROR: key " << k
                    << " is not a recognized option!";
          badkeyfound = true;
        }
      }
      if ( badkeyfound ) {
        std::cout << "\n-- Allowable keys: "
                  << "matrices, repeats, showdist, measure, rootdir\n\n";
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
        MeasureT measure = obtain_matrix_measure( matrixPath, type_ );
        matrixListing.push_back( NameMeasurePairT( matrixPath, measure ) );
        allMatrixPaths.push_back( matrixPath );
      }
      sort_name_measure_pairs( matrixListing );
      allMatrixPaths.resize( numMatrices_ );

      for ( int i = 0; i < numMatrices_; ++i ) {
        allMatrixPaths[i] = matrixListing[i].first;
      }

      // color matrices into pods
      std::vector<MeasureT> podMeasures( comm.numRanks, 0 );
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
            MeasureT measure = obtain_matrix_measure( p, type_ );
            printf( "%s%0.1e%s%s\n", "    ", measure, "   ", p.c_str() );
          }
        }
      }
      comm.barrier();  // to get a contiguous list per pod
    }

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

        SpectrumT s;
        DVecT A = mat.mat_;
        int N = mat.nrows_;
        DVecT left, right, real( N ), imag( N ), work;
        int lwork, info;
        double wkopt;
        const char nchar = 'N';

        // first query dgeev for the optimal workspace size
        lwork = -1;
        dgeev_( &nchar, &nchar, &N, A.data(), &N, real.data(), imag.data(),
                left.data(), &N, right.data(), &N, &wkopt, &lwork, &info );
        // allocate the workspace and compute the eigenvalues
        lwork = (std::size_t)wkopt;
        work = DVecT( lwork );
        dgeev_( &nchar, &nchar, &N, A.data(), &N, real.data(), imag.data(),
                left.data(), &N, right.data(), &N, work.data(), &lwork, &info );

        if ( firstRunOfPod ) {
          for ( int i = 0; i < N; ++i ) {
            s.push_back( EigenvalueT( real[i], imag[i] ) );
          }
          sort_spectrum( s, "LM" );
          spectra_.push_back( s );
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

  void write_eigenvalues() {
    int podIdx = 0;
    for ( const auto& spectrum : spectra_ ) {
      write_spectrum( spectrum, std::ofstream( pod[podIdx] + ".spectrum" ) );
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
