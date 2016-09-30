/*
 * The MIT License
 *
 * Copyright (c) 2016 Mike Hansen
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


#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <complex>

#include <mpi.h>

#include "types.hpp"
#include "mpi-util.hpp"
#include "spectrum-util.hpp"
#include "square-matrix.hpp"
#include "lapack-wrapper.hpp"
#include "command-line-parser.hpp"


namespace paladin{};
using namespace paladin;


int main( int argc, char *argv[] )
{


  int numRanks, myRank;

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &numRanks );
  MPI_Comm_rank( MPI_COMM_WORLD, &myRank );

  const int rootRank = 0;


  std::vector<std::string> allFlockPaths, pod;
  std::vector<int> podColors;
  int flockSize;
  bool showPods = false;
  int numRuns = 1;

  if( myRank == rootRank ){

    // parse command line options
    CommandLineParser clp( argc, argv );

    std::string flockPath = clp.getValue( "-flock", "no-flock-provided!" );
    std::string measureString = clp.getValue( "-load-measure", "nnz" );
    numRuns = std::stoi( clp.getValue( "-timing-repeats", "1" ) );
    showPods = clp.checkExists( "-show-pods" );

    LoadPredictor type = string_to_measure_type( measureString );

    int numMatrices = count_nonempty_lines( flockPath );

    std::cout << '\n' << "--------------------------------" << '\n' << '\n';
    std::cout << "MPI ranks     : " << numRanks << '\n' << '\n';
    std::cout << "Flock listing : " << flockPath << '\n' << '\n';
    std::cout << "Matrices      : " << numMatrices << '\n' << '\n';
    std::cout << "Measure type  : " << measure_type_description( type ) << '\n' << '\n';
    std::cout << "Timing runs   : " << numRuns << '\n' << '\n';
    std::cout << "--------------------------------" << '\n' << '\n';

    // populate and sort the flock
    std::ifstream flockFile( flockPath );
    std::vector<NameMeasurePairT> flock;
    while( flockFile ){
      std::string matrixPath;
      std::getline( flockFile, matrixPath );
      if( matrixPath.empty() ) break;
      MeasureT measure = obtain_matrix_measure( matrixPath, type );
      flock.push_back( NameMeasurePairT( matrixPath, measure ) );
      allFlockPaths.push_back( matrixPath );
    }
    sort_name_measure_pairs( flock );
    flockSize = flock.size();
    allFlockPaths.resize( flockSize );
    for( int i=0; i<flockSize; ++i )
      allFlockPaths[i] = flock[i].first;

    // color matrices into pods
    std::vector<MeasureT> podMeasures( numRanks, 0 );
    for( size_t i=0; i<flock.size(); ++i ){
      int minPodIdx = std::distance( podMeasures.begin(), std::min_element( podMeasures.begin(), podMeasures.end() ) );
      podMeasures[minPodIdx] += flock[i].second;
      podColors.push_back( minPodIdx );
    }

    for( int rank=0; rank<numRanks; ++rank )
      std::cout << "rank " << rank << " load measure = " << podMeasures[rank] << '\n';

    // print predicted load imbalance
    const double minPodIdx = static_cast<double>( podMeasures[ std::distance( podMeasures.begin(),
                                                                              std::min_element( podMeasures.begin(),
                                                                                                podMeasures.end() ) ) ] );
    const double maxPodIdx = static_cast<double>( podMeasures[ std::distance( podMeasures.begin(),
                                                                              std::max_element( podMeasures.begin(),
                                                                                                podMeasures.end() ) ) ] );
    const double lifPred = ( maxPodIdx - minPodIdx ) / minPodIdx;
    std::cout << "predicted load imbalance = " << 100 * lifPred << "%" << '\n' << '\n';
  }

  // distribute flock size, resize storage vectors on leaf ranks
  MPI_Bcast( &flockSize, 1, MPI_INT, 0, MPI_COMM_WORLD );
  if( myRank != rootRank ){
    allFlockPaths.resize( flockSize );
    podColors.resize( flockSize );
  }

  // distribute coloring vector
  MPI_Bcast( &podColors[0], flockSize, MPI_INT, 0, MPI_COMM_WORLD );

  // distribute the flock vector, one string at a time...
  for( int i=0; i<flockSize; ++i ){
    std::string str = allFlockPaths[i];
    broadcast_string( str, myRank, numRanks );
    allFlockPaths[i] = str;
  }

  // assign pods by color
  for( int i=0; i<flockSize; ++i )
    if( podColors[i] == myRank )
      pod.push_back( allFlockPaths[i] );



  // distribute show pods condition
  MPI_Bcast( &showPods, 1, MPI_INT, 0, MPI_COMM_WORLD );
  // print matrix counts per pod and maybe the pods themselves
  for( int rank=0; rank<numRanks; ++rank ){
    if( myRank == rank ){
      std::cout << "rank " << myRank << ": " << pod.size() << " matrices" << '\n';
      if( showPods )
        for( size_t i=0; i<pod.size(); ++i )
          std::cout << "    " << pod[i] << '\n';
    }
    MPI_Barrier( MPI_COMM_WORLD ); // to get a contiguous list per pod
  }

  // distribute number of timing repeats
  MPI_Bcast( &numRuns, 1, MPI_INT, 0, MPI_COMM_WORLD );

  // run simultaneous eigenvalue calculations
  std::chrono::time_point<std::chrono::system_clock> start, end, startRead, endRead;
  double localReadRunTime = 0.0;
  start = std::chrono::system_clock::now();

  for( int num=0; num<numRuns; ++num ){
    for( size_t idx=0; idx<pod.size(); ++idx ){
      startRead = std::chrono::system_clock::now();
      SquareMatrix mat = read_matrix_from_mm_file( pod[idx] );
      endRead = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = endRead - startRead;
      localReadRunTime += elapsed_seconds.count();
      SpectrumT s;
      DVecT A = mat.mat_;
      int N = mat.nrows_;
      DVecT left;
      DVecT right;
      DVecT real( N );
      DVecT imag( N );
      int lwork;
      DVecT work;
      double wkopt;
      int info;
      const char nchar = 'N';
      lwork = -1;
      dgeev_( &nchar, &nchar, &N, A.data(), &N, real.data(), imag.data(), left.data(), &N, right.data(), &N, &wkopt, &lwork, &info );
      lwork = (size_t) wkopt;
      work = DVecT( lwork );
      dgeev_( &nchar, &nchar, &N, A.data(), &N, real.data(), imag.data(), left.data(), &N, right.data(), &N, work.data(), &lwork, &info );
      for( int i=0; i<N; ++i )
        s.push_back( EigenvalueT( real[i], imag[i] ) );
      sort_spectrum( s, "LM" );

      std::ostringstream strsMyRank;
      strsMyRank << myRank;
      std::string myRankStr = strsMyRank.str();

      std::ostringstream strsIdx;
      strsIdx << idx;
      std::string idxStr = strsIdx.str();

      //      std::string output( "eigs_" + myRankStr + "_" + idxStr + ".eig" );
      //      std::ofstream ofs;
      //      ofs.open( output, std::ofstream::out | std::ofstream::app );
      //      show_spectrum( s );
    }
  }

  end = std::chrono::system_clock::now();


  // display timings
  std::chrono::duration<double> elapsed_seconds = end - start;
  double localRunTime = elapsed_seconds.count() / static_cast<double>( numRuns );
  localReadRunTime = localReadRunTime / static_cast<double>( numRuns );

  std::cout << "rank " << myRank << ": total cputime = " << localRunTime << " s averaged over " << numRuns << " runs" << '\n'
      << "cputime on matrix input: " << localReadRunTime << " s (" << 100.0 * localReadRunTime / localRunTime << "%)" << '\n' << '\n';

  double minRunTime = 0.0;
  double maxRunTime = 0.0;
  double minReadRunTime = 0.0;
  double maxReadRunTime = 0.0;
  MPI_Reduce( &localReadRunTime, &minReadRunTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
  MPI_Reduce( &localReadRunTime, &maxReadRunTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
  MPI_Reduce( &localRunTime, &minRunTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
  MPI_Reduce( &localRunTime, &maxRunTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

  if( myRank == rootRank ){
    std::cout << "min run  time over " << numRanks << " ranks = " << minRunTime     << '\n'
        << "max run  time over " << numRanks << " ranks = " << maxRunTime     << '\n'
        << "min read time over " << numRanks << " ranks = " << minReadRunTime << '\n'
        << "max read time over " << numRanks << " ranks = " << maxReadRunTime << '\n';
    const double lifComp = ( maxRunTime - minRunTime ) / minRunTime;
    std::cout << '\n' << "measured load imbalance = " << 100 * lifComp << "%" << '\n';
  }




  MPI_Finalize();
  return 0;
}
