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
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "types.hpp"
#include "mpi-util.hpp"
#include "spectrum-util.hpp"
#include "square-matrix.hpp"
#include "lapack-wrapper.hpp"
#include "command-line-parser.hpp"


namespace paladin{};
using namespace paladin;



struct MpiComm
{
  int numRanks, myRank, rootRank;
  bool amRoot;

  MpiComm( int& argc, char *argv[] )
  {
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &numRanks );
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    rootRank = 0;
    amRoot = ( myRank == rootRank );
  }

  ~MpiComm()
  {
    MPI_Finalize();
  }


  void barrier() const { MPI_Barrier( MPI_COMM_WORLD );}
};



/**
   * @class Paladin
   * @brief PArallel LApack DIstributor for many eigenvalue calculations
   */
class Paladin
{
protected:
  std::vector<std::string> allFlockPaths, pod;
  std::vector<int> podColors;
  std::vector<SpectrumT> spectra_;
  std::vector<std::string> spectraOutputFileNames_;
  double localRunTime_, localReadRunTime_, localEigRunTime_;
  std::string flockPath_;
  int numRuns_, flockSize_;
  bool showPods_;
  LoadPredictor type_;



  void print_header()
  {
    std::cout << '\n' << '\n';
    std::cout << " ____   ____  _       ____  ___    ____  ____  " << '\n'
              << "|    \\ /    || |     /    ||   \\  |    ||    \\ " << '\n'
              << "|  o  |  o  || |    |  o  ||    \\  |  | |  _  |" << '\n'
              << "|   _/|     || |___ |     ||  D  | |  | |  |  |" << '\n'
              << "|  |  |  _  ||     ||  _  ||     | |  | |  |  |" << '\n'
              << "|  |  |  |  ||     ||  |  ||     | |  | |  |  |" << '\n'
              << "|__|  |__|__||_____||__|__||_____||____||__|__|" << '\n' << '\n'
              << "PArallel LApack DIstributor for many eigenvalue calculations" << '\n'
              << "Copyright (c) 2016 Mike Hansen" << '\n'
              << "------------------------------------------------" << '\n';
    std::cout << '\n' << '\n';
  }


public:

  Paladin( int& argc, char *argv[], const MpiComm& comm )
{
    if( comm.amRoot ){
      CommandLineParser clp( argc, argv );
      flockPath_ = clp.getValue( "-flock", "no-flock-provided!" );
      std::string measureString = clp.getValue( "-load-measure", "nnz" );
      numRuns_ = std::stoi( clp.getValue( "-timing-repeats", "1" ) );
      showPods_ = clp.checkExists( "-show-pods" );
      type_ = string_to_measure_type( measureString );
      flockSize_ = count_nonempty_lines( flockPath_ );

      print_header();
      std::cout << "MPI ranks     : " << comm.numRanks << '\n';
      std::cout << "Flock listing : " << flockPath_ << '\n';
      std::cout << "Matrices      : " << flockSize_ << '\n';
      std::cout << "Measure type  : " << measure_type_description( type_ ) << '\n';
      std::cout << "Timing runs   : " << numRuns_ << '\n' << '\n';

      // populate and sort the flock
      std::ifstream flockFile( flockPath_ );
      std::vector<NameMeasurePairT> flock;
      while( flockFile ){
        std::string matrixPath;
        std::getline( flockFile, matrixPath );
        if( matrixPath.empty() ) break;
        MeasureT measure = obtain_matrix_measure( matrixPath, type_ );
        flock.push_back( NameMeasurePairT( matrixPath, measure ) );
        allFlockPaths.push_back( matrixPath );
      }
      sort_name_measure_pairs( flock );
      allFlockPaths.resize( flockSize_ );
      for( int i=0; i<flockSize_; ++i ){
        allFlockPaths[i] = flock[i].first;
      }

      // color matrices into pods
      std::vector<MeasureT> podMeasures( comm.numRanks, 0 );
      for( size_t i=0; i<flock.size(); ++i ){
        int minPodIdx = std::distance( podMeasures.begin(), std::min_element( podMeasures.begin(), podMeasures.end() ) );
        podMeasures[minPodIdx] += flock[i].second;
        podColors.push_back( minPodIdx );
      }
    }

    // distribute flock size, resize storage vectors on leaf ranks
    MPI_Bcast( &flockSize_, 1, MPI_INT, 0, MPI_COMM_WORLD );
    if( !comm.amRoot ){
      allFlockPaths.resize( flockSize_ );
      podColors.resize( flockSize_ );
    }

    // distribute coloring vector
    MPI_Bcast( &podColors[0], flockSize_, MPI_INT, 0, MPI_COMM_WORLD );

    // distribute the flock vector, one string at a time...
    for( int i=0; i<flockSize_; ++i ){
      std::string str = allFlockPaths[i];
      broadcast_string( str, comm.myRank, comm.numRanks );
      allFlockPaths[i] = str;
    }

    // assign pods by color
    for( int i=0; i<flockSize_; ++i ){
      if( podColors[i] == comm.myRank ){
        pod.push_back( allFlockPaths[i] );
      }
    }

    // distribute show pods condition
    MPI_Bcast( &showPods_, 1, MPI_INT, 0, MPI_COMM_WORLD );
    // distribute measure type in case we show that
    MPI_Bcast( &type_, 1, MPI_INT, 0, MPI_COMM_WORLD );

    // print matrix counts per pod and maybe the pods themselves
    for( int rank=0; rank<comm.numRanks; ++rank ){
      if( comm.myRank == rank ){
        std::cout << "rank " << comm.myRank << ": " << pod.size() << " matrices" << '\n';
        if( showPods_ ){
          for( size_t i=0; i<pod.size(); ++i ){
            MeasureT measure = obtain_matrix_measure( pod[i], type_ );
            printf( "%s%0.1e%s%s\n", "    ", measure, "   ", pod[i].c_str() );
          }
        }
      }
      comm.barrier(); // to get a contiguous list per pod
    }

    // distribute number of timing repeats
    MPI_Bcast( &numRuns_, 1, MPI_INT, 0, MPI_COMM_WORLD );
}



  void compute()
  {
    std::chrono::time_point<std::chrono::system_clock> start, end, startRead, endRead;
    localReadRunTime_ = 0.0;
    start = std::chrono::system_clock::now();

    for( size_t idx=0; idx<pod.size(); ++idx ){
      bool firstRunOfPod = true;
      for( int num=0; num<numRuns_; ++num ){
        startRead = std::chrono::system_clock::now();
        SquareMatrix mat = read_matrix_from_mm_file( pod[idx] );
        endRead = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = endRead - startRead;
        localReadRunTime_ += elapsed_seconds.count();
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
        std::string spectrumOutputFileName = "";
        if( firstRunOfPod ){
          spectra_.push_back( s );
          spectraOutputFileNames_.push_back( spectrumOutputFileName );
        }
        firstRunOfPod = false;

        //        std::ostringstream strsMyRank;
        //        strsMyRank << comm.myRank;
        //        std::string myRankStr = strsMyRank.str();
        //        std::ostringstream strsIdx;
        //        strsIdx << idx;
        //        std::string idxStr = strsIdx.str();
        //      std::string output( "eigs_" + comm.myRankStr + "_" + idxStr + ".eig" );
        //      std::ofstream ofs;
        //      ofs.open( output, std::ofstream::out | std::ofstream::app );
        //      show_spectrum( s );
      }
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    localRunTime_ = elapsed_seconds.count() / static_cast<double>( numRuns_ );
    localReadRunTime_ = localReadRunTime_ / static_cast<double>( numRuns_ );
    localEigRunTime_ = localRunTime_ - localReadRunTime_;
  }

  void write_eigenvalues() {} // TODO: write me!

  void display_timings( const MpiComm& comm )
  {
    if( comm.amRoot ){
      std::cout << '\n' << "Timings averaged over " << numRuns_ << " runs:" << '\n';
      std::cout << "                      runtimes (s)"   << '\t' << '\t' << '\t' << "  fractions" << '\n';
      std::cout <<         "        -----------------------------------------" << '\t' << "-------------" << '\n';
      std::cout << "rank" << '\t' << "read" << '\t' << '\t' << "eig." << '\t' << '\t' << "total" << '\t' << '\t' << "read" << '\t' << "eig." << '\n';
    }
    comm.barrier();

    const double readFrac = localReadRunTime_ / localRunTime_;
    const double eigFrac  = localEigRunTime_  / localRunTime_;

    for( int rank=0; rank<comm.numRanks; ++rank ){
      if( comm.myRank == rank ){
        printf( "%i\t%0.3e\t%0.3e\t%0.3e\t%0.2f\t%0.2f\n", comm.myRank, localReadRunTime_, localEigRunTime_, localRunTime_, readFrac, eigFrac );
      }
      comm.barrier();
    }

    double minRunTime = 0.0;
    double maxRunTime = 0.0;
    double totalRunTime = 0.0;
    double minReadRunTime = 0.0;
    double maxReadRunTime = 0.0;
    double totalReadRunTime = 0.0;
    double minEigRunTime = 0.0;
    double maxEigRunTime = 0.0;
    double totalEigRunTime = 0.0;
    MPI_Reduce( &localReadRunTime_, &minReadRunTime  , 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &localReadRunTime_, &maxReadRunTime  , 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &localReadRunTime_, &totalReadRunTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &localEigRunTime_, &minEigRunTime  , 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &localEigRunTime_, &maxEigRunTime  , 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &localEigRunTime_, &totalEigRunTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &localRunTime_, &minRunTime  , 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &localRunTime_, &maxRunTime  , 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &localRunTime_, &totalRunTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    if( comm.amRoot ){
      const double avgReadRunTime = totalReadRunTime / static_cast<double>( comm.numRanks );
      const double avgEigRunTime  = totalEigRunTime  / static_cast<double>( comm.numRanks );
      const double avgRunTime     = totalRunTime     / static_cast<double>( comm.numRanks );
      const double avgReadFrac    = avgReadRunTime / avgRunTime;
      const double avgEigFrac     = avgEigRunTime  / avgRunTime;

      std::cout << "---------------------------------------------------------------------" << '\n';
      printf( "%s\t%0.3e\t%0.3e\t%0.3e\n", "min", minReadRunTime, minEigRunTime, minRunTime );
      printf( "%s\t%0.3e\t%0.3e\t%0.3e\t%0.2f\t%0.2f\n", "avg", avgReadRunTime, avgEigRunTime, avgRunTime, avgReadFrac, avgEigFrac );
      printf( "%s\t%0.3e\t%0.3e\t%0.3e\n", "max", maxReadRunTime, maxEigRunTime, maxRunTime );

      const double lifComp = maxRunTime / avgRunTime - 1.0;
      std::cout << '\n' << "load imbalance = " << 100 * lifComp << "%" << '\n';
    }
  }

};





int main( int argc, char *argv[] )
{

  MpiComm comm( argc, argv );

  Paladin roland( argc, argv, comm );
  roland.compute();
  roland.display_timings( comm );

  return 0;
}
