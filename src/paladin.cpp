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
#include <mpi.h>
#include <algorithm>
#include <complex>


// types
typedef std::vector<double> DVecT;
typedef std::vector<std::string> StrVecT;
typedef double MeasureT;
typedef std::pair<std::string, MeasureT> NameMeasurePairT;
typedef std::complex<double> EigT;
typedef std::vector<EigT> SpectrumT;
// --------------------------------


// to send a string from root to all other ranks
void broadcast_string( std::string& str, const int myRank, const int numRanks, const int rootRank = 0 )
{
  if( myRank == rootRank ){
    for( int sendToRank=0; sendToRank<numRanks; ++sendToRank )
      if( sendToRank != rootRank )
        MPI_Send( str.c_str(), str.size(), MPI_CHAR, sendToRank, 0, MPI_COMM_WORLD );
  }
  else{
    MPI_Status status;
    int size;
    MPI_Probe( rootRank, 0, MPI_COMM_WORLD, &status );
    MPI_Get_count( &status, MPI_CHAR, &size );

    std::vector<char> tmp( size );
    MPI_Recv( tmp.data(), size, MPI_CHAR, rootRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
    str.assign( tmp.begin(), tmp.end() );
  }
}
// --------------------------------

// measure type enumeration, string conversions
enum class MeasureType
{
  DIMENSION,
  NUMNONZEROS,
  DIMCUBED
};

MeasureType string_to_measure_type( const std::string& str )
{
  MeasureType type;
  if     ( str == "dim" ) type = MeasureType::DIMENSION;
  else if( str == "nnz" ) type = MeasureType::NUMNONZEROS;
  else if( str == "dcb" ) type = MeasureType::DIMCUBED;
  else                    type = MeasureType::DIMCUBED;
  return type;
}

std::string measure_type_description( const MeasureType& type )
{
  std::string str;
  if     ( type == MeasureType::DIMENSION )   str = "matrix dimension";
  else if( type == MeasureType::NUMNONZEROS ) str = "number of nonzeros";
  else if( type == MeasureType::DIMCUBED )    str = "matrix dimension cubed";
  else                                        str = "matrix dimension";
  return str;
}
// --------------------------------

// for measure-name sorting
struct MeasureComparator
{
  inline bool operator()( const NameMeasurePairT& left, const NameMeasurePairT& right )
  {
    return ( left.second > right.second );
  }
};
// --------------------------------

// for sorting and writing spectra
struct EigenSorter{
  virtual ~EigenSorter() {}
  virtual bool operator() ( const EigT& L, const EigT& R ) const { return ( std::norm( L ) > std::norm( R ) );}
};

struct LReal : public EigenSorter { bool operator() ( const EigT& L, const EigT& R ) const { return ( L.real() > R.real() );} };
struct SReal : public EigenSorter { bool operator() ( const EigT& L, const EigT& R ) const { return ( L.real() < R.real() );} };
struct LImag : public EigenSorter { bool operator() ( const EigT& L, const EigT& R ) const { return ( L.imag() > R.imag() );} };
struct SImag : public EigenSorter { bool operator() ( const EigT& L, const EigT& R ) const { return ( L.imag() < R.imag() );} };
struct LMag  : public EigenSorter { bool operator() ( const EigT& L, const EigT& R ) const { return ( std::norm( L ) > std::norm( R ) );} };
struct SMag  : public EigenSorter { bool operator() ( const EigT& L, const EigT& R ) const { return ( std::norm( L ) < std::norm( R ) );} };


void sort_spectrum( SpectrumT& s, const std::string sortStr )
{
  EigenSorter* sorter;
  if     ( sortStr == "SM" ) sorter = new SMag;
  else if( sortStr == "LM" ) sorter = new LMag;
  else if( sortStr == "SR" ) sorter = new SReal;
  else if( sortStr == "LR" ) sorter = new LReal;
  else if( sortStr == "SI" ) sorter = new SImag;
  else if( sortStr == "LI" ) sorter = new LImag;
  else                       sorter = new LMag;
  std::sort( s.begin(), s.end(), std::ref( *sorter ) );
}

void show_spectrum( const SpectrumT& s )
{
  for( size_t i=0; i<s.size(); ++i )
    std::cout << "lambda_" << i << " = " << s[i].real() << " + " << s[i].imag() << "j" << std::endl;
}

void write_spectrum( const SpectrumT& s, std::ofstream& realout, std::ofstream& imagout )
{
  for( size_t i=0; i<s.size(); ++i ){
    realout << s[i].real() << std::endl;
    imagout << s[i].imag() << std::endl;
  }
}
// --------------------------------


// container for a square matrix, wrapping a DVecT and supporting 1-based indexing
// 1-based indexing because matrix market input uses 1-based indexing
struct SquareMatrix
{
  const int nrows_, nelem_;
  DVecT mat_;

  SquareMatrix( const int nrows ) : nrows_( nrows ), nelem_( nrows*nrows ), mat_( nelem_, 0.0 ) {}

  // one-based indexing
  int idx1D( const int i, const int j ) const { return (i-1)*nrows_+(j-1);}
  double& operator()( const int idx )      { return mat_[idx-1];}
  double& operator()( const int i, const int j ) { return mat_[idx1D(i,j)];}

  // setting based on a 'triple'
  void set_element( const int i, const int j, const double aij ) { (*this)( i, j ) = aij;}

  void print() const
  {
    std::cout << "[";
    for( int i=0; i<nrows_; ++i ){
      for( int j=0; j<nrows_; ++j )
        std::cout << " " << mat_[i*nrows_+j] << ",";
      if( i==nrows_-1 ) std::cout << " ]" << std::endl;
      else              std::cout << std::endl;
    }
  }
};
// --------------------------------


// splits a string by whitespace
StrVecT split( const std::string& str )
{
  StrVecT split;
  std::stringstream stream( str );
  std::string buffer;
  while( stream >> buffer )
    split.push_back( buffer );
  return split;
}
// --------------------------------


// read a matrix market file into a dense SquareMatrix, compute its measure, etc.
// matrix market header must only be one line long!
void read_header( std::ifstream& file )
{
  std::string head;
  std::getline( file, head );
  StrVecT Header = split( head );
  assert( Header[0] == "MatrixMarket" || Header[0] == "%%MatrixMarket" );
  assert( Header[1] == "matrix"       );
  assert( Header[2] == "coordinate"   );
  assert( Header[3] == "real"         );
  assert( Header[4] == "general"      );
}

SquareMatrix allocate_matrix( std::ifstream& file )
{
  std::string line2;
  std::getline( file, line2 );
  StrVecT Line2 = split( line2 );

  int nrows = std::stoi( Line2[0] );
  int ncols = std::stoi( Line2[1] );
  int nnz   = std::stoi( Line2[2] );

  assert( nrows == ncols );
  assert( nnz <= nrows * ncols );

  return SquareMatrix( nrows );
}

void read_matrix( SquareMatrix& A, std::ifstream& file )
{
  while( file ){
    std::string contents;
    std::getline( file, contents );
    if( contents.empty() ) break;
    StrVecT stuff = split( contents );
    const int i = std::stoi( stuff[0] );
    const int j = std::stoi( stuff[1] );
    const double      v = std::stod( stuff[2] );
    A( i, j ) = v;
  }
}

int count_nonempty_lines( std::ifstream file )
{
  int nlines = 0;
  std::string line;
  while ( std::getline( file, line ) && line != "" ) ++nlines;
  return nlines;
}

int count_nonempty_lines( const std::string& filename )
{
  return count_nonempty_lines( std::ifstream( filename.c_str() ) );
}

int read_matrix_dimension_from_mm_file( const std::string& matrixPath )
{
  std::ifstream file( matrixPath );
  if( !file ){std::cerr << "Couldn't open the file, exiting!" << std::endl; exit(1);}
  std::string line;
  std::getline( file, line );
  std::getline( file, line );

  StrVecT Line = split( line );
  int nrows = std::stoi( Line[0] );
  int ncols = std::stoi( Line[1] );

  assert( nrows == ncols );
  return nrows;
}

int read_matrix_nnzeros_from_mm_file( const std::string& matrixPath )
{
  std::ifstream file( matrixPath );
  if( !file ){std::cerr << "Couldn't open the file, exiting!" << std::endl; exit(1);}
  std::string line;
  std::getline( file, line );
  std::getline( file, line );

  StrVecT Line = split( line );
  int nrows = std::stoi( Line[0] );
  int ncols = std::stoi( Line[1] );
  int nnz   = std::stoi( Line[2] );

  assert( nrows == ncols );
  return nnz;
}

MeasureT obtain_matrix_measure( const std::string& matrixPath,
                                const MeasureType& type = MeasureType::DIMENSION )
{
  MeasureT measure = 0;
  MeasureT dim = ( MeasureT ) read_matrix_dimension_from_mm_file( matrixPath );
  MeasureT nnz = ( MeasureT ) read_matrix_nnzeros_from_mm_file( matrixPath );
  switch( type ){
    case MeasureType::DIMENSION:   measure = dim;             break;
    case MeasureType::NUMNONZEROS: measure = nnz;             break;
    case MeasureType::DIMCUBED:    measure = dim * dim * dim; break;
    default:
      std::cerr << "Invalid measure type given, exiting!" << std::endl;
      exit(1);
  }
  return measure;
}

SquareMatrix read_matrix_from_mm_file( const std::string& matrixPath )
{
  std::ifstream file( matrixPath );
  if( !file ){std::cerr << "Couldn't open the file, exiting!" << std::endl; exit(1);}
  read_header( file );
  SquareMatrix A = allocate_matrix( file );
  read_matrix( A, file );
  return A;
}
// --------------------------------



/*
 * dgeev_ - full spectral decomposition of a general, double-precision matrix
 */
extern "C" void dgeev_(
    const char *doleft,  /* 'N' for no left eigenmatrix and 'V' for computing it, char reference */
    const char *doright, /* 'N' for no right eigenmatrix and 'V' for computing it, char reference */
    int *N,              /* number of rows in A, integer reference */
    double *A,           /* A matrix in column-major form, double array */
    int *lda,            /* leading dimension of A, integer reference */
    double *real,        /* eigenvalue real parts, double array */
    double *imag,        /* eigenvalue imaginary parts, double array */
    double *leftmat,     /* left eigenmatrix, double array */
    int *ldleft,         /* leading dim. of left eig matrix, integer reference */
    double *rightmat,    /* right eigenmatrix, double array */
    int *ldright,        /* leading dim. of right eig matrix, integer reference */
    double *work,        /* work array size, double reference */
    int *lwork,          /* work array, double array */
    int *info            /* info structure, integer reference */ );
// --------------------------------



// parsing command line options
class CommandLineParser{
protected:
  StrVecT pieces;
public:
  CommandLineParser( const int argc, char *argv[] )
{
    for( int i=1; i<argc; ++i )
      this->pieces.push_back( std::string( argv[i] ) );
}


  bool checkExists( const std::string& name ) const
  {
    return std::find( pieces.begin(), pieces.end(), name ) != pieces.end();
  }

  std::string getValue( const std::string& name, const std::string& defaultStr ) const
  {
    if( this->checkExists( name ) ){
      std::vector<std::string>::const_iterator itr;
      itr = std::find( pieces.begin(), pieces.end(), name );
      if (itr != pieces.end() && ++itr != pieces.end()){
        return *itr;
      }
      else return defaultStr;
    }
    else return defaultStr;
  }
};
// --------------------------------










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

    MeasureType type = string_to_measure_type( measureString );

    int numMatrices = count_nonempty_lines( flockPath );

    std::cout << std::endl << "--------------------------------" << std::endl << std::endl;
    std::cout << "MPI ranks     : " << numRanks << std::endl << std::endl;
    std::cout << "Flock listing : " << flockPath << std::endl << std::endl;
    std::cout << "Matrices      : " << numMatrices << std::endl << std::endl;
    std::cout << "Measure type  : " << measure_type_description( type ) << std::endl << std::endl;
    std::cout << "Timing runs   : " << numRuns << std::endl << std::endl;
    std::cout << "--------------------------------" << std::endl << std::endl;

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
    std::sort( flock.begin(), flock.end(), MeasureComparator() );
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
      std::cout << "rank " << rank << " load measure = " << podMeasures[rank] << std::endl;

    // print predicted load imbalance
    const double minPodIdx = static_cast<double>( podMeasures[ std::distance( podMeasures.begin(),
                                                                              std::min_element( podMeasures.begin(),
                                                                                                podMeasures.end() ) ) ] );
    const double maxPodIdx = static_cast<double>( podMeasures[ std::distance( podMeasures.begin(),
                                                                              std::max_element( podMeasures.begin(),
                                                                                                podMeasures.end() ) ) ] );
    const double lifPred = ( maxPodIdx - minPodIdx ) / minPodIdx;
    std::cout << "predicted load imbalance = " << 100 * lifPred << "%" << std::endl << std::endl;
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
      std::cout << "rank " << myRank << ": " << pod.size() << " matrices" << std::endl;
      if( showPods )
        for( size_t i=0; i<pod.size(); ++i )
          std::cout << "    " << pod[i] << std::endl;
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
        s.push_back( EigT( real[i], imag[i] ) );
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

  std::cout << "rank " << myRank << ": total cputime = " << localRunTime << " s averaged over " << numRuns << " runs" << std::endl
            << "cputime on matrix input: " << localReadRunTime << " s (" << 100.0 * localReadRunTime / localRunTime << "%)" << std::endl << std::endl;

  double minRunTime = 0.0;
  double maxRunTime = 0.0;
  double minReadRunTime = 0.0;
  double maxReadRunTime = 0.0;
  MPI_Reduce( &localReadRunTime, &minReadRunTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
  MPI_Reduce( &localReadRunTime, &maxReadRunTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
  MPI_Reduce( &localRunTime, &minRunTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
  MPI_Reduce( &localRunTime, &maxRunTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

  if( myRank == rootRank ){
    std::cout << "min run  time over " << numRanks << " ranks = " << minRunTime     << std::endl
              << "max run  time over " << numRanks << " ranks = " << maxRunTime     << std::endl
              << "min read time over " << numRanks << " ranks = " << minReadRunTime << std::endl
              << "max read time over " << numRanks << " ranks = " << maxReadRunTime << std::endl;
    const double lifComp = ( maxRunTime - minRunTime ) / minRunTime;
    std::cout << std::endl << "measured load imbalance = " << 100 * lifComp << "%" << std::endl;
  }




  MPI_Finalize();
  return 0;
}
