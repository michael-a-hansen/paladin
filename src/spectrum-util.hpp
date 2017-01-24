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

#ifndef SRC_SPECTRUM_UTIL_HPP_
#define SRC_SPECTRUM_UTIL_HPP_

namespace paladin {

using EigenvalueT = std::complex<double>;
using SpectrumT = std::vector<EigenvalueT>;
using ComplexVecT = std::vector<std::complex<double> >;
using RLPairT = std::pair<ComplexVecT, ComplexVecT>;
using EigenPairT = std::pair<EigenvalueT, RLPairT>;

/**
 * @class SpectralDecomposition
 *
 * Stores and sorts eigenvalues, eigenvectors
 */
class SpectralDecomposition {
 protected:
  std::vector<std::complex<double> > eigenvalues_;
  std::vector<std::vector<std::complex<double> > > leftEigenvectors_, rightEigenvectors_;

  const double vectorWriteThreshold_ = 1.0e-3;

  /**
   * @class LMagComplexValue
   * Comparator class for sorting complex values by largest magnitude
   */
  struct LMagComplexValue {
    bool operator()( const EigenvalueT& L, const EigenvalueT& R ) const {
      return ( std::norm( L ) < std::norm( R ) );
    }
  };

 public:
  /**
   * @brief construct a SpectralDecomposition object
   * @param N the size of the matrix
   *
   */
  SpectralDecomposition( const int N, const double threshold = 1.0e-3 )
      : vectorWriteThreshold_( threshold ) {
    eigenvalues_.reserve( N );
    leftEigenvectors_.reserve( N );
    rightEigenvectors_.reserve( N );
  }

  /**
   * @brief append an eigenvalue to the spectrum
   * @param eigenvalue the eigenvalue to append
   */
  void append_eigenvalue( const std::complex<double> eigenvalue ) {
    eigenvalues_.push_back( eigenvalue );
  }

  /**
   * @brief append a left eigenvector to the decomposition
   * @param eigenvectorL the left eigenvector to append
   */
  void append_left_eigenvector( const std::vector<std::complex<double> > eigenvectorL ) {
    leftEigenvectors_.push_back( eigenvectorL );
  }

  /**
   * @brief append a right eigenvector to the decomposition
   * @param eigenvectorR the right eigenvector to append
   */
  void append_right_eigenvector( const std::vector<std::complex<double> > eigenvectorR ) {
    rightEigenvectors_.push_back( eigenvectorR );
  }

  /**
   * @brief print the spectrum to stdout
   */
  void show_eigenvalues() const {
    int idx = -1;
    for ( const auto& eig : eigenvalues_ ) {
      ++idx;
      std::cout << "lambda_" << idx << " = " << eig.real() << " + " << eig.imag() << "j" << '\n';
    }
  }

  /**
   * @brief write the spectrum to a file
   * @param eigsout ofstream reference to the file
   *
   * This writes the spectrum to a file in two space-separated columns,
   * with the real part in the first and imaginary part in the second
   *
   * realpart_0 imagpart_0
   * realpart_1 imagpart_1
   * .
   * .
   * .
   */
  void write_eigenvalues( std::ofstream&& eigsout ) const {
    const int N = static_cast<int>( eigenvalues_.size() );
    for ( int i = 0; i < N; ++i ) {
      eigsout << eigenvalues_[i].real() << " " << eigenvalues_[i].imag() << '\n';
    }
  }

  /**
   * @brief write the eigenvectors to a file
   * @param leftout ofstream reference to the file for left eigenvectors
   *
   * This writes out the eigenvectors as a matrix market format file.
   * All elements whose magnitude is within a threshold of the largest element
   * magnitude are written.
   */
  void write_left_eigenvectors( std::ofstream&& leftout ) const {
    const int N = static_cast<int>( eigenvalues_.size() );

    int nnzL = 0;
    for ( int i = 0; i < N; ++i ) {
      double maxL = std::norm( *std::max_element(
          leftEigenvectors_[i].begin(), leftEigenvectors_[i].end(), LMagComplexValue() ) );
      for ( int j = 0; j < N; ++j ) {
        const std::complex<double> leij = leftEigenvectors_[i][j];
        if ( std::norm( leij ) > vectorWriteThreshold_ * maxL ) {
          nnzL++;
        }
      }
    }
    leftout << "%%MatrixMarket matrix coordinate complex general\n";
    leftout << N << " " << N << " " << nnzL << '\n';
    for ( int i = 0; i < N; ++i ) {
      double maxElement = std::norm( *std::max_element(
          leftEigenvectors_[i].begin(), leftEigenvectors_[i].end(), LMagComplexValue() ) );
      for ( int j = 0; j < N; ++j ) {
        const std::complex<double> leij = leftEigenvectors_[i][j];
        if ( std::norm( leij ) > vectorWriteThreshold_ * maxElement ) {
          leftout << j + 1 << " " << i + 1 << " " << leij.real() << " " << leij.imag() << '\n';
        }
      }
    }
  }

  /**
   * @brief write the eigenvectors to a file
   * @param rightout ofstream reference to the file for right eigenvectors
   *
   * This writes out the eigenvectors as a matrix market format file.
   * All elements whose magnitude is within a threshold of the largest element
   * magnitude are written.
   */
  void write_right_eigenvectors( std::ofstream&& rightout ) const {
    const int N = static_cast<int>( eigenvalues_.size() );

    int nnzR = 0;
    for ( int i = 0; i < N; ++i ) {
      double maxR = std::norm( *std::max_element(
          rightEigenvectors_[i].begin(), rightEigenvectors_[i].end(), LMagComplexValue() ) );
      for ( int j = 0; j < N; ++j ) {
        const std::complex<double> reij = rightEigenvectors_[i][j];
        if ( std::norm( reij ) > vectorWriteThreshold_ * maxR ) {
          nnzR++;
        }
      }
    }
    rightout << "%%MatrixMarket matrix coordinate complex general\n";
    rightout << N << " " << N << " " << nnzR << '\n';
    for ( int i = 0; i < N; ++i ) {
      double maxElement = std::norm( *std::max_element(
          rightEigenvectors_[i].begin(), rightEigenvectors_[i].end(), LMagComplexValue() ) );
      for ( int j = 0; j < N; ++j ) {
        const std::complex<double> reij = rightEigenvectors_[i][j];
        if ( std::norm( reij ) > vectorWriteThreshold_ * maxElement ) {
          rightout << j + 1 << " " << i + 1 << " " << reij.real() << " " << reij.imag() << '\n';
        }
      }
    }
  }
};

}  // namespace paladin

#endif /* SRC_SPECTRUM_UTIL_HPP_ */
