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

#include <make-unique.hpp>

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
  std::vector<std::vector<std::complex<double> > > leftEigenvectors_,
      rightEigenvectors_;

  const double vectorWriteThreshold_ = 1e-3;

  /**
   * @class EigenSorter
   * Base class for sorting eigenpairs by real/imaginary parts or magnitude
   */
  struct EigenSorter {
    virtual ~EigenSorter() {}
    virtual bool operator()( const EigenPairT& L, const EigenPairT& R ) const {
      return ( std::norm( L.first ) > std::norm( R.first ) );
    }
  };

  /**
   * @class LReal
   * Comparator class for sorting eigenpairs by largest real part
   */
  struct LReal : public EigenSorter {
    bool operator()( const EigenPairT& L, const EigenPairT& R ) const {
      return ( L.first.real() > R.first.real() );
    }
  };

  /**
   * @class SReal
   * Comparator class for sorting eigenpairs by smallest real part
   */
  struct SReal : public EigenSorter {
    bool operator()( const EigenPairT& L, const EigenPairT& R ) const {
      return ( L.first.real() < R.first.real() );
    }
  };

  /**
   * @class LImag
   * Comparator class for sorting eigenpairs by largest imaginary part
   */
  struct LImag : public EigenSorter {
    bool operator()( const EigenPairT& L, const EigenPairT& R ) const {
      return ( L.first.imag() > R.first.imag() );
    }
  };

  /**
   * @class SImag
   * Comparator class for sorting eigenpairs by smallest imaginary part
   */
  struct SImag : public EigenSorter {
    bool operator()( const EigenPairT& L, const EigenPairT& R ) const {
      return ( L.first.imag() < R.first.imag() );
    }
  };

  /**
   * @class LMag
   * Comparator class for sorting eigenpairs by largest magnitude
   */
  struct LMag : public EigenSorter {
    bool operator()( const EigenPairT& L, const EigenPairT& R ) const {
      return ( std::norm( L.first ) > std::norm( R.first ) );
    }
  };

  /**
   * @class SMag
   * Comparator class for sorting eigenpairs by smallest magnitude
   */
  struct SMag : public EigenSorter {
    bool operator()( const EigenPairT& L, const EigenPairT& R ) const {
      return ( std::norm( L.first ) < std::norm( R.first ) );
    }
  };

  /**
   * @class LMagComplexValue
   * Comparator class for sorting complex values by largest magnitude
   */
  struct LMagComplexValue : public EigenSorter {
    bool operator()( const EigenvalueT& L, const EigenvalueT& R ) const {
      return ( std::norm( L ) > std::norm( R ) );
    }
  };

 public:
  /**
   * @brief construct a SpectralDecomposition object
   * @param N the size of the matrix
   *
   */
  SpectralDecomposition( const int N, const double threshold = 1e-3 )
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
  void append_left_eigenvector(
      const std::vector<std::complex<double> > eigenvectorL ) {
    leftEigenvectors_.push_back( eigenvectorL );
  }

  /**
   * @brief append a right eigenvector to the decomposition
   * @param eigenvectorR the right eigenvector to append
   */
  void append_right_eigenvector(
      const std::vector<std::complex<double> > eigenvectorR ) {
    rightEigenvectors_.push_back( eigenvectorR );
  }

  /**
   * @brief sorting a spectrum in a variety of ways
   * @param s reference to a spectrum object which is sorted in place
   * @param sortStr a string indicating the type of sort (see below)
   *
   * Sort strings:
   * - SM: smallest magnitude
   * - LM: largest magnitude
   * - SR: smallest real part
   * - LR: largest real part
   * - SI: smallest imaginary part
   * - LI: largest imaginary part
   */
  void sort( const std::string sortStr,
             const bool doLeft,
             const bool doRight ) {
    std::unique_ptr<EigenSorter> sorter;
    if ( sortStr == "SM" ) {
      sorter = paladin::make_unique<SMag>();
    } else if ( sortStr == "LM" ) {
      sorter = paladin::make_unique<LMag>();
    } else if ( sortStr == "SR" ) {
      sorter = paladin::make_unique<SReal>();
    } else if ( sortStr == "LR" ) {
      sorter = paladin::make_unique<LReal>();
    } else if ( sortStr == "SI" ) {
      sorter = paladin::make_unique<SImag>();
    } else if ( sortStr == "LI" ) {
      sorter = paladin::make_unique<LImag>();
    } else {
      sorter = paladin::make_unique<LMag>();
    }

    std::vector<EigenPairT> eigenpairs;
    for ( std::size_t i = 0; i < eigenvalues_.size(); ++i ) {
      if ( doLeft && doRight ) {
        eigenpairs.push_back( std::make_pair(
            eigenvalues_[i],
            std::make_pair( leftEigenvectors_[i], rightEigenvectors_[i] ) ) );
      } else if ( doLeft ) {
        std::vector<std::complex<double> > null( eigenvalues_.size() );
        eigenpairs.push_back( std::make_pair(
            eigenvalues_[i], std::make_pair( leftEigenvectors_[i], null ) ) );
      } else if ( doRight ) {
        std::vector<std::complex<double> > null( eigenvalues_.size() );
        eigenpairs.push_back( std::make_pair(
            eigenvalues_[i], std::make_pair( null, rightEigenvectors_[i] ) ) );
      } else {
        std::vector<std::complex<double> > null( eigenvalues_.size() );
        eigenpairs.push_back(
            std::make_pair( eigenvalues_[i], std::make_pair( null, null ) ) );
      }
    }
    std::sort( eigenpairs.begin(), eigenpairs.end(), *sorter );

    for ( std::size_t i = 0; i < eigenvalues_.size(); ++i ) {
      eigenvalues_[i] = eigenpairs[i].first;
      if ( doLeft && doRight ) {
        rightEigenvectors_[i] = eigenpairs[i].second.first;
        leftEigenvectors_[i] = eigenpairs[i].second.second;
      } else if ( doLeft ) {
        leftEigenvectors_[i] = eigenpairs[i].second.second;
      } else if ( doRight ) {
        rightEigenvectors_[i] = eigenpairs[i].second.first;
      } else {
      }
    }
  }

  /**
   * @brief print the spectrum to stdout
   */
  void show_eigenvalues() const {
    int idx = -1;
    for ( const auto& eig : eigenvalues_ ) {
      ++idx;
      std::cout << "lambda_" << idx << " = " << eig.real() << " + "
                << eig.imag() << "j" << '\n';
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
    for ( const auto& eig : eigenvalues_ ) {
      eigsout << eig.real() << " " << eig.imag() << '\n';
    }
  }

  /**
   * @brief write the eigenvectors to a file
   * @param leftout ofstream reference to the file for left eigenvectors
   * @param rightout ofstream reference to the file for right eigenvectors
   *
   * Document me when I'm working!
   */
  void write_eigenvectors( std::ofstream&& leftout,
                           std::ofstream&& rightout,
                           bool writeLeft,
                           bool writeRight ) const {
    const int N = static_cast<int>( eigenvalues_.size() );

    if ( writeLeft ) {
      leftout << "%%MatrixMarket matrix coordinate complex general\n";
      // iterate through eigenvectors
      for ( int i = 0; i < N; ++i ) {
        // iterate through each eigenvector
        double maxElement = std::norm( *std::max_element(
            leftEigenvectors_[i].begin(), leftEigenvectors_[i].end(),
            LMagComplexValue() ) );
        for ( int j = 0; j < N; ++j ) {
          const std::complex<double> leij = leftEigenvectors_[i][j];
          // only print if absolute value beats threshold
          if ( std::norm( leij ) > vectorWriteThreshold_ * maxElement ) {
            // write element of eigenvector matrix to file
            leftout << j + 1 << " " << i + 1 << " " << leij.real() << " "
                    << leij.imag() << '\n';
          }
        }
      }
    }
    if ( writeRight ) {
      rightout << "%%MatrixMarket matrix coordinate complex general\n";
      // iterate through eigenvectors
      for ( int i = 0; i < N; ++i ) {
        // iterate through each eigenvector
        double maxElement = std::norm( *std::max_element(
            rightEigenvectors_[i].begin(), rightEigenvectors_[i].end(),
            LMagComplexValue() ) );
        for ( int j = 0; j < N; ++j ) {
          const std::complex<double> reij = rightEigenvectors_[i][j];
          // only print if absolute value beats threshold
          if ( std::norm( reij ) > vectorWriteThreshold_ * maxElement ) {
            // write element of eigenvector matrix to file
            rightout << j + 1 << " " << i + 1 << " " << reij.real() << " "
                     << reij.imag() << '\n';
          }
        }
      }
    }
  }
};

}  // namespace paladin

#endif /* SRC_SPECTRUM_UTIL_HPP_ */
