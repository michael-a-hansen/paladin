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

/**
 *  \file   square-matrix.hpp
 *  \date   Sep 29, 2016
 *  \author mike
 */

#ifndef SQUARE_MATRIX_HPP_
#define SQUARE_MATRIX_HPP_

#include <cassert>
#include "types.hpp"

namespace paladin {

/**
 * @brief splits a string by whitespace
 * @param str string to be split by whitespace
 * @result a vector of strings formed by splitting the input
 */
StrVecT split(const std::string& str) {
  StrVecT split;
  std::stringstream stream(std::move(str));
  std::string buffer;
  while (stream >> buffer) split.push_back(buffer);
  return split;
}

/**
 * @brief count lines in a file
 * @param file ifstream for the file
 * @result number of nonempty lines
 */
int count_nonempty_lines(std::ifstream& file) {
  int nlines = 0;
  std::string line;
  while (std::getline(file, line) && line != "") ++nlines;
  return nlines;
}

/**
 * @brief count lines in a file
 * @param filename the file name as a string
 * @result number of lines
 */
int count_nonempty_lines(const std::string& filename) {
  std::ifstream stream(filename);
  return count_nonempty_lines(stream);
}

/**
 * @class SquareMatrix
 * @brief container for a vector that supports 2D 1-based indexing
 */
struct SquareMatrix {
  const int nrows_, nnz_, nelem_;
  DVecT mat_;

  /**
   * @brief constructing a SquareMatrix and allocating storage
   * @param nrows number of rows in the matrix
   * @param nnz number of nonzero elements in the matrix
   */
  SquareMatrix(const int nrows, const int nnz)
      : nrows_(nrows), nnz_(nnz), nelem_(nrows * nrows), mat_(nelem_, 0.0) {}

  /**
   * @brief parentheses operator for accessing a square matrix element with a 2D
   * index
   * @param i row index (1-based)
   * @param j column index (1-based)
   * @result reference to the element at the i-th row and j-th column
   */
  double& operator()(const int i, const int j) {
    return mat_[(i - 1) * nrows_ + (j - 1)];
  }

  /**
   * @brief print the matrix to stdout in square form
   */
  void print() const {
    std::cout << "[";
    for (int i = 0; i < nrows_; ++i) {
      for (int j = 0; j < nrows_; ++j)
        std::cout << " " << mat_[i * nrows_ + j] << ",";
      if (i == nrows_ - 1)
        std::cout << " ]" << '\n';
      else
        std::cout << '\n';
    }
  }

  /**
   * @brief read the matrix from a matrix market file
   * @param file ifstream for the matrix market file
   */
  void read_matrix(std::ifstream& file) {
    const int maxIdxDigits = 5;
    const int maxDoubleDigits = 16;
    const int idxSize = maxIdxDigits + 1;
    const int doubleSize = maxDoubleDigits + 16;
    char ic[idxSize], jc[idxSize], vc[doubleSize];

    for (int i = 0; i < nnz_; ++i) {
      file.getline(ic, idxSize, ' ');
      file.getline(jc, idxSize, ' ');
      file.getline(vc, doubleSize, '\n');
      const int ii = atoi(ic);
      const int jj = atoi(jc);
      const double vv = atof(vc);
      mat_[(ii - 1) * nrows_ + (jj - 1)] = vv;
    }
  }
};

/**
 * @brief read the header (first line) of the matrix market exchange format
 * @param file ifstream for the matrix market file
 *
 * IMPORTANT: matrix market header must only be one line!
 */
void read_header(std::ifstream& file) {
  std::string head;
  std::getline(file, head);
  StrVecT Header = split(head);
  assert(Header[0] == "MatrixMarket" || Header[0] == "%%MatrixMarket");
  assert(Header[1] == "matrix");
  assert(Header[2] == "coordinate");
  assert(Header[3] == "real");
  assert(Header[4] == "general");
}

/**
 * @brief allocate storage for a matrix
 * @param file ifstream for the matrix market file
 * @result a SquareMatrix with enough space allocated
 */
SquareMatrix allocate_matrix(std::ifstream& file) {
  std::string line2;
  std::getline(file, line2);
  StrVecT Line2 = split(line2);

  int nrows = std::stoi(Line2[0]);
  int nnz = std::stoi(Line2[2]);

  assert(nrows == std::stoi(Line2[1]));
  assert(nnz <= nrows * std::stoi(Line2[1]));

  return SquareMatrix(nrows, nnz);
}

/**
 * @brief read a SquareMatrix from a matrix market file
 * @param matrixPath the path to the matrix as a string
 * @result a SquareMatrix populated with the matrix market data
 */
SquareMatrix read_matrix_from_mm_file(const std::string& matrixPath) {
  std::ifstream file(matrixPath);
  if (!file) {
    std::cerr << "Couldn't open the matrix file, " << matrixPath << ", exiting!"
              << '\n';
    exit(1);
  }
  read_header(file);
  SquareMatrix A = allocate_matrix(file);
  A.read_matrix(file);
  return A;
}

/**
 * @brief find the matrix dimension from a matrix market file
 * @param matrixPath the path to the matrix as a string
 * @result matrix dimension (number of rows)
 */
int read_matrix_dimension_from_mm_file(const std::string& matrixPath) {
  std::ifstream file(matrixPath);
  if (!file) {
    std::cerr << "Couldn't open the matrix file, " << matrixPath << ", exiting!"
              << '\n';
    exit(1);
  }
  std::string line;
  std::getline(file, line);
  std::getline(file, line);

  StrVecT Line = split(line);
  int nrows = std::stoi(Line[0]);

  assert(nrows == std::stoi(Line[1]));
  return nrows;
}

/**
 * @brief find the number of nonzeros of a matrix from a matrix market file
 * @param matrixPath the path to the matrix as a string
 * @result matrix nonzero count
 */
int read_matrix_nnzeros_from_mm_file(const std::string& matrixPath) {
  std::ifstream file(matrixPath);
  if (!file) {
    std::cerr << "Couldn't open the matrix file, " << matrixPath << ", exiting!"
              << '\n';
    exit(1);
  }
  std::string line;
  std::getline(file, line);
  std::getline(file, line);

  StrVecT Line = split(line);
  int nnz = std::stoi(Line[2]);
  assert(Line[0] == Line[1]);

  return nnz;
}

/**
 * @brief compute the predicted load measure of a matrix
 * @param matrixPath the path to the matrix as a string
 * @param type LoadPredictor type used to estimate the matrix load
 * @result predicted load measure
 */
MeasureT obtain_matrix_measure(
    const std::string& matrixPath,
    const LoadPredictor& type = LoadPredictor::DIMENSION) {
  MeasureT measure = 0;
  MeasureT dim = (MeasureT)read_matrix_dimension_from_mm_file(matrixPath);
  MeasureT nnz = (MeasureT)read_matrix_nnzeros_from_mm_file(matrixPath);
  try {
    switch (type) {
      case LoadPredictor::DIMENSION:
        measure = dim;
        break;
      case LoadPredictor::NUMNONZEROS:
        measure = nnz;
        break;
      case LoadPredictor::DIMCUBED:
        measure = dim * dim * dim;
        break;
      case LoadPredictor::NNZDIMSQRD:
        measure = dim * dim * nnz;
        break;
      case LoadPredictor::SPARSITY:
        measure = nnz / dim / dim;
        break;
      case LoadPredictor::SPARSENCUBED:
        measure = nnz * dim;
        break;
      default:
        throw std::invalid_argument("Invalid load type.");
    }
  } catch (const std::invalid_argument& badarg) {
    std::cerr << "Invalid argument: " << badarg.what() << '\n';
  }
  return measure;
}

}  // namespace paladin

#endif /* SQUARE_MATRIX_HPP_ */
