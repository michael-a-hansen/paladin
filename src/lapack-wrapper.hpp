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

#ifndef SRC_LAPACK_WRAPPER_HPP_
#define SRC_LAPACK_WRAPPER_HPP_

namespace paladin {

/*
 * @brief full spectral decomposition of a general, double-precision matrix
 * @param doleft 'N' for no left eigenmatrix and 'V' for computing it, char ref.
 * @param doright 'N' for no right eigenmatrix and 'V' for computing it, char ref.
 * @param N number of rows in A, integer reference
 * @param A a matrix in column-major form, double array
 * @param lda leading dimension of A, integer reference
 * @param real eigenvalue real parts, double array
 * @param imag eigenvalue imaginary parts, double array
 * @param leftmat left eigenmatrix, double array
 * @param ldleft leading dim. of left eig matrix, integer reference
 * @param rightmat right eigenmatrix, double array
 * @param ldright leading dim. of right eig matrix, integer reference
 * @param work work array size, double reference
 * @param lwork work array, double array
 * @param info info structure, integer reference
 */
extern "C" void dgeev_( const char* doleft,
                        const char* doright,
                        int* N,
                        double* A,
                        int* lda,
                        double* real,
                        double* imag,
                        double* leftmat,
                        int* ldleft,
                        double* rightmat,
                        int* ldright,
                        double* work,
                        int* lwork,
                        int* info );

}  // namespace paladin

#endif /* SRC_LAPACK_WRAPPER_HPP_ */
