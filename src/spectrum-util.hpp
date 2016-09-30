/*
 * The MIT License
 *
 * Copyright (c) 2016 The University of Utah
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
 *  \file   spectrum-util.hpp
 *  \date   Sep 29, 2016
 *  \author mike
 */

#ifndef SPECTRUM_UTIL_HPP_
#define SPECTRUM_UTIL_HPP_

#include "types.hpp"

namespace paladin
{




  /**
   * @brief sorting a spectrum in a variety of ways
   * @param s reference to a spectrum object which is sorted in place
   * @param sortStr a string indicating the type of sort (see below)
   *
   * Sort strings:
   * - SM: smallest by magnitude
   * - LM: largest by magnitude
   * - SR: smallest real part
   * - LR: largest real part
   * - SI: smallest imaginary part
   * - LI: largest imaginary part
   */
  void sort_spectrum( SpectrumT& s, const std::string sortStr )
  {
    struct EigenSorter{
      virtual ~EigenSorter() {}
      virtual bool operator() ( const EigenvalueT& L, const EigenvalueT& R ) const { return ( std::norm( L ) > std::norm( R ) );}
    };
    struct LReal : public EigenSorter { bool operator() ( const EigenvalueT& L, const EigenvalueT& R ) const { return ( L.real() > R.real() );} };
    struct SReal : public EigenSorter { bool operator() ( const EigenvalueT& L, const EigenvalueT& R ) const { return ( L.real() < R.real() );} };
    struct LImag : public EigenSorter { bool operator() ( const EigenvalueT& L, const EigenvalueT& R ) const { return ( L.imag() > R.imag() );} };
    struct SImag : public EigenSorter { bool operator() ( const EigenvalueT& L, const EigenvalueT& R ) const { return ( L.imag() < R.imag() );} };
    struct LMag  : public EigenSorter { bool operator() ( const EigenvalueT& L, const EigenvalueT& R ) const { return ( std::norm( L ) > std::norm( R ) );} };
    struct SMag  : public EigenSorter { bool operator() ( const EigenvalueT& L, const EigenvalueT& R ) const { return ( std::norm( L ) < std::norm( R ) );} };

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

  /**
   * @brief print the spectrum to stdout
   * @param s const reference to a spectrum object
   */
  void show_spectrum( const SpectrumT& s )
  {
    for( size_t i=0; i<s.size(); ++i )
      std::cout << "lambda_" << i << " = " << s[i].real() << " + " << s[i].imag() << "j" << '\n';
  }

  /**
   * @brief write the spectrum to two files, for the real and imaginary parts
   * @param s const reference to a spectrum object
   * @param realout ofstream reference to a file for the real part
   * @param imagout ofstream reference to a file for the imaginary part
   */
  void write_spectrum( const SpectrumT& s, std::ofstream& realout, std::ofstream& imagout )
  {
    for( size_t i=0; i<s.size(); ++i ){
      realout << s[i].real() << '\n';
      imagout << s[i].imag() << '\n';
    }
  }

} // namespace paladin


#endif /* SPECTRUM_UTIL_HPP_ */
