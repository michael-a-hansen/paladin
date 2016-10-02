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
 *  \file   types.hpp
 *  \date   Sep 29, 2016
 *  \author mike
 */

#ifndef TYPES_HPP_
#define TYPES_HPP_


namespace paladin
{


  typedef std::vector<double> DVecT;
  typedef std::vector<std::string> StrVecT;

  typedef double MeasureT;
  typedef std::pair<std::string, MeasureT> NameMeasurePairT;

  typedef std::complex<double> EigenvalueT;
  typedef std::vector<EigenvalueT> SpectrumT;



  /**
   * @enum LoadPredictor
   * @brief Types of predicting a matrix load
   */
  enum class LoadPredictor
  {
    DIMENSION,
    NUMNONZEROS,
    DIMCUBED,
    NNZDIMSQRD,
    SPARSITY,
    SPARSENCUBED
  };

  /**
   * @brief reading abbreviated input strings into LoadPredictor types
   * @param str a std::string of abbreviated format specifying a LoadPredictor type
   * @result the LoadPredictor type corresponding to the input string
   */
  LoadPredictor string_to_measure_type( const std::string& str )
  {
    LoadPredictor type;
    try{
      if     ( str == "dim" ) type = LoadPredictor::DIMENSION;
      else if( str == "nnz" ) type = LoadPredictor::NUMNONZEROS;
      else if( str == "dcb" ) type = LoadPredictor::DIMCUBED;
      else if( str == "zds" ) type = LoadPredictor::NNZDIMSQRD;
      else if( str == "sps" ) type = LoadPredictor::SPARSITY;
      else if( str == "spc" ) type = LoadPredictor::SPARSENCUBED;
      else throw std::invalid_argument( "Invalid load string." );
    }
    catch( const std::invalid_argument& badarg ) {
      std::cerr << "Invalid argument: " << badarg.what() << '\n';
    }
    return type;
  }

  /**
   * @brief converting a LoadPredictor type to a descriptive string
   * @param type the LoadPredictor type
   * @result the (string) description of the input LoadPredictor type
   */
  std::string measure_type_description( const LoadPredictor& type )
  {
    std::string str;
    try{
      if     ( type == LoadPredictor::DIMENSION )    str = "matrix dimension";
      else if( type == LoadPredictor::NUMNONZEROS )  str = "number of nonzeros";
      else if( type == LoadPredictor::DIMCUBED )     str = "matrix dimension cubed";
      else if( type == LoadPredictor::NNZDIMSQRD )   str = "number nonzeros * dim squared";
      else if( type == LoadPredictor::SPARSITY )     str = "matrix sparsity = nnz / dim / dim";
      else if( type == LoadPredictor::SPARSENCUBED ) str = "sparsity * dim cubed";
      else throw std::invalid_argument( "Invalid load type." );
    }
    catch( const std::invalid_argument& badarg ) {
      std::cerr << "Invalid argument: " << badarg.what() << '\n';
    }
    return str;
  }

  /**
   * @brief sorting a vector of NameMeasurePairs by the load measure
   * @param flock reference to a vector of NameMeasurePairs which is sorted in place
   */
  void sort_name_measure_pairs( std::vector<NameMeasurePairT>& flock )
  {
    struct MeasureComparator{
      inline bool operator()( const NameMeasurePairT& left, const NameMeasurePairT& right ){
        return ( left.second > right.second );
      }
    };
    std::sort( flock.begin(), flock.end(), MeasureComparator() );
  }


} // namespace paladin


#endif /* TYPES_HPP_ */
