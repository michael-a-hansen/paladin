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
 *  \file   command-line-parser.hpp
 *  \date   Sep 29, 2016
 *  \author mike
 */

#ifndef COMMAND_LINE_PARSER_HPP_
#define COMMAND_LINE_PARSER_HPP_

#include <vector>
#include <string>
#include "types.hpp"

namespace paladin
{

  /**
   * @class CommandLineParser
   * @brief tool for parsing the command line
   *
   * usage: mpirun -np 4 executable-name -option1 arg1 -option2 arg2 ...
   */
  class CommandLineParser{
  protected:
    StrVecT pieces;
  public:

    /**
     * @brief constructing a CommandLineParser from argc, argv
     * @param argc number of command line arguments, obtained through main()
     * @param argv char** of command line arguments, obtained through main()
     */
    CommandLineParser( const int argc, char *argv[] )
  {
      for( int i=1; i<argc; ++i )
        this->pieces.push_back( std::string( argv[i] ) );
  }


    /**
     * @brief check that an option was provided in the command line
     * @param name name of the option
     * @result true/false if the option was passed/not
     */
    bool checkExists( const std::string& name ) const
    {
      return std::find( pieces.begin(), pieces.end(), name ) != pieces.end();
    }

    /**
     * @brief get the value of an option provided in the command line
     * @param name name of the option
     * @param name default value for the option
     * @result the value provided for the option
     */
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

} // namespace paladin


#endif /* COMMAND_LINE_PARSER_HPP_ */
