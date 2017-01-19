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

#ifndef SRC_COMMAND_LINE_PARSER_HPP_
#define SRC_COMMAND_LINE_PARSER_HPP_

#include <string>
#include <types.hpp>
#include <vector>

namespace paladin {

/**
 * @brief split a string by a delimiter
 * @param str the string to split
 * @param delimiter the delimiter [default is a space]
 */
std::vector<std::string> split_string( const std::string& str,
                                       char delimiter = ' ' ) {
  std::stringstream ss;
  ss.str( str );
  std::string item;
  std::vector<std::string> elems;
  while ( std::getline( ss, item, delimiter ) ) {
    if ( !item.empty() ) {
      elems.push_back( item );
    }
  }
  return elems;
}

/**
 * @class CommandLineParser
 * @brief tool for parsing the command line
 *
 * usage: executable-name --key1=value1 --key2=value2 ...
 *
 * - Two dashes must be used
 * - An equals sign must be used
 * - If a key is written multiple times, only the first is used
 */
class CommandLineParser {
 protected:
  std::vector<std::string> keys_;
  std::vector<std::string> vals_;

 public:
  /**
   * @brief constructing a CommandLineParser from argc, argv
   * @param argc number of command line arguments, obtained through main()
   * @param argv char** of command line arguments, obtained through main()
   */
  CommandLineParser( const int argc, char* argv[] ) {
    for ( int i = 1; i < argc; ++i ) {
      std::string inputPiece = static_cast<std::string>( argv[i] );
      inputPiece.erase( 0, 2 );
      std::vector<std::string> pieces = split_string( inputPiece, '=' );
      keys_.push_back( pieces[0] );
      if ( pieces.size() > 1 ) {
        vals_.push_back( pieces[1] );
      } else {
        vals_.push_back( "n/a" );
      }
    }
  }

  /**
   * @brief check that an option was provided in the command line
   * @param name name of the option
   * @result true/false if the option was passed/not
   */
  bool checkExists( const std::string& name ) const {
    return std::find( keys_.begin(), keys_.end(), name ) != keys_.end();
  }

  /**
   * @brief get the value of an option provided in the command line
   * @param name name of the option
   * @param name default value for the option
   * @result the value provided for the option
   */
  std::string getValue( const std::string& name,
                        const std::string& defaultStr ) const {
    if ( checkExists( name ) ) {
      auto itr = std::find( keys_.begin(), keys_.end(), name );
      int idx = itr - keys_.begin();
      if ( itr != keys_.end() ) {
        return vals_[idx];
      } else
        return defaultStr;
    } else
      return defaultStr;
  }

  /**
   * @brief obtain the keys
   */
  std::vector<std::string> get_keys() { return keys_; }

  /**
   * @brief obtain the values
   */
  std::vector<std::string> get_values() { return vals_; }
};

}  // namespace paladin

#endif /* SRC_COMMAND_LINE_PARSER_HPP_ */
