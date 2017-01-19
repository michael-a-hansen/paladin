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

#ifndef SRC_HEADER_HPP_
#define SRC_HEADER_HPP_

#include <iostream>

namespace paladin {

void print_header() {
  std::cout << '\n' << '\n';
  std::cout
      << " ____   ____  _       ____  ___    ____  ____  " << '\n'
      << "|    \\ /    || |     /    ||   \\  |    ||    \\ " << '\n'
      << "|  o  |  o  || |    |  o  ||    \\  |  | |  _  |" << '\n'
      << "|   _/|     || |___ |     ||  D  | |  | |  |  |" << '\n'
      << "|  |  |  _  ||     ||  _  ||     | |  | |  |  |" << '\n'
      << "|  |  |  |  ||     ||  |  ||     | |  | |  |  |" << '\n'
      << "|__|  |__|__||_____||__|__||_____||____||__|__|" << '\n'
      << '\n'
      << "PArallel LApack DIstributor for many serial eigendecompositions."
      << '\n'
      << "Copyright (c) 2016, 2017 Michael A. Hansen" << '\n'
      << "------------------------------------------------" << '\n';
  std::cout << '\n' << '\n';
}

}  // namespace paladin

#endif /* SRC_HEADER_HPP_ */
