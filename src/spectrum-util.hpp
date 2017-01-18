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
 *  \file   spectrum-util.hpp
 *  \date   Sep 29, 2016
 *  \author mike
 */

#ifndef SPECTRUM_UTIL_HPP_
#define SPECTRUM_UTIL_HPP_

#include "types.hpp"

namespace paladin {


struct EigenSorter {
  virtual ~EigenSorter() {}
  virtual bool operator()(const EigenvalueT& L, const EigenvalueT& R) const {
    return (std::norm(L) > std::norm(R));
  }
};
struct LReal : public EigenSorter {
  bool operator()(const EigenvalueT& L, const EigenvalueT& R) const {
    return (L.real() > R.real());
  }
};
struct SReal : public EigenSorter {
  bool operator()(const EigenvalueT& L, const EigenvalueT& R) const {
    return (L.real() < R.real());
  }
};
struct LImag : public EigenSorter {
  bool operator()(const EigenvalueT& L, const EigenvalueT& R) const {
    return (L.imag() > R.imag());
  }
};
struct SImag : public EigenSorter {
  bool operator()(const EigenvalueT& L, const EigenvalueT& R) const {
    return (L.imag() < R.imag());
  }
};
struct LMag : public EigenSorter {
  bool operator()(const EigenvalueT& L, const EigenvalueT& R) const {
    return (std::norm(L) > std::norm(R));
  }
};
struct SMag : public EigenSorter {
  bool operator()(const EigenvalueT& L, const EigenvalueT& R) const {
    return (std::norm(L) < std::norm(R));
  }
};

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
void sort_spectrum(SpectrumT& s, const std::string sortStr) {

  std::unique_ptr<EigenSorter> sorter;
  if (sortStr == "SM") {
    sorter = paladin::make_unique<SMag>();
  } else if (sortStr == "LM") {
    sorter = paladin::make_unique<LMag>();
  } else if (sortStr == "SR") {
    sorter = paladin::make_unique<SReal>();
  } else if (sortStr == "LR") {
    sorter = paladin::make_unique<LReal>();
  } else if (sortStr == "SI") {
    sorter = paladin::make_unique<SImag>();
  } else if (sortStr == "LI") {
    sorter = paladin::make_unique<LImag>();
  } else {
    sorter = paladin::make_unique<LMag>();
  }
  std::sort(s.begin(), s.end(), *sorter);
}

/**
 * @brief print the spectrum to stdout
 * @param s const reference to a spectrum object
 */
void show_spectrum(const SpectrumT& s) {
  int idx = -1;
  for (auto eig : s) {
    ++idx;
    std::cout << "lambda_" << idx << " = " << eig.real() << " + " << eig.imag()
              << "j" << '\n';
  }
}

/**
 * @brief write the spectrum to a file
 * @param s const reference to a spectrum object
 * @param realout ofstream reference to the file
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
void write_spectrum(const SpectrumT& s, std::ofstream&& eigsout) {
  for (auto eig : s) {
    eigsout << eig.real() << " " << eig.imag() << '\n';
  }
}

}  // namespace paladin

#endif /* SPECTRUM_UTIL_HPP_ */
