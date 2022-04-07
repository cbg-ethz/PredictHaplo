// PredictHaplo-Paired: a program for estimating haplotypes from "next
// generation sequencing" reads
//
// Copyright (C) 2014 Volker Roth
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <phaplo/binary.hpp>

namespace phaplo {
std::string binary(int number, std::stringstream &strs) {
  int remainder;

  if (number <= 1) {
    strs << number;
    return strs.str();
  }

  remainder = number % 2;
  binary(number >> 1, strs);
  strs << remainder;

  return strs.str();
}
}
