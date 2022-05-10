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

#include <bitset>
#include <limits>

namespace phaplo {

/**
 * @brief Converts an integer to its binary representation.
 */
std::bitset<std::numeric_limits<unsigned int>::digits>
binary(unsigned int number);

/**
 * @brief Returns the number of bits required to represent `bits`.
 */
std::size_t used_bits(std::bitset<std::numeric_limits<unsigned int>::digits> bits);

} // namespace phaplo