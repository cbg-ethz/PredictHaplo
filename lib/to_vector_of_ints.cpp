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

#include <phaplo/to_vector_of_ints.hpp>
#include <algorithm>

namespace phaplo {

std::vector<int> to_vector_of_ints(const std::string &in) {
  auto r = std::vector<int>(in.size());
  std::transform(in.begin(), in.end(), r.begin(), [](const auto c) {
    return c - '0';
  });
  return r;
}

} // namespace phaplo