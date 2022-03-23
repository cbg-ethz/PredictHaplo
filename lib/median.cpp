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

#include <algorithm>
#include <limits>
#include <phaplo/median.hpp>

namespace phaplo {

double median(std::vector<int> v) noexcept {
  if (v.empty())
    return std::numeric_limits<double>::quiet_NaN();

  const auto n = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + n, v.end());

  if (v.size() % 2 == 1)
    return v[n];
  std::nth_element(v.begin(), v.begin() + n - 1, v.begin() + n);
  return 0.5 * (v[n] + v[n - 1]);
}

} // namespace phaplo