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

#pragma once

#include <vector>

namespace phaplo {
/**
 * @brief Returns the median of the values in `v`.
 * @remark Returns NaN for an empty vector `v`.
 * @remark Reorders `v`.
 */
double median(std::vector<int> v) noexcept;
} // namespace phaplo