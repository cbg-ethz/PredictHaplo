// PredictHaplo-Paired: a program for estimating haplotypes from "next
// generation sequencing" reads
//
// Copyright (C) 2022 ETH Zurich
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

#include <cstdlib>

namespace phaplo {
struct ParseStatistics {
  struct {
    /**
     * @brief The total number of discovered reads.
     */
    std::size_t total_count;

    /**
     * @brief Number of reads that were discarded because they are unpaired.
     */
    std::size_t unpaired;

    /**
     * @brief Number of reads that were discarded because they are unmapped.
     */
    std::size_t unmapped;

    /**
     * @brief Number of reads that were discarded because of an unsupported
     * attribute (bit >= 0x100).
     */
    std::size_t unsupported_attribute;
  } reads;

  struct {
    /**
     * @brief The total number of discovered pairs.
     */
    std::size_t total_count;

    /**
     * @brief Number of pairs that were discarded because the sequence is too
     * short.
     */
    std::size_t sequence_too_short;

    /**
     * @brief Number of pairs that were discarded because the gap fraction is
     * too high.
     */
    std::size_t gap_fraction_too_high;

    /**
     * @brief Number of pairs that were discarded because the alignment score
     * fraction is too low.
     */
    std::size_t align_score_fraction_too_low;
  } pairs;
};
} // namespace phaplo