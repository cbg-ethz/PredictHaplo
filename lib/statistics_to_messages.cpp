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

#include "phaplo/statistics_to_messages.hpp"
#include "phaplo/ParseStatistics.hpp"
#include <iomanip>
#include <string>
#include <vector>

namespace phaplo {

namespace {
std::string to_percentage(const std::size_t value, const std::size_t total) {
  std::stringstream stream;
  stream << std::setprecision(1) << std::fixed;
  stream << 100. * static_cast<double>(value) / static_cast<double>(total)
         << "\%";
  return stream.str();
}
} // namespace

std::vector<std::string>
statistics_to_messages(const ParseStatistics &statistics) noexcept {
  auto r = std::vector<std::string>();

  const auto to_reads_percentage = [&](const std::size_t x) {
    return to_percentage(x, statistics.reads.total_count);
  };
  const auto to_pairs_percentage = [&](const std::size_t x) {
    return to_percentage(x, statistics.pairs.total_count);
  };

  if (statistics.reads.unsupported_attribute) {
    r.push_back(to_reads_percentage(statistics.reads.unsupported_attribute) +
                " of the reads were discarded because of an unsupported "
                "attribute.");
  }

  if (statistics.reads.unmapped) {
    r.push_back(to_reads_percentage(statistics.reads.unmapped) +
                " of the reads were discarded because they are unmapped.");
  }

  if (statistics.reads.unpaired) {
    r.push_back(to_reads_percentage(statistics.reads.unpaired) +
                " of the reads were discarded because they are unpaired.");
  }

  if (statistics.pairs.sequence_too_short) {
    r.push_back(to_pairs_percentage(statistics.pairs.sequence_too_short) +
                " of the read pairs were discarded because the sequence is too "
                "short. The flag \"--min_length\" can be used to configure the "
                "minimum length.");
  }

  if (statistics.pairs.gap_fraction_too_high) {
    r.push_back(to_pairs_percentage(statistics.pairs.gap_fraction_too_high) +
                " of the read pairs were discarded because the gap fraction "
                "is too high. The flag \"--max_gap_fraction\" can be used to "
                "configure the maximum gap fraction.");
  }

  if (statistics.pairs.align_score_fraction_too_low) {
    r.push_back(
        to_pairs_percentage(statistics.pairs.align_score_fraction_too_low) +
        " of the read pairs were discarded because the alignment score "
        "fraction (scaled value of TAG \"AS\") is too low. The flag "
        "\"--min_align_score_fraction\" can be used to configure the minimum "
        "alignment score fraction.");
  }

  return r;
}
} // namespace phaplo