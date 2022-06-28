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

#include "phaplo/ParseStatistics.hpp"
#include "phaplo/statistics_to_messages.hpp"
#include <gtest/gtest.h>
#include <vector>

namespace phaplo {
namespace {
struct statistics_to_messages_Test : ::testing::Test {};

TEST_F(statistics_to_messages_Test, returns_no_messages_for_empty_statistics) {
  const auto statistics = ParseStatistics();
  EXPECT_TRUE(statistics_to_messages(statistics).empty());
}

TEST_F(statistics_to_messages_Test,
       returns_correct_message_for_unmapped_reads) {
  auto statistics = ParseStatistics();
  statistics.reads.total_count = 1000;
  statistics.reads.unmapped = 501;

  const auto messages = statistics_to_messages(statistics);
  ASSERT_EQ(messages.size(), 1);
  EXPECT_EQ(messages.front(),
            "50.1\% of the reads were discarded because they are unmapped.");
}

TEST_F(statistics_to_messages_Test,
       returns_correct_message_for_unpaired_reads) {
  auto statistics = ParseStatistics();
  statistics.reads.total_count = 1000;
  statistics.reads.unpaired = 501;

  const auto messages = statistics_to_messages(statistics);
  ASSERT_EQ(messages.size(), 1);
  EXPECT_EQ(messages.front(),
            "50.1\% of the reads were discarded because they are unpaired.");
}

TEST_F(statistics_to_messages_Test,
       returns_correct_message_for_unsupported_attributes) {
  auto statistics = ParseStatistics();
  statistics.reads.total_count = 1000;
  statistics.reads.unsupported_attribute = 501;

  const auto messages = statistics_to_messages(statistics);
  ASSERT_EQ(messages.size(), 1);
  EXPECT_EQ(messages.front(),
            "50.1\% of the reads were discarded because of an unsupported "
            "attribute.");
}

TEST_F(statistics_to_messages_Test,
       returns_correct_message_for_too_short_sequences) {
  auto statistics = ParseStatistics();
  statistics.pairs.total_count = 1000;
  statistics.pairs.sequence_too_short = 501;

  const auto messages = statistics_to_messages(statistics);
  ASSERT_EQ(messages.size(), 1);
  EXPECT_EQ(messages.front(),
            "50.1\% of the read pairs were discarded because the sequence is "
            "too short. The flag \"--min_length\" can be used to configure the "
            "minimum length.");
}

TEST_F(statistics_to_messages_Test,
       returns_correct_message_for_too_high_gap_fraction) {
  auto statistics = ParseStatistics();
  statistics.pairs.total_count = 1000;
  statistics.pairs.gap_fraction_too_high = 501;

  const auto messages = statistics_to_messages(statistics);
  ASSERT_EQ(messages.size(), 1);
  EXPECT_EQ(messages.front(),
            "50.1% of the read pairs were discarded because the gap fraction "
            "is too high. The flag \"--max_gap_fraction\" can be used to "
            "configure the maximum gap fraction.");
}

TEST_F(statistics_to_messages_Test,
       returns_correct_message_for_too_low_align_score_fraction) {
  auto statistics = ParseStatistics();
  statistics.pairs.total_count = 1000;
  statistics.pairs.align_score_fraction_too_low = 501;

  const auto messages = statistics_to_messages(statistics);
  ASSERT_EQ(messages.size(), 1);
  EXPECT_EQ(messages.front(),
            "50.1\% of the read pairs were discarded because the alignment "
            "score fraction (scaled value of TAG \"AS\") is too low. The flag "
            "\"--min_align_score_fraction\" can be used to configure the "
            "minimum alignment score fraction.");
}

TEST_F(statistics_to_messages_Test, reports_multiple_error_messages) {
  auto a = ParseStatistics();
  a.reads.total_count = 2000;
  a.reads.unmapped = 100;

  auto b = ParseStatistics();
  b.pairs.total_count = 1000;
  b.pairs.align_score_fraction_too_low = 501;

  auto combined = ParseStatistics();
  combined.reads.total_count = 2000;
  combined.reads.unmapped = 100;
  combined.pairs.total_count = 1000;
  combined.pairs.align_score_fraction_too_low = 501;

  const auto messages_a = statistics_to_messages(a);
  const auto messages_b = statistics_to_messages(b);
  const auto messages_combined = statistics_to_messages(combined);

  ASSERT_EQ(messages_a.size(), 1);
  ASSERT_EQ(messages_b.size(), 1);
  ASSERT_EQ(messages_combined.size(), 2);

  EXPECT_EQ(messages_combined.front(), messages_a.front());
  EXPECT_EQ(messages_combined.back(), messages_b.front());
}

} // namespace
} // namespace phaplo
