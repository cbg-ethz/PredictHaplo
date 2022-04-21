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

#include <cmath>
#include <gtest/gtest.h>
#include <phaplo/Exception.hpp>

namespace phaplo {
namespace {
struct Exception_Test : ::testing::Test {};

TEST_F(Exception_Test, what_returns_message_passed_to_constructor) {
  const auto message = std::string("abc");
  EXPECT_STREQ(Exception(message).what(), message.c_str());
}

TEST_F(Exception_Test, unknown_error_has_id_1) {
  const auto error = Error(ErrorCode::unknown);
  EXPECT_EQ(error.id(), 1);
}

TEST_F(Exception_Test, unknown_error_has_message_informing_user) {
  const auto error = Error(ErrorCode::unknown);
  EXPECT_STREQ(error.what(), "An unspecified error occured.");
}

TEST_F(Exception_Test, no_sam_file_error_has_id_2) {
  const auto error = Error(ErrorCode::no_sam_file);
  EXPECT_EQ(error.id(), 2);
}

TEST_F(Exception_Test, no_sam_file_has_message_informing_user) {
  const auto error = Error(ErrorCode::no_sam_file);
  EXPECT_STREQ(error.what(),
               "Please provide the path to the SAM file via \"--sam\".");
}

TEST_F(Exception_Test, multiple_sam_files_error_has_id_3) {
  const auto error = Error(ErrorCode::multiple_sam_files);
  EXPECT_EQ(error.id(), 3);
}

TEST_F(Exception_Test, multiple_sam_files_has_message_informing_user) {
  const auto error = Error(ErrorCode::multiple_sam_files);
  EXPECT_STREQ(error.what(),
               "Please provide only a single SAM file via \"--sam\".");
}

TEST_F(Exception_Test, no_valid_reads_error_has_id_4) {
  const auto error = Error(ErrorCode::no_valid_reads);
  EXPECT_EQ(error.id(), 4);
}

TEST_F(Exception_Test, no_valid_reads_has_message_informing_user) {
  const auto error = Error(ErrorCode::no_valid_reads);
  EXPECT_STREQ(error.what(), "No valid reads were discovered.");
}

TEST_F(Exception_Test, parse_sam_failed_error_has_id_5) {
  const auto error = Error(ErrorCode::parse_sam_failed);
  EXPECT_EQ(error.id(), 5);
}

TEST_F(Exception_Test, parse_sam_failed_has_message_informing_user) {
  const auto error = Error(ErrorCode::parse_sam_failed);
  EXPECT_STREQ(error.what(), "Parsing the SAM file failed.");
}

TEST_F(Exception_Test, no_reference_file_error_has_id_6) {
  const auto error = Error(ErrorCode::no_reference_file);
  EXPECT_EQ(error.id(), 6);
}

TEST_F(Exception_Test, no_reference_file_has_message_informing_user) {
  const auto error = Error(ErrorCode::no_reference_file);
  EXPECT_STREQ(error.what(), "Please provide the path to the reference "
                             "sequence file via \"--reference\".");
}

} // namespace
} // namespace phaplo