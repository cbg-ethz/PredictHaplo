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

#include <gtest/gtest.h>
#include <phaplo/binary.hpp>
#include <vector>

namespace phaplo {
namespace {
struct binary_Test : ::testing::Test {};

TEST_F(binary_Test, binary_of_0_is_only_unset_bits) {
  EXPECT_FALSE(binary(0).any());
}

TEST_F(binary_Test,
       binary_of_positive_number_is_binary_representation_of_number) {
  const auto result = binary(6);

  ASSERT_GE(result.size(), 3);
  EXPECT_EQ(result[0], false);
  EXPECT_EQ(result[1], true);
  EXPECT_EQ(result[2], true);

  EXPECT_EQ(result, decltype(result)("110"));
}

TEST_F(binary_Test, used_bits_returns_0_for_0) {
  EXPECT_EQ(used_bits(binary(0)), 0);
}

TEST_F(binary_Test, used_bits_returns_1_for_1) {
  EXPECT_EQ(used_bits(binary(1)), 1);
}

TEST_F(binary_Test, used_bits_returns_3_for_6) {
  EXPECT_EQ(used_bits(binary(6)), 3);
}


} // namespace
} // namespace phaplo