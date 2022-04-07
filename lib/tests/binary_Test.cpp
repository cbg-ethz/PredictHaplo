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
struct binary_Test : ::testing::Test {
  std::stringstream stream;
};

TEST_F(binary_Test, binary_of_0_is_0) { EXPECT_EQ(binary(0, stream), "0"); }

TEST_F(binary_Test,
       binary_of_negative_number_is_string_representation_of_number) {
  EXPECT_EQ(binary(-1, stream), "-1");
}

TEST_F(binary_Test,
       binary_of_positive_number_is_binary_representation_of_number) {
  EXPECT_EQ(binary(6, stream), "110");
}

} // namespace
} // namespace phaplo