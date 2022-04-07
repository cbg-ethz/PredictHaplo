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
#include <phaplo/to_vector_of_ints.hpp>

namespace phaplo {
namespace {
struct to_vector_of_ints_Test : ::testing::Test {};

TEST_F(to_vector_of_ints_Test, converts_empty_string_to_empty_vector) {
  EXPECT_EQ(to_vector_of_ints(""), std::vector<int>{});
}

TEST_F(to_vector_of_ints_Test, converts_1_to_vector_of_1) {
  EXPECT_EQ(to_vector_of_ints("1"), std::vector<int>{1});
}

TEST_F(to_vector_of_ints_Test, converts_123_to_vector_of_1_2_3) {
  const auto reference = std::vector<int>{1, 2, 3};
  EXPECT_EQ(to_vector_of_ints("123"), reference);
}

} // namespace
} // namespace phaplo