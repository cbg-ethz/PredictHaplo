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
#include <phaplo/median.hpp>
#include <vector>

namespace phaplo {
namespace {
struct median_Test : ::testing::Test {};

TEST_F(median_Test, median_of_empty_vector_is_nan) {
  EXPECT_TRUE(std::isnan(median({})));
}

TEST_F(median_Test, median_of_single_element_vector_is_element) {
  EXPECT_DOUBLE_EQ(median({2}), 2.);
}

TEST_F(median_Test, median_of_two_element_vector_is_mean) {
  EXPECT_DOUBLE_EQ(median({2, 3}), 2.5);
}

TEST_F(median_Test, median_of_three_element_vector_is_center_element) {
  EXPECT_DOUBLE_EQ(median({2, 3, 4}), 3.);
}

TEST_F(median_Test, median_is_robust_to_outliers) {
  EXPECT_DOUBLE_EQ(median({2, 3, 99}), 3.);
}
} // namespace
} // namespace phaplo