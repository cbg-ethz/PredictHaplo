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
#include <phaplo/is_alignment_match.hpp>

namespace phaplo {
namespace {
struct is_alignment_match_Test : ::testing::Test {};

TEST_F(is_alignment_match_Test, M_represents_alignment_match) {
  EXPECT_TRUE(is_alignment_match('M'));
}

TEST_F(is_alignment_match_Test, X_represents_alignment_match) {
  EXPECT_TRUE(is_alignment_match('X'));
}

TEST_F(is_alignment_match_Test, equal_sign_represents_alignment_match) {
  EXPECT_TRUE(is_alignment_match('='));
}

TEST_F(is_alignment_match_Test, D_does_not_represent_alignment_match) {
  EXPECT_FALSE(is_alignment_match('D'));
}

TEST_F(is_alignment_match_Test, m_does_not_represent_alignment_match) {
  EXPECT_FALSE(is_alignment_match('m'));
}
} // namespace
} // namespace phaplo