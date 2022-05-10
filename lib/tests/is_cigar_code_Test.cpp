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
#include <phaplo/is_cigar_code.hpp>

namespace phaplo {
namespace {
struct is_cigar_code_Test : ::testing::Test {};

TEST_F(is_cigar_code_Test, M_is_cigar_code) { EXPECT_TRUE(is_cigar_code('M')); }

TEST_F(is_cigar_code_Test, I_is_cigar_code) { EXPECT_TRUE(is_cigar_code('I')); }

TEST_F(is_cigar_code_Test, D_is_cigar_code) { EXPECT_TRUE(is_cigar_code('D')); }

TEST_F(is_cigar_code_Test, N_is_cigar_code) { EXPECT_TRUE(is_cigar_code('N')); }

TEST_F(is_cigar_code_Test, S_is_cigar_code) { EXPECT_TRUE(is_cigar_code('S')); }

TEST_F(is_cigar_code_Test, H_is_cigar_code) { EXPECT_TRUE(is_cigar_code('H')); }

TEST_F(is_cigar_code_Test, P_is_cigar_code) { EXPECT_TRUE(is_cigar_code('P')); }

TEST_F(is_cigar_code_Test, X_is_cigar_code) { EXPECT_TRUE(is_cigar_code('X')); }

TEST_F(is_cigar_code_Test, equal_sign_is_cigar_code) {
  EXPECT_TRUE(is_cigar_code('='));
}

TEST_F(is_cigar_code_Test, m_is_not_cigar_code) {
  EXPECT_FALSE(is_cigar_code('m'));
}

TEST_F(is_cigar_code_Test, 0_is_not_cigar_code) {
  EXPECT_FALSE(is_cigar_code('0'));
}

TEST_F(is_cigar_code_Test, 1_is_not_cigar_code) {
  EXPECT_FALSE(is_cigar_code('1'));
}

TEST_F(is_cigar_code_Test, 2_is_not_cigar_code) {
  EXPECT_FALSE(is_cigar_code('2'));
}

TEST_F(is_cigar_code_Test, 3_is_not_cigar_code) {
  EXPECT_FALSE(is_cigar_code('3'));
}

TEST_F(is_cigar_code_Test, 4_is_not_cigar_code) {
  EXPECT_FALSE(is_cigar_code('4'));
}

TEST_F(is_cigar_code_Test, 5_is_not_cigar_code) {
  EXPECT_FALSE(is_cigar_code('5'));
}

TEST_F(is_cigar_code_Test, 6_is_not_cigar_code) {
  EXPECT_FALSE(is_cigar_code('6'));
}

TEST_F(is_cigar_code_Test, 7_is_not_cigar_code) {
  EXPECT_FALSE(is_cigar_code('7'));
}

TEST_F(is_cigar_code_Test, 8_is_not_cigar_code) {
  EXPECT_FALSE(is_cigar_code('8'));
}

TEST_F(is_cigar_code_Test, 9_is_not_cigar_code) {
  EXPECT_FALSE(is_cigar_code('9'));
}

} // namespace
} // namespace phaplo