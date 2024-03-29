# Copyright (C) 2022 ETH Zurich
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

add_library(phaplo
  Exception.cpp
  ParseStatistics.cpp
  binary.cpp
  is_alignment_match.cpp
  is_cigar_code.cpp
  median.cpp
  to_vector_of_ints.cpp
  statistics_to_messages.cpp
)
target_include_directories(phaplo
  PUBLIC
  include
)
target_compile_features(phaplo
  PUBLIC
  cxx_std_17
)

add_subdirectory(tests)