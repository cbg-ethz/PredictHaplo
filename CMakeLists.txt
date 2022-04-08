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

project(predicthaplo VERSION 2.1.4)
enable_testing()

find_package(LAPACK REQUIRED)

add_library(scythe INTERFACE)
target_include_directories(scythe
  INTERFACE
  SCYTHE/include
)
target_link_libraries(scythe
  INTERFACE
  ${LAPACK_LIBRARIES}
)
target_compile_definitions(scythe
  INTERFACE
  SCYTHE_LAPACK
  SCYTHE_DEBUG=0
  HAVE_TRUNC
)
target_compile_features(scythe
  INTERFACE
  cxx_std_11
)

add_subdirectory(lib)

add_executable(predicthaplo
  PredictHaplo_externAlign.cpp
)
target_link_libraries(predicthaplo
  PRIVATE
  scythe
  phaplo
)
target_compile_features(predicthaplo
  PRIVATE
  cxx_std_17
)
