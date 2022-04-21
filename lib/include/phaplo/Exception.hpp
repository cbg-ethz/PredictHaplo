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

#pragma once

#include <stdexcept>

namespace phaplo {

/**
 * @brief Base class of all phaplo exceptions.
 */
struct Exception : std::runtime_error {
  explicit Exception(const std::string &message);
};

/**
 * @brief Error codes used by phaplo.
 */
enum class ErrorCode : int {
  unknown = 1,
  no_sam_file,
  multiple_sam_files,
  no_valid_reads,
  parse_sam_failed,
  no_reference_file,
  unsupported_flag,
};

/**
 * @brief Wraps an ErrorCode in an Exception.
 */
class Error : public Exception {
public:
  explicit Error(const ErrorCode code);
  int id() const noexcept;

private:
  int id_;
};
} // namespace phaplo
