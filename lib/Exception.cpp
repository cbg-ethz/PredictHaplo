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

#include <phaplo/Exception.hpp>
#include <type_traits>

namespace phaplo {

namespace {
std::string message(const ErrorCode code) noexcept {
  switch (code) {
  case ErrorCode::no_sam_file:
    return "Please provide the path to the SAM file via \"--sam\".";
  case ErrorCode::multiple_sam_files:
    return "Please provide only a single SAM file via \"--sam\".";
  case ErrorCode::no_valid_reads:
    return "No valid reads were discovered.";
  case ErrorCode::unknown:
    return "An unspecified error occured.";
  }

  return message(ErrorCode::unknown);
}
} // namespace

Exception::Exception(const std::string &message)
    : std::runtime_error(message) {}

Error::Error(const ErrorCode code)
    : Exception(message(code)),
      id_{static_cast<std::underlying_type_t<ErrorCode>>(code)} {}

int Error::id() const noexcept { return id_; }

} // namespace phaplo