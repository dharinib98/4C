// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_utils_parameter_list.hpp"

#include "4C_utils_string.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Tuple.hpp>

#include <array>

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  void bool_parameter(std::string const& paramName, std::string const& value,
      std::string const& docString, SectionSpecs& section_specs)
  {
    // Check that value in lower-case is either "true" or "false", case insensitive
    const std::string value_lower = Core::Utils::to_lower(value);
    FOUR_C_ASSERT_ALWAYS(value_lower == "yes" || value_lower == "no",
        "Invalid value for boolean parameter '" + paramName + "': " + value);
    bool sane_default = value_lower == "yes";

    section_specs.specs.emplace_back(Core::IO::InputSpecBuilders::entry<bool>(
        paramName, {.description = docString, .default_value = sane_default}));
  }


  void int_parameter(std::string const& paramName, int const value, std::string const& docString,
      SectionSpecs& section_specs)
  {
    section_specs.specs.emplace_back(Core::IO::InputSpecBuilders::entry<int>(
        paramName, {.description = docString, .default_value = value}));
  }


  void double_parameter(std::string const& paramName, double const& value,
      std::string const& docString, SectionSpecs& section_specs)
  {
    section_specs.specs.emplace_back(Core::IO::InputSpecBuilders::entry<double>(
        paramName, {.description = docString, .default_value = value}));
  }


  void string_parameter(std::string const& paramName, std::string const& value,
      std::string const& docString, SectionSpecs& section_specs,
      std::vector<std::string> const& validParams)
  {
    if (validParams.empty())
    {
      section_specs.specs.emplace_back(Core::IO::InputSpecBuilders::entry<std::string>(
          paramName, {.description = docString, .default_value = value}));
    }
    else
    {
      section_specs.specs.emplace_back(Core::IO::InputSpecBuilders::selection<std::string>(
          paramName, validParams, {.description = docString, .default_value = value}));
    }
  }

}  // namespace Core::Utils

FOUR_C_NAMESPACE_CLOSE
