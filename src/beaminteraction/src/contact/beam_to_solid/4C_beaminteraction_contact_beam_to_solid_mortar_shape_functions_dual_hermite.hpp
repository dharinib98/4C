// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_CONTACT_BEAM_TO_SOLID_MORTAR_SHAPE_FUNCTIONS_DUAL_HERMITE_HPP
#define FOUR_C_BEAMINTERACTION_CONTACT_BEAM_TO_SOLID_MORTAR_SHAPE_FUNCTIONS_DUAL_HERMITE_HPP

#include "4C_config.hpp"

#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_shape_functions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace BeamInteraction
{
  struct t_hermite_dual
      : public GeometryPair::ElementDiscretizationBase<Core::FE::CellType::line2, 2>
  {
  };
}  // namespace BeamInteraction

namespace GeometryPair
{
  template <>
  struct ShapeFunctionData<BeamInteraction::t_hermite_dual>
  {
    double ref_length_;
  };

  template <>
  struct EvaluateShapeFunction<BeamInteraction::t_hermite_dual>
  {
    template <typename V, typename T>
    static void evaluate(V& N, const T& xi,
        const ShapeFunctionData<BeamInteraction::t_hermite_dual>& shape_function_data)
    {
      switch (BeamInteraction::t_hermite_dual::discretization_)
      {
        case Core::FE::CellType::line2:
        {
          const T xi2 = xi * xi;
          const T xi3 = xi2 * xi;

          N(0) = -35.0 / 4.0 * xi3 + 15.0 / 4.0 * xi2 + 15.0 / 4.0 * xi - 3.0 / 4.0;

          N(1) = (35.0 / 2.0 * xi3 - 15.0 / 4.0 * xi2 - 10.0 * xi + 5.0 / 4.0) /
                 (shape_function_data.ref_length_);

          N(2) = 35.0 / 4.0 * xi3 + 15.0 / 4.0 * xi2 - 15.0 / 4.0 * xi - 3.0 / 4.0;

          N(3) = (-35.0 / 2.0 * xi3 - 15.0 / 4.0 * xi2 + 10.0 * xi + 5.0 / 4.0) /
                 (shape_function_data.ref_length_);
          break;
        }
        default:
          FOUR_C_THROW("distype unknown\n");
          break;
      }

      return;
    }
  };

  template <>
  struct PrintElementData<BeamInteraction::t_hermite_dual>
  {
    template <typename ScalarType>
    static void print(const ElementData<BeamInteraction::t_hermite_dual, ScalarType>& element_data,
        std::ostream& out)
    {
      constexpr auto max_precision{std::numeric_limits<double>::digits10 + 1};
      out << std::setprecision(max_precision);
      out << "\nElement reference length: " << element_data.shape_function_data_.ref_length_;
      out << "\nElement state vector: ";
      element_data.element_position_.print(out);
    }
  };
}  // namespace GeometryPair

namespace BeamInteraction
{
  template <typename Mortar, typename BeamElementData>
  inline void set_beam_to_solid_mortar_shape_function_data(
      GeometryPair::ElementData<Mortar, double>& element_data_lambda,
      const BeamElementData& beam_element_data)
  {
    // Default: no additional mortar shape-function data needed.
  }

  template <typename BeamElementData>
  inline void set_beam_to_solid_mortar_shape_function_data(
      GeometryPair::ElementData<t_hermite_dual, double>& element_data_lambda,
      const BeamElementData& beam_element_data)
  {
    element_data_lambda.shape_function_data_.ref_length_ =
        beam_element_data.shape_function_data_.ref_length_;
  }

  template <typename Mortar, typename V, typename T>
  inline void evaluate_beam_to_solid_mortar_shape_function(
      V& N, const T& eta, const GeometryPair::ElementData<Mortar, double>& mortar_element_data)
  {
    GeometryPair::EvaluateShapeFunction<Mortar>::evaluate(
        N, eta, mortar_element_data.shape_function_data_);
  }
}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE

#endif