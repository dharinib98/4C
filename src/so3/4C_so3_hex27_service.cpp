// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_so3_hex27.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  return Center Coords in Reference System                            |
 *----------------------------------------------------------------------*/
std::vector<double> Discret::Elements::SoHex27::soh27_element_center_refe_coords()
{
  // update element geometry
  Core::LinAlg::Matrix<NUMNOD_SOH27, NUMDIM_SOH27> xrefe;  // material coord. of element
  for (int i = 0; i < NUMNOD_SOH27; ++i)
  {
    const auto& x = nodes()[i]->x();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  const Core::FE::CellType distype = shape();
  Core::LinAlg::Matrix<NUMNOD_SOH27, 1> funct;
  // Element midpoint at r=s=t=0.0
  Core::FE::shape_function_3d(funct, 0.0, 0.0, 0.0, distype);
  Core::LinAlg::Matrix<1, NUMDIM_SOH27> midpoint;
  // midpoint.multiply('T','N',1.0,funct,xrefe,0.0);
  midpoint.multiply_tn(funct, xrefe);
  std::vector<double> centercoords(3);
  centercoords[0] = midpoint(0, 0);
  centercoords[1] = midpoint(0, 1);
  centercoords[2] = midpoint(0, 2);
  return centercoords;
}

FOUR_C_NAMESPACE_CLOSE
