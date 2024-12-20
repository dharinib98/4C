// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_largerotations.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
unsigned int Core::LargeRotations::numbering_trafo(const unsigned int j, const unsigned int numnode)
{
  // Node numbering j=1,...,num_node() according to Crisfield 1999:
  // LINE2  1---2
  // LINE3  1---2---3
  // LINE4  1---2---3---4
  // LINE5  1---2---3---4---5

  // Storage position i=0,...,num_node()-1 of nodal quantities applied in 4C:
  // LINE2  (1,2)
  // LINE3  (1,3,2)
  // LINE4  (1,4,2,3)
  // LINE5  (1,5,2,3,4)

  // Initialization
  unsigned int i = 0;

  // Transformation
  if (j == 1)
    i = 0;
  else if (j == numnode)
    i = 1;
  else
    i = j;

  return i;
}

FOUR_C_NAMESPACE_CLOSE
