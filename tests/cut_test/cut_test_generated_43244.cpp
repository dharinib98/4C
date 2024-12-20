// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_combintersection.hpp"
#include "4C_cut_levelsetintersection.hpp"
#include "4C_cut_meshintersection.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_tetmeshintersection.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.hpp"

void test_generated_43244()
{
  Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  std::vector<double> lsvs(8);
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.00666667;
    nids.push_back(35242);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35243);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.008;
    nids.push_back(-6349);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35243);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35259);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.008;
    nids.push_back(-6349);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.00666667;
    nids.push_back(35258);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.00666667;
    nids.push_back(35242);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.008;
    nids.push_back(-6349);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35243);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35244);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6350);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35244);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35260);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6350);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35260);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35259);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6350);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35259);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35243);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6350);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35244);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35245);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6351);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35245);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35261);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6351);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35261);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35260);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6351);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35260);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35244);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6351);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35245);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35246);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6352);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35246);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35262);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6352);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35262);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35261);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6352);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35261);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35245);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6352);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35246);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35247);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6353);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35247);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35263);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6353);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35263);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35262);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6353);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35262);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35246);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6353);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35259);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35260);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6365);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35260);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35276);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6365);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35276);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35275);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6365);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35260);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35261);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6366);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35261);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35277);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6366);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35277);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35276);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6366);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35276);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35260);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6366);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35261);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35262);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6367);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35262);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35278);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6367);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35278);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35277);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6367);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35277);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35261);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6367);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35262);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35263);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6368);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35263);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35279);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6368);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35279);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35278);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6368);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35278);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35262);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6368);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35275);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35276);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6380);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35276);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.218391;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35292);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6380);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35276);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35277);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6381);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35277);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.218391;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35293);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6381);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.218391;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35292);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35276);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6381);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35277);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35278);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6382);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35278);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.218391;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35294);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6382);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.218391;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35293);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35277);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6382);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35278);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35279);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6383);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35279);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.218391;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35295);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6383);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.218391;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35294);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35278);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6383);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.945652;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(31471);
    tri3_xyze(0, 1) = 0.956522;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32223);
    tri3_xyze(0, 2) = 0.951087;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5221);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.956522;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32223);
    tri3_xyze(0, 1) = 0.956522;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32207);
    tri3_xyze(0, 2) = 0.951087;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5221);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.945652;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(31487);
    tri3_xyze(0, 1) = 0.956522;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32239);
    tri3_xyze(0, 2) = 0.951087;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5223);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.956522;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32239);
    tri3_xyze(0, 1) = 0.956522;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32223);
    tri3_xyze(0, 2) = 0.951087;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5223);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.956522;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32223);
    tri3_xyze(0, 1) = 0.945652;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(31471);
    tri3_xyze(0, 2) = 0.951087;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5223);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.956522;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32255);
    tri3_xyze(0, 1) = 0.956522;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32239);
    tri3_xyze(0, 2) = 0.951087;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5225);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.956522;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32239);
    tri3_xyze(0, 1) = 0.945652;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(31487);
    tri3_xyze(0, 2) = 0.951087;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5225);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.142304;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35179);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.142304;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35180);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6290);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.142304;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35180);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35196);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6290);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35196);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35195);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6290);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.142304;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35180);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.142304;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35181);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6291);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.142304;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35181);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35197);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6291);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35197);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35196);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6291);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35196);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.142304;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35180);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6291);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.142304;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35181);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.142304;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35182);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6292);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.142304;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35182);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35198);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6292);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35198);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35197);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6292);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35197);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.142304;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35181);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6292);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.142304;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35182);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.142304;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35183);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6293);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.142304;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35183);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35199);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6293);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35199);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35198);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6293);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35198);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.142304;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35182);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6293);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35195);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35196);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6305);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35196);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35212);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6305);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35212);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35211);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6305);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35196);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35197);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6306);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35197);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35213);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6306);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35213);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35212);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6306);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35212);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35196);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6306);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35197);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35198);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6307);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35198);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35214);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6307);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35214);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35213);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6307);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35213);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35197);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6307);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35198);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35199);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6308);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35199);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35215);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6308);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35215);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35214);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6308);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35214);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35198);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6308);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.956522;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32223);
    tri3_xyze(0, 1) = 0.956522;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32239);
    tri3_xyze(0, 2) = 0.961957;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5345);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.956522;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32239);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32991);
    tri3_xyze(0, 2) = 0.961957;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5345);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32991);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32975);
    tri3_xyze(0, 2) = 0.961957;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5345);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32975);
    tri3_xyze(0, 1) = 0.956522;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32223);
    tri3_xyze(0, 2) = 0.961957;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5345);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.956522;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32239);
    tri3_xyze(0, 1) = 0.956522;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32255);
    tri3_xyze(0, 2) = 0.961957;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5347);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33007);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32991);
    tri3_xyze(0, 2) = 0.961957;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5347);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32991);
    tri3_xyze(0, 1) = 0.956522;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32239);
    tri3_xyze(0, 2) = 0.961957;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5347);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33759);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33775);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5593);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33775);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34527);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5593);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34527);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34511);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5593);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34511);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33759);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5593);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33775);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.218391;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33791);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5595);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.218391;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34543);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34527);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5595);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34527);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33775);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5595);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.956522;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32207);
    tri3_xyze(0, 1) = 0.956522;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32223);
    tri3_xyze(0, 2) = 0.961957;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5343);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.956522;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32223);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32975);
    tri3_xyze(0, 2) = 0.961957;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5343);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32975);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32959);
    tri3_xyze(0, 2) = 0.961957;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5343);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32959);
    tri3_xyze(0, 1) = 0.956522;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32207);
    tri3_xyze(0, 2) = 0.961957;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5343);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = 0.012;
    nids.push_back(35203);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = 0.00933333;
    nids.push_back(35204);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.0106667;
    nids.push_back(-6312);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = 0.00933333;
    nids.push_back(35204);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.00933333;
    nids.push_back(35220);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.0106667;
    nids.push_back(-6312);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.00933333;
    nids.push_back(35220);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.012;
    nids.push_back(35219);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.0106667;
    nids.push_back(-6312);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.012;
    nids.push_back(35219);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = 0.012;
    nids.push_back(35203);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.0106667;
    nids.push_back(-6312);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = 0.00933333;
    nids.push_back(35204);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = 0.00666667;
    nids.push_back(35205);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.008;
    nids.push_back(-6313);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = 0.00666667;
    nids.push_back(35205);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.00666667;
    nids.push_back(35221);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.008;
    nids.push_back(-6313);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.00666667;
    nids.push_back(35221);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.00933333;
    nids.push_back(35220);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.008;
    nids.push_back(-6313);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.00933333;
    nids.push_back(35220);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = 0.00933333;
    nids.push_back(35204);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.008;
    nids.push_back(-6313);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = 0.00666667;
    nids.push_back(35205);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = 0.004;
    nids.push_back(35206);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.00533333;
    nids.push_back(-6314);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = 0.004;
    nids.push_back(35206);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.004;
    nids.push_back(35222);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.00533333;
    nids.push_back(-6314);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.004;
    nids.push_back(35222);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.00666667;
    nids.push_back(35221);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.00533333;
    nids.push_back(-6314);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.00666667;
    nids.push_back(35221);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = 0.00666667;
    nids.push_back(35205);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.00533333;
    nids.push_back(-6314);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = 0.004;
    nids.push_back(35206);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = 0.00133333;
    nids.push_back(35207);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.00266667;
    nids.push_back(-6315);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = 0.00133333;
    nids.push_back(35207);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.00133333;
    nids.push_back(35223);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.00266667;
    nids.push_back(-6315);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.00133333;
    nids.push_back(35223);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.004;
    nids.push_back(35222);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.00266667;
    nids.push_back(-6315);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.004;
    nids.push_back(35222);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = 0.004;
    nids.push_back(35206);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 0.00266667;
    nids.push_back(-6315);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = 0.00133333;
    nids.push_back(35207);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.00133333;
    nids.push_back(35208);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 3.57165e-14;
    nids.push_back(-6316);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.00133333;
    nids.push_back(35208);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.00133333;
    nids.push_back(35224);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 3.57165e-14;
    nids.push_back(-6316);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.00133333;
    nids.push_back(35224);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.00133333;
    nids.push_back(35223);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 3.57165e-14;
    nids.push_back(-6316);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.00133333;
    nids.push_back(35223);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = 0.00133333;
    nids.push_back(35207);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = 3.57165e-14;
    nids.push_back(-6316);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.00133333;
    nids.push_back(35208);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.004;
    nids.push_back(35209);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.00266667;
    nids.push_back(-6317);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.004;
    nids.push_back(35209);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.004;
    nids.push_back(35225);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.00266667;
    nids.push_back(-6317);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.004;
    nids.push_back(35225);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.00133333;
    nids.push_back(35224);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.00266667;
    nids.push_back(-6317);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.00133333;
    nids.push_back(35224);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.00133333;
    nids.push_back(35208);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.00266667;
    nids.push_back(-6317);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.004;
    nids.push_back(35209);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.00666667;
    nids.push_back(35210);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.00533333;
    nids.push_back(-6318);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.00666667;
    nids.push_back(35210);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.00666667;
    nids.push_back(35226);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.00533333;
    nids.push_back(-6318);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.00666667;
    nids.push_back(35226);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.004;
    nids.push_back(35225);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.00533333;
    nids.push_back(-6318);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.004;
    nids.push_back(35225);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.004;
    nids.push_back(35209);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.00533333;
    nids.push_back(-6318);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.00666667;
    nids.push_back(35210);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35211);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.008;
    nids.push_back(-6319);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35211);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35227);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.008;
    nids.push_back(-6319);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35227);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.00666667;
    nids.push_back(35226);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.008;
    nids.push_back(-6319);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.00666667;
    nids.push_back(35226);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.00666667;
    nids.push_back(35210);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.008;
    nids.push_back(-6319);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35211);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35212);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6320);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35212);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35228);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6320);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35228);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35227);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6320);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35227);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35211);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6320);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35212);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35213);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6321);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35213);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35229);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6321);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35229);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35228);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6321);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35228);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35212);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6321);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35213);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35214);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6322);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35214);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35230);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6322);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35230);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35229);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6322);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35229);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35213);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6322);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35214);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35215);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6323);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35215);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35231);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6323);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35231);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35230);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6323);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35230);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35214);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6323);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.012;
    nids.push_back(35219);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.00933333;
    nids.push_back(35220);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.0106667;
    nids.push_back(-6327);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.00933333;
    nids.push_back(35220);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.00933333;
    nids.push_back(35236);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.0106667;
    nids.push_back(-6327);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.00933333;
    nids.push_back(35236);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.012;
    nids.push_back(35235);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.0106667;
    nids.push_back(-6327);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.012;
    nids.push_back(35235);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.012;
    nids.push_back(35219);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.0106667;
    nids.push_back(-6327);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.00933333;
    nids.push_back(35220);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.00666667;
    nids.push_back(35221);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.008;
    nids.push_back(-6328);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.00666667;
    nids.push_back(35221);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.00666667;
    nids.push_back(35237);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.008;
    nids.push_back(-6328);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.00666667;
    nids.push_back(35237);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.00933333;
    nids.push_back(35236);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.008;
    nids.push_back(-6328);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.00933333;
    nids.push_back(35236);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.00933333;
    nids.push_back(35220);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.008;
    nids.push_back(-6328);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.00666667;
    nids.push_back(35221);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.004;
    nids.push_back(35222);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.00533333;
    nids.push_back(-6329);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.004;
    nids.push_back(35222);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.004;
    nids.push_back(35238);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.00533333;
    nids.push_back(-6329);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.004;
    nids.push_back(35238);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.00666667;
    nids.push_back(35237);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.00533333;
    nids.push_back(-6329);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.00666667;
    nids.push_back(35237);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.00666667;
    nids.push_back(35221);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.00533333;
    nids.push_back(-6329);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.004;
    nids.push_back(35222);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.00133333;
    nids.push_back(35223);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.00266667;
    nids.push_back(-6330);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.00133333;
    nids.push_back(35223);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.00133333;
    nids.push_back(35239);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.00266667;
    nids.push_back(-6330);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.00133333;
    nids.push_back(35239);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.004;
    nids.push_back(35238);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.00266667;
    nids.push_back(-6330);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.004;
    nids.push_back(35238);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.004;
    nids.push_back(35222);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 0.00266667;
    nids.push_back(-6330);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = 0.00133333;
    nids.push_back(35223);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.00133333;
    nids.push_back(35224);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 3.66125e-14;
    nids.push_back(-6331);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.00133333;
    nids.push_back(35224);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.00133333;
    nids.push_back(35240);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 3.66125e-14;
    nids.push_back(-6331);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.00133333;
    nids.push_back(35240);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.00133333;
    nids.push_back(35239);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 3.66125e-14;
    nids.push_back(-6331);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.00133333;
    nids.push_back(35239);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = 0.00133333;
    nids.push_back(35223);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = 3.66125e-14;
    nids.push_back(-6331);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.00133333;
    nids.push_back(35224);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.004;
    nids.push_back(35225);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.00266667;
    nids.push_back(-6332);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.004;
    nids.push_back(35225);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.004;
    nids.push_back(35241);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.00266667;
    nids.push_back(-6332);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.004;
    nids.push_back(35241);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.00133333;
    nids.push_back(35240);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.00266667;
    nids.push_back(-6332);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.00133333;
    nids.push_back(35240);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.00133333;
    nids.push_back(35224);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.00266667;
    nids.push_back(-6332);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.004;
    nids.push_back(35225);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.00666667;
    nids.push_back(35226);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.00533333;
    nids.push_back(-6333);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.00666667;
    nids.push_back(35226);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.00666667;
    nids.push_back(35242);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.00533333;
    nids.push_back(-6333);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.00666667;
    nids.push_back(35242);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.004;
    nids.push_back(35241);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.00533333;
    nids.push_back(-6333);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.004;
    nids.push_back(35241);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.004;
    nids.push_back(35225);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.00533333;
    nids.push_back(-6333);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.00666667;
    nids.push_back(35226);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35227);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.008;
    nids.push_back(-6334);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35227);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35243);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.008;
    nids.push_back(-6334);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35243);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.00666667;
    nids.push_back(35242);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.008;
    nids.push_back(-6334);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.00666667;
    nids.push_back(35242);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.00666667;
    nids.push_back(35226);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.008;
    nids.push_back(-6334);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35227);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35228);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6335);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35228);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35244);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6335);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35244);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35243);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6335);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.00933333;
    nids.push_back(35243);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.00933333;
    nids.push_back(35227);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.0106667;
    nids.push_back(-6335);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35228);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35229);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6336);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35229);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35245);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6336);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35245);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35244);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6336);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.012;
    nids.push_back(35244);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.012;
    nids.push_back(35228);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.0133333;
    nids.push_back(-6336);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35229);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35230);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6337);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35230);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35246);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6337);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35246);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35245);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6337);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.0146667;
    nids.push_back(35245);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.0146667;
    nids.push_back(35229);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.016;
    nids.push_back(-6337);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35230);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35231);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6338);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35231);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35247);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6338);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35247);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35246);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6338);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.0173333;
    nids.push_back(35246);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.0173333;
    nids.push_back(35230);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.0186667;
    nids.push_back(-6338);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.012;
    nids.push_back(35235);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.00933333;
    nids.push_back(35236);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 0.0106667;
    nids.push_back(-6342);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.00933333;
    nids.push_back(35236);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = 0.00933333;
    nids.push_back(35252);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 0.0106667;
    nids.push_back(-6342);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.00933333;
    nids.push_back(35236);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.00666667;
    nids.push_back(35237);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 0.008;
    nids.push_back(-6343);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.00666667;
    nids.push_back(35237);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = 0.00666667;
    nids.push_back(35253);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 0.008;
    nids.push_back(-6343);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = 0.00933333;
    nids.push_back(35252);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.00933333;
    nids.push_back(35236);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 0.008;
    nids.push_back(-6343);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.00666667;
    nids.push_back(35237);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.004;
    nids.push_back(35238);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 0.00533333;
    nids.push_back(-6344);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.004;
    nids.push_back(35238);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = 0.004;
    nids.push_back(35254);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 0.00533333;
    nids.push_back(-6344);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = 0.00666667;
    nids.push_back(35253);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.00666667;
    nids.push_back(35237);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 0.00533333;
    nids.push_back(-6344);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.004;
    nids.push_back(35238);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.00133333;
    nids.push_back(35239);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 0.00266667;
    nids.push_back(-6345);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.00133333;
    nids.push_back(35239);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = 0.00133333;
    nids.push_back(35255);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 0.00266667;
    nids.push_back(-6345);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = 0.004;
    nids.push_back(35254);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.004;
    nids.push_back(35238);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 0.00266667;
    nids.push_back(-6345);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = 0.00133333;
    nids.push_back(35239);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.00133333;
    nids.push_back(35240);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 3.53636e-14;
    nids.push_back(-6346);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.00133333;
    nids.push_back(35240);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.00133333;
    nids.push_back(35256);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 3.53636e-14;
    nids.push_back(-6346);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = 0.00133333;
    nids.push_back(35255);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = 0.00133333;
    nids.push_back(35239);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = 3.53636e-14;
    nids.push_back(-6346);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.00133333;
    nids.push_back(35240);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.004;
    nids.push_back(35241);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.00266667;
    nids.push_back(-6347);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.004;
    nids.push_back(35241);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.004;
    nids.push_back(35257);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.00266667;
    nids.push_back(-6347);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.00133333;
    nids.push_back(35256);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.00133333;
    nids.push_back(35240);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.00266667;
    nids.push_back(-6347);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.004;
    nids.push_back(35241);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.00666667;
    nids.push_back(35242);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.00533333;
    nids.push_back(-6348);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.00666667;
    nids.push_back(35242);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.00666667;
    nids.push_back(35258);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.00533333;
    nids.push_back(-6348);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.004;
    nids.push_back(35257);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.004;
    nids.push_back(35241);
    tri3_xyze(0, 2) = 1;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.00533333;
    nids.push_back(-6348);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34511);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34527);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5715);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34527);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35279);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5715);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33743);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33727);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5467);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33727);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32975);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5467);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32991);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33007);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5469);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33007);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33759);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5469);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33759);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33743);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5469);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33743);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32991);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5469);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33023);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33775);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5471);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33775);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33759);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5471);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33759);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33007);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5471);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32943);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33695);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5461);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33695);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.142304;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33679);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5461);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32959);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33711);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5463);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33711);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33695);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5463);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33695);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32943);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5463);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32959);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32975);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5465);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32975);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33727);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5465);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33727);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33711);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5465);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33711);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32959);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5465);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32975);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(32991);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5467);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.967391;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(32991);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33743);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5467);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.218391;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33791);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33775);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5473);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33775);
    tri3_xyze(0, 1) = 0.967391;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33023);
    tri3_xyze(0, 2) = 0.972826;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5473);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.142304;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33679);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33695);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5583);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33695);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34447);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5583);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34447);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.142304;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34431);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5583);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.142304;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34431);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.142304;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33679);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5583);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33695);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33711);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5585);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33711);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34463);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5585);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34463);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34447);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5585);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34447);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33695);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5585);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33711);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33727);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5587);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33727);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34479);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5587);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34479);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34463);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5587);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34463);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33711);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5587);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33727);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33743);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5589);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33743);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34495);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5589);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34495);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34479);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5589);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34479);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33727);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5589);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33743);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33759);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5591);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.978261;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(33759);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34511);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5591);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34511);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34495);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5591);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34495);
    tri3_xyze(0, 1) = 0.978261;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(33743);
    tri3_xyze(0, 2) = 0.983696;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5591);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.142304;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34431);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34447);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5705);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34447);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35199);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5705);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35199);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.142304;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35183);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5705);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.142304;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35183);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.142304;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34431);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.147739;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5705);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34447);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34463);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5707);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34463);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35215);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5707);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35215);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35199);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5707);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.153174;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35199);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.153174;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34447);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.158609;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5707);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34463);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34479);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5709);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34479);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35231);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5709);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35231);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35215);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5709);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.164043;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35215);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.164043;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34463);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.169478;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5709);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34479);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34495);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5711);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34495);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35247);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5711);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35247);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35231);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5711);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.174913;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35231);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.174913;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34479);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.180348;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5711);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34495);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34511);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5713);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34511);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35263);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5713);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35263);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35247);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5713);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.185783;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35247);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.185783;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34495);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.191217;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5713);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35279);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35263);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5715);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.196652;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35263);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.196652;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34511);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.202087;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5715);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.98913;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(34527);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.218391;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34543);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5717);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.218391;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35295);
    tri3_xyze(0, 1) = 1;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(35279);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5717);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1;
    tri3_xyze(1, 0) = 0.207522;
    tri3_xyze(2, 0) = -0.02;
    nids.push_back(35279);
    tri3_xyze(0, 1) = 0.98913;
    tri3_xyze(1, 1) = 0.207522;
    tri3_xyze(2, 1) = -0.02;
    nids.push_back(34527);
    tri3_xyze(0, 2) = 0.994565;
    tri3_xyze(1, 2) = 0.212957;
    tri3_xyze(2, 2) = -0.02;
    nids.push_back(-5717);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.977778;
    hex8_xyze(1, 0) = 0.166667;
    hex8_xyze(2, 0) = -0.0543478;
    nids.push_back(376614);
    hex8_xyze(0, 1) = 1;
    hex8_xyze(1, 1) = 0.166667;
    hex8_xyze(2, 1) = -0.0543478;
    nids.push_back(376616);
    hex8_xyze(0, 2) = 1;
    hex8_xyze(1, 2) = 0.188889;
    hex8_xyze(2, 2) = -0.0543478;
    nids.push_back(376978);
    hex8_xyze(0, 3) = 0.977778;
    hex8_xyze(1, 3) = 0.188889;
    hex8_xyze(2, 3) = -0.0543478;
    nids.push_back(376976);
    hex8_xyze(0, 4) = 0.977778;
    hex8_xyze(1, 4) = 0.166667;
    hex8_xyze(2, 4) = -0.0326087;
    nids.push_back(409556);
    hex8_xyze(0, 5) = 1;
    hex8_xyze(1, 5) = 0.166667;
    hex8_xyze(2, 5) = -0.0326087;
    nids.push_back(409558);
    hex8_xyze(0, 6) = 1;
    hex8_xyze(1, 6) = 0.188889;
    hex8_xyze(2, 6) = -0.0326087;
    nids.push_back(409920);
    hex8_xyze(0, 7) = 0.977778;
    hex8_xyze(1, 7) = 0.188889;
    hex8_xyze(2, 7) = -0.0326087;
    nids.push_back(409918);

    intersection.add_element(39194, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.977778;
    hex8_xyze(1, 0) = 0.144444;
    hex8_xyze(2, 0) = -0.0326087;
    nids.push_back(409194);
    hex8_xyze(0, 1) = 1;
    hex8_xyze(1, 1) = 0.144444;
    hex8_xyze(2, 1) = -0.0326087;
    nids.push_back(409196);
    hex8_xyze(0, 2) = 1;
    hex8_xyze(1, 2) = 0.166667;
    hex8_xyze(2, 2) = -0.0326087;
    nids.push_back(409558);
    hex8_xyze(0, 3) = 0.977778;
    hex8_xyze(1, 3) = 0.166667;
    hex8_xyze(2, 3) = -0.0326087;
    nids.push_back(409556);
    hex8_xyze(0, 4) = 0.977778;
    hex8_xyze(1, 4) = 0.144444;
    hex8_xyze(2, 4) = -0.0108696;
    nids.push_back(442136);
    hex8_xyze(0, 5) = 1;
    hex8_xyze(1, 5) = 0.144444;
    hex8_xyze(2, 5) = -0.0108696;
    nids.push_back(442138);
    hex8_xyze(0, 6) = 1;
    hex8_xyze(1, 6) = 0.166667;
    hex8_xyze(2, 6) = -0.0108696;
    nids.push_back(442500);
    hex8_xyze(0, 7) = 0.977778;
    hex8_xyze(1, 7) = 0.166667;
    hex8_xyze(2, 7) = -0.0108696;
    nids.push_back(442498);

    intersection.add_element(43154, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.955556;
    hex8_xyze(1, 0) = 0.166667;
    hex8_xyze(2, 0) = -0.0326087;
    nids.push_back(409554);
    hex8_xyze(0, 1) = 0.977778;
    hex8_xyze(1, 1) = 0.166667;
    hex8_xyze(2, 1) = -0.0326087;
    nids.push_back(409556);
    hex8_xyze(0, 2) = 0.977778;
    hex8_xyze(1, 2) = 0.188889;
    hex8_xyze(2, 2) = -0.0326087;
    nids.push_back(409918);
    hex8_xyze(0, 3) = 0.955556;
    hex8_xyze(1, 3) = 0.188889;
    hex8_xyze(2, 3) = -0.0326087;
    nids.push_back(409916);
    hex8_xyze(0, 4) = 0.955556;
    hex8_xyze(1, 4) = 0.166667;
    hex8_xyze(2, 4) = -0.0108696;
    nids.push_back(442496);
    hex8_xyze(0, 5) = 0.977778;
    hex8_xyze(1, 5) = 0.166667;
    hex8_xyze(2, 5) = -0.0108696;
    nids.push_back(442498);
    hex8_xyze(0, 6) = 0.977778;
    hex8_xyze(1, 6) = 0.188889;
    hex8_xyze(2, 6) = -0.0108696;
    nids.push_back(442860);
    hex8_xyze(0, 7) = 0.955556;
    hex8_xyze(1, 7) = 0.188889;
    hex8_xyze(2, 7) = -0.0108696;
    nids.push_back(442858);

    intersection.add_element(43243, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.977778;
    hex8_xyze(1, 0) = 0.166667;
    hex8_xyze(2, 0) = -0.0326087;
    nids.push_back(409556);
    hex8_xyze(0, 1) = 1;
    hex8_xyze(1, 1) = 0.166667;
    hex8_xyze(2, 1) = -0.0326087;
    nids.push_back(409558);
    hex8_xyze(0, 2) = 1;
    hex8_xyze(1, 2) = 0.188889;
    hex8_xyze(2, 2) = -0.0326087;
    nids.push_back(409920);
    hex8_xyze(0, 3) = 0.977778;
    hex8_xyze(1, 3) = 0.188889;
    hex8_xyze(2, 3) = -0.0326087;
    nids.push_back(409918);
    hex8_xyze(0, 4) = 0.977778;
    hex8_xyze(1, 4) = 0.166667;
    hex8_xyze(2, 4) = -0.0108696;
    nids.push_back(442498);
    hex8_xyze(0, 5) = 1;
    hex8_xyze(1, 5) = 0.166667;
    hex8_xyze(2, 5) = -0.0108696;
    nids.push_back(442500);
    hex8_xyze(0, 6) = 1;
    hex8_xyze(1, 6) = 0.188889;
    hex8_xyze(2, 6) = -0.0108696;
    nids.push_back(442862);
    hex8_xyze(0, 7) = 0.977778;
    hex8_xyze(1, 7) = 0.188889;
    hex8_xyze(2, 7) = -0.0108696;
    nids.push_back(442860);

    intersection.add_element(43244, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 1;
    hex8_xyze(1, 0) = 0.166667;
    hex8_xyze(2, 0) = -0.0326087;
    nids.push_back(409558);
    hex8_xyze(0, 1) = 1.02222;
    hex8_xyze(1, 1) = 0.166667;
    hex8_xyze(2, 1) = -0.0326087;
    nids.push_back(409560);
    hex8_xyze(0, 2) = 1.02222;
    hex8_xyze(1, 2) = 0.188889;
    hex8_xyze(2, 2) = -0.0326087;
    nids.push_back(409922);
    hex8_xyze(0, 3) = 1;
    hex8_xyze(1, 3) = 0.188889;
    hex8_xyze(2, 3) = -0.0326087;
    nids.push_back(409920);
    hex8_xyze(0, 4) = 1;
    hex8_xyze(1, 4) = 0.166667;
    hex8_xyze(2, 4) = -0.0108696;
    nids.push_back(442500);
    hex8_xyze(0, 5) = 1.02222;
    hex8_xyze(1, 5) = 0.166667;
    hex8_xyze(2, 5) = -0.0108696;
    nids.push_back(442502);
    hex8_xyze(0, 6) = 1.02222;
    hex8_xyze(1, 6) = 0.188889;
    hex8_xyze(2, 6) = -0.0108696;
    nids.push_back(442864);
    hex8_xyze(0, 7) = 1;
    hex8_xyze(1, 7) = 0.188889;
    hex8_xyze(2, 7) = -0.0108696;
    nids.push_back(442862);

    intersection.add_element(43245, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.977778;
    hex8_xyze(1, 0) = 0.188889;
    hex8_xyze(2, 0) = -0.0326087;
    nids.push_back(409918);
    hex8_xyze(0, 1) = 1;
    hex8_xyze(1, 1) = 0.188889;
    hex8_xyze(2, 1) = -0.0326087;
    nids.push_back(409920);
    hex8_xyze(0, 2) = 1;
    hex8_xyze(1, 2) = 0.211111;
    hex8_xyze(2, 2) = -0.0326087;
    nids.push_back(410282);
    hex8_xyze(0, 3) = 0.977778;
    hex8_xyze(1, 3) = 0.211111;
    hex8_xyze(2, 3) = -0.0326087;
    nids.push_back(410280);
    hex8_xyze(0, 4) = 0.977778;
    hex8_xyze(1, 4) = 0.188889;
    hex8_xyze(2, 4) = -0.0108696;
    nids.push_back(442860);
    hex8_xyze(0, 5) = 1;
    hex8_xyze(1, 5) = 0.188889;
    hex8_xyze(2, 5) = -0.0108696;
    nids.push_back(442862);
    hex8_xyze(0, 6) = 1;
    hex8_xyze(1, 6) = 0.211111;
    hex8_xyze(2, 6) = -0.0108696;
    nids.push_back(443224);
    hex8_xyze(0, 7) = 0.977778;
    hex8_xyze(1, 7) = 0.211111;
    hex8_xyze(2, 7) = -0.0108696;
    nids.push_back(443222);

    intersection.add_element(43334, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.977778;
    hex8_xyze(1, 0) = 0.166667;
    hex8_xyze(2, 0) = -0.0108696;
    nids.push_back(442498);
    hex8_xyze(0, 1) = 1;
    hex8_xyze(1, 1) = 0.166667;
    hex8_xyze(2, 1) = -0.0108696;
    nids.push_back(442500);
    hex8_xyze(0, 2) = 1;
    hex8_xyze(1, 2) = 0.188889;
    hex8_xyze(2, 2) = -0.0108696;
    nids.push_back(442862);
    hex8_xyze(0, 3) = 0.977778;
    hex8_xyze(1, 3) = 0.188889;
    hex8_xyze(2, 3) = -0.0108696;
    nids.push_back(442860);
    hex8_xyze(0, 4) = 0.977778;
    hex8_xyze(1, 4) = 0.166667;
    hex8_xyze(2, 4) = 0.0108696;
    nids.push_back(475440);
    hex8_xyze(0, 5) = 1;
    hex8_xyze(1, 5) = 0.166667;
    hex8_xyze(2, 5) = 0.0108696;
    nids.push_back(475442);
    hex8_xyze(0, 6) = 1;
    hex8_xyze(1, 6) = 0.188889;
    hex8_xyze(2, 6) = 0.0108696;
    nids.push_back(475804);
    hex8_xyze(0, 7) = 0.977778;
    hex8_xyze(1, 7) = 0.188889;
    hex8_xyze(2, 7) = 0.0108696;
    nids.push_back(475802);

    intersection.add_element(47294, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  intersection.cut_test_cut(
      true, Cut::VCellGaussPts_DirectDivergence, Cut::BCellGaussPts_Tessellation);
  intersection.cut_finalize(
      true, Cut::VCellGaussPts_DirectDivergence, Cut::BCellGaussPts_Tessellation, false, true);
}
