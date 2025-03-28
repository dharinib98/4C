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

void test_generated_7022()
{
  Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  std::vector<double> lsvs(8);
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.10169279791156897452;
    tri3_xyze(1, 0) = 0.1016927979115599262;
    tri3_xyze(2, 0) = 0.82408768370493257205;
    nids.push_back(7849);
    tri3_xyze(0, 1) = 0.099313748259548931263;
    tri3_xyze(1, 1) = 0.09931374825953985519;
    tri3_xyze(2, 1) = 0.82938926261447065702;
    nids.push_back(7851);
    tri3_xyze(0, 2) = 0.1117970706775816836;
    tri3_xyze(1, 2) = 0.085784909507892970204;
    tri3_xyze(2, 2) = 0.82673847315970205862;
    nids.push_back(-1131);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.099313748259548931263;
    tri3_xyze(1, 0) = 0.09931374825953985519;
    tri3_xyze(2, 0) = 0.82938926261447065702;
    nids.push_back(7851);
    tri3_xyze(0, 1) = 0.12163400383955559614;
    tri3_xyze(1, 1) = 0.070225424859374732667;
    tri3_xyze(2, 1) = 0.8293892626144715452;
    nids.push_back(8101);
    tri3_xyze(0, 2) = 0.1117970706775816836;
    tri3_xyze(1, 2) = 0.085784909507892970204;
    tri3_xyze(2, 2) = 0.82673847315970205862;
    nids.push_back(-1131);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.099313748259548931263;
    tri3_xyze(1, 0) = 0.09931374825953985519;
    tri3_xyze(2, 0) = 0.82938926261447065702;
    nids.push_back(7851);
    tri3_xyze(0, 1) = 0.096483611104761266009;
    tri3_xyze(1, 1) = 0.096483611104752328713;
    tri3_xyze(2, 1) = 0.83422735529628144757;
    nids.push_back(7853);
    tri3_xyze(0, 2) = 0.10889979276944010411;
    tri3_xyze(1, 2) = 0.083561749977300828207;
    tri3_xyze(2, 2) = 0.83180830895537649639;
    nids.push_back(-1132);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.096483611104761266009;
    tri3_xyze(1, 0) = 0.096483611104752328713;
    tri3_xyze(2, 0) = 0.83422735529628144757;
    nids.push_back(7853);
    tri3_xyze(0, 1) = 0.11816780787389465079;
    tri3_xyze(1, 1) = 0.068224215685536396259;
    tri3_xyze(2, 1) = 0.83422735529628233575;
    nids.push_back(8103);
    tri3_xyze(0, 2) = 0.10889979276944010411;
    tri3_xyze(1, 2) = 0.083561749977300828207;
    tri3_xyze(2, 2) = 0.83180830895537649639;
    nids.push_back(-1132);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.12163400383955559614;
    tri3_xyze(1, 0) = 0.070225424859374732667;
    tri3_xyze(2, 0) = 0.8293892626144715452;
    nids.push_back(8101);
    tri3_xyze(0, 1) = 0.099313748259548931263;
    tri3_xyze(1, 1) = 0.09931374825953985519;
    tri3_xyze(2, 1) = 0.82938926261447065702;
    nids.push_back(7851);
    tri3_xyze(0, 2) = 0.10889979276944010411;
    tri3_xyze(1, 2) = 0.083561749977300828207;
    tri3_xyze(2, 2) = 0.83180830895537649639;
    nids.push_back(-1132);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.096483611104761266009;
    tri3_xyze(1, 0) = 0.096483611104752328713;
    tri3_xyze(2, 0) = 0.83422735529628144757;
    nids.push_back(7853);
    tri3_xyze(0, 1) = 0.093247019400778799469;
    tri3_xyze(1, 1) = 0.093247019400770181363;
    tri3_xyze(2, 1) = 0.83852566213863666267;
    nids.push_back(7855);
    tri3_xyze(0, 2) = 0.10552556179077209053;
    tri3_xyze(1, 2) = 0.080972611483694403867;
    tri3_xyze(2, 2) = 0.83637650871745949921;
    nids.push_back(-1133);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.093247019400778799469;
    tri3_xyze(1, 0) = 0.093247019400770181363;
    tri3_xyze(2, 0) = 0.83852566213863666267;
    nids.push_back(7855);
    tri3_xyze(0, 1) = 0.11420380878365363198;
    tri3_xyze(1, 1) = 0.065935599743718723009;
    tri3_xyze(2, 1) = 0.8385256621386373288;
    nids.push_back(8105);
    tri3_xyze(0, 2) = 0.10552556179077209053;
    tri3_xyze(1, 2) = 0.080972611483694403867;
    tri3_xyze(2, 2) = 0.83637650871745949921;
    nids.push_back(-1133);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.11816780787389465079;
    tri3_xyze(1, 0) = 0.068224215685536396259;
    tri3_xyze(2, 0) = 0.83422735529628233575;
    nids.push_back(8103);
    tri3_xyze(0, 1) = 0.096483611104761266009;
    tri3_xyze(1, 1) = 0.096483611104752328713;
    tri3_xyze(2, 1) = 0.83422735529628144757;
    nids.push_back(7853);
    tri3_xyze(0, 2) = 0.10552556179077209053;
    tri3_xyze(1, 2) = 0.080972611483694403867;
    tri3_xyze(2, 2) = 0.83637650871745949921;
    nids.push_back(-1133);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.093247019400778799469;
    tri3_xyze(1, 0) = 0.093247019400770181363;
    tri3_xyze(2, 0) = 0.83852566213863666267;
    nids.push_back(7855);
    tri3_xyze(0, 1) = 0.0896550161322196143;
    tri3_xyze(1, 1) = 0.08965501613221139865;
    tri3_xyze(2, 1) = 0.84221639627494804259;
    nids.push_back(7857);
    tri3_xyze(0, 2) = 0.10172759137977901989;
    tri3_xyze(1, 2) = 0.078058326287794216336;
    tri3_xyze(2, 2) = 0.84037102920679274121;
    nids.push_back(-1134);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0896550161322196143;
    tri3_xyze(1, 0) = 0.08965501613221139865;
    tri3_xyze(2, 0) = 0.84221639627494804259;
    nids.push_back(7857);
    tri3_xyze(0, 1) = 0.10980452120246404768;
    tri3_xyze(1, 1) = 0.063395669874476562322;
    tri3_xyze(2, 1) = 0.84221639627494881974;
    nids.push_back(8107);
    tri3_xyze(0, 2) = 0.10172759137977901989;
    tri3_xyze(1, 2) = 0.078058326287794216336;
    tri3_xyze(2, 2) = 0.84037102920679274121;
    nids.push_back(-1134);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.11420380878365363198;
    tri3_xyze(1, 0) = 0.065935599743718723009;
    tri3_xyze(2, 0) = 0.8385256621386373288;
    nids.push_back(8105);
    tri3_xyze(0, 1) = 0.093247019400778799469;
    tri3_xyze(1, 1) = 0.093247019400770181363;
    tri3_xyze(2, 1) = 0.83852566213863666267;
    nids.push_back(7855);
    tri3_xyze(0, 2) = 0.10172759137977901989;
    tri3_xyze(1, 2) = 0.078058326287794216336;
    tri3_xyze(2, 2) = 0.84037102920679274121;
    nids.push_back(-1134);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0896550161322196143;
    tri3_xyze(1, 0) = 0.08965501613221139865;
    tri3_xyze(2, 0) = 0.84221639627494804259;
    nids.push_back(7857);
    tri3_xyze(0, 1) = 0.085764249336388140299;
    tri3_xyze(1, 1) = 0.085764249336379549948;
    tri3_xyze(2, 1) = 0.84524135262314836758;
    nids.push_back(7859);
    tri3_xyze(0, 2) = 0.097565777798640077045;
    tri3_xyze(1, 2) = 0.074864854408048964185;
    tri3_xyze(2, 2) = 0.84372887444904864918;
    nids.push_back(-1135);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.085764249336388140299;
    tri3_xyze(1, 0) = 0.085764249336379549948;
    tri3_xyze(2, 0) = 0.84524135262314836758;
    nids.push_back(7859);
    tri3_xyze(0, 1) = 0.1050393245234885059;
    tri3_xyze(1, 1) = 0.060644482289128352759;
    tri3_xyze(2, 1) = 0.84524135262314914474;
    nids.push_back(8109);
    tri3_xyze(0, 2) = 0.097565777798640077045;
    tri3_xyze(1, 2) = 0.074864854408048964185;
    tri3_xyze(2, 2) = 0.84372887444904864918;
    nids.push_back(-1135);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.1050393245234885059;
    tri3_xyze(1, 0) = 0.060644482289128352759;
    tri3_xyze(2, 0) = 0.84524135262314914474;
    nids.push_back(8109);
    tri3_xyze(0, 1) = 0.10980452120246404768;
    tri3_xyze(1, 1) = 0.063395669874476562322;
    tri3_xyze(2, 1) = 0.84221639627494881974;
    nids.push_back(8107);
    tri3_xyze(0, 2) = 0.097565777798640077045;
    tri3_xyze(1, 2) = 0.074864854408048964185;
    tri3_xyze(2, 2) = 0.84372887444904864918;
    nids.push_back(-1135);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.10980452120246404768;
    tri3_xyze(1, 0) = 0.063395669874476562322;
    tri3_xyze(2, 0) = 0.84221639627494881974;
    nids.push_back(8107);
    tri3_xyze(0, 1) = 0.0896550161322196143;
    tri3_xyze(1, 1) = 0.08965501613221139865;
    tri3_xyze(2, 1) = 0.84221639627494804259;
    nids.push_back(7857);
    tri3_xyze(0, 2) = 0.097565777798640077045;
    tri3_xyze(1, 2) = 0.074864854408048964185;
    tri3_xyze(2, 2) = 0.84372887444904864918;
    nids.push_back(-1135);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.085764249336388140299;
    tri3_xyze(1, 0) = 0.085764249336379549948;
    tri3_xyze(2, 0) = 0.84524135262314836758;
    nids.push_back(7859);
    tri3_xyze(0, 1) = 0.081636078729885278027;
    tri3_xyze(1, 1) = 0.081636078729876770943;
    tri3_xyze(2, 1) = 0.84755282581460522628;
    nids.push_back(7861);
    tri3_xyze(0, 2) = 0.093105755333676681884;
    tri3_xyze(1, 2) = 0.071442558803689978908;
    tri3_xyze(2, 2) = 0.84639708921887713;
    nids.push_back(-1136);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.081636078729885278027;
    tri3_xyze(1, 0) = 0.081636078729876770943;
    tri3_xyze(2, 0) = 0.84755282581460522628;
    nids.push_back(7861);
    tri3_xyze(0, 1) = 0.09998336874494478943;
    tri3_xyze(1, 1) = 0.05772542485937524892;
    tri3_xyze(2, 1) = 0.84755282581460589242;
    nids.push_back(8111);
    tri3_xyze(0, 2) = 0.093105755333676681884;
    tri3_xyze(1, 2) = 0.071442558803689978908;
    tri3_xyze(2, 2) = 0.84639708921887713;
    nids.push_back(-1136);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.09998336874494478943;
    tri3_xyze(1, 0) = 0.05772542485937524892;
    tri3_xyze(2, 0) = 0.84755282581460589242;
    nids.push_back(8111);
    tri3_xyze(0, 1) = 0.1050393245234885059;
    tri3_xyze(1, 1) = 0.060644482289128352759;
    tri3_xyze(2, 1) = 0.84524135262314914474;
    nids.push_back(8109);
    tri3_xyze(0, 2) = 0.093105755333676681884;
    tri3_xyze(1, 2) = 0.071442558803689978908;
    tri3_xyze(2, 2) = 0.84639708921887713;
    nids.push_back(-1136);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.1050393245234885059;
    tri3_xyze(1, 0) = 0.060644482289128352759;
    tri3_xyze(2, 0) = 0.84524135262314914474;
    nids.push_back(8109);
    tri3_xyze(0, 1) = 0.085764249336388140299;
    tri3_xyze(1, 1) = 0.085764249336379549948;
    tri3_xyze(2, 1) = 0.84524135262314836758;
    nids.push_back(7859);
    tri3_xyze(0, 2) = 0.093105755333676681884;
    tri3_xyze(1, 2) = 0.071442558803689978908;
    tri3_xyze(2, 2) = 0.84639708921887713;
    nids.push_back(-1136);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.081636078729885278027;
    tri3_xyze(1, 0) = 0.081636078729876770943;
    tri3_xyze(2, 0) = 0.84755282581460522628;
    nids.push_back(7861);
    tri3_xyze(0, 1) = 0.077335608029225483273;
    tri3_xyze(1, 1) = 0.077335608029217170478;
    tri3_xyze(2, 1) = 0.84911436253628203019;
    nids.push_back(7863);
    tri3_xyze(0, 2) = 0.088417861203449368834;
    tri3_xyze(1, 2) = 0.067845411120778492675;
    tri3_xyze(2, 2) = 0.8483335941754439613;
    nids.push_back(-1137);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.077335608029225483273;
    tri3_xyze(1, 0) = 0.077335608029217170478;
    tri3_xyze(2, 0) = 0.84911436253628203019;
    nids.push_back(7863);
    tri3_xyze(0, 1) = 0.094716389309741882974;
    tri3_xyze(1, 1) = 0.05468453286464477342;
    tri3_xyze(2, 1) = 0.84911436253628269633;
    nids.push_back(8113);
    tri3_xyze(0, 2) = 0.088417861203449368834;
    tri3_xyze(1, 2) = 0.067845411120778492675;
    tri3_xyze(2, 2) = 0.8483335941754439613;
    nids.push_back(-1137);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.094716389309741882974;
    tri3_xyze(1, 0) = 0.05468453286464477342;
    tri3_xyze(2, 0) = 0.84911436253628269633;
    nids.push_back(8113);
    tri3_xyze(0, 1) = 0.09998336874494478943;
    tri3_xyze(1, 1) = 0.05772542485937524892;
    tri3_xyze(2, 1) = 0.84755282581460589242;
    nids.push_back(8111);
    tri3_xyze(0, 2) = 0.088417861203449368834;
    tri3_xyze(1, 2) = 0.067845411120778492675;
    tri3_xyze(2, 2) = 0.8483335941754439613;
    nids.push_back(-1137);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.09998336874494478943;
    tri3_xyze(1, 0) = 0.05772542485937524892;
    tri3_xyze(2, 0) = 0.84755282581460589242;
    nids.push_back(8111);
    tri3_xyze(0, 1) = 0.081636078729885278027;
    tri3_xyze(1, 1) = 0.081636078729876770943;
    tri3_xyze(2, 1) = 0.84755282581460522628;
    nids.push_back(7861);
    tri3_xyze(0, 2) = 0.088417861203449368834;
    tri3_xyze(1, 2) = 0.067845411120778492675;
    tri3_xyze(2, 2) = 0.8483335941754439613;
    nids.push_back(-1137);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.077335608029225483273;
    tri3_xyze(1, 0) = 0.077335608029217170478;
    tri3_xyze(2, 0) = 0.84911436253628203019;
    nids.push_back(7863);
    tri3_xyze(0, 1) = 0.072930658226334904581;
    tri3_xyze(1, 1) = 0.072930658226326466886;
    tri3_xyze(2, 1) = 0.84990133642126131086;
    nids.push_back(7865);
    tri3_xyze(0, 2) = 0.083576026298804145798;
    tri3_xyze(1, 2) = 0.064130140527105736203;
    tri3_xyze(2, 2) = 0.84950784947877200359;
    nids.push_back(-1138);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.072930658226334904581;
    tri3_xyze(1, 0) = 0.072930658226326466886;
    tri3_xyze(2, 0) = 0.84990133642126131086;
    nids.push_back(7865);
    tri3_xyze(0, 1) = 0.089321449629914353996;
    tri3_xyze(1, 1) = 0.051569762988234575662;
    tri3_xyze(2, 1) = 0.84990133642126197699;
    nids.push_back(8115);
    tri3_xyze(0, 2) = 0.083576026298804145798;
    tri3_xyze(1, 2) = 0.064130140527105736203;
    tri3_xyze(2, 2) = 0.84950784947877200359;
    nids.push_back(-1138);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.089321449629914353996;
    tri3_xyze(1, 0) = 0.051569762988234575662;
    tri3_xyze(2, 0) = 0.84990133642126197699;
    nids.push_back(8115);
    tri3_xyze(0, 1) = 0.094716389309741882974;
    tri3_xyze(1, 1) = 0.05468453286464477342;
    tri3_xyze(2, 1) = 0.84911436253628269633;
    nids.push_back(8113);
    tri3_xyze(0, 2) = 0.083576026298804145798;
    tri3_xyze(1, 2) = 0.064130140527105736203;
    tri3_xyze(2, 2) = 0.84950784947877200359;
    nids.push_back(-1138);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.094716389309741882974;
    tri3_xyze(1, 0) = 0.05468453286464477342;
    tri3_xyze(2, 0) = 0.84911436253628269633;
    nids.push_back(8113);
    tri3_xyze(0, 1) = 0.077335608029225483273;
    tri3_xyze(1, 1) = 0.077335608029217170478;
    tri3_xyze(2, 1) = 0.84911436253628203019;
    nids.push_back(7863);
    tri3_xyze(0, 2) = 0.083576026298804145798;
    tri3_xyze(1, 2) = 0.064130140527105736203;
    tri3_xyze(2, 2) = 0.84950784947877200359;
    nids.push_back(-1138);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.072930658226334904581;
    tri3_xyze(1, 0) = 0.072930658226326466886;
    tri3_xyze(2, 0) = 0.84990133642126131086;
    nids.push_back(7865);
    tri3_xyze(0, 1) = 0.068490698010994513489;
    tri3_xyze(1, 1) = 0.068490698010986075794;
    tri3_xyze(2, 1) = 0.84990133642126153291;
    nids.push_back(7867);
    tri3_xyze(0, 2) = 0.078656609248559272718;
    tri3_xyze(1, 2) = 0.060355339059328980134;
    tri3_xyze(2, 2) = 0.84990133642126175495;
    nids.push_back(-1139);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.068490698010994513489;
    tri3_xyze(1, 0) = 0.068490698010986075794;
    tri3_xyze(2, 0) = 0.84990133642126153291;
    nids.push_back(7867);
    tri3_xyze(0, 1) = 0.083883631126993332683;
    tri3_xyze(1, 1) = 0.048430237011768788313;
    tri3_xyze(2, 1) = 0.84990133642126208802;
    nids.push_back(8117);
    tri3_xyze(0, 2) = 0.078656609248559272718;
    tri3_xyze(1, 2) = 0.060355339059328980134;
    tri3_xyze(2, 2) = 0.84990133642126175495;
    nids.push_back(-1139);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.083883631126993332683;
    tri3_xyze(1, 0) = 0.048430237011768788313;
    tri3_xyze(2, 0) = 0.84990133642126208802;
    nids.push_back(8117);
    tri3_xyze(0, 1) = 0.089321449629914353996;
    tri3_xyze(1, 1) = 0.051569762988234575662;
    tri3_xyze(2, 1) = 0.84990133642126197699;
    nids.push_back(8115);
    tri3_xyze(0, 2) = 0.078656609248559272718;
    tri3_xyze(1, 2) = 0.060355339059328980134;
    tri3_xyze(2, 2) = 0.84990133642126175495;
    nids.push_back(-1139);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.089321449629914353996;
    tri3_xyze(1, 0) = 0.051569762988234575662;
    tri3_xyze(2, 0) = 0.84990133642126197699;
    nids.push_back(8115);
    tri3_xyze(0, 1) = 0.072930658226334904581;
    tri3_xyze(1, 1) = 0.072930658226326466886;
    tri3_xyze(2, 1) = 0.84990133642126131086;
    nids.push_back(7865);
    tri3_xyze(0, 2) = 0.078656609248559272718;
    tri3_xyze(1, 2) = 0.060355339059328980134;
    tri3_xyze(2, 2) = 0.84990133642126175495;
    nids.push_back(-1139);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.068490698010994513489;
    tri3_xyze(1, 0) = 0.068490698010986075794;
    tri3_xyze(2, 0) = 0.84990133642126153291;
    nids.push_back(7867);
    tri3_xyze(0, 1) = 0.064085748208103893164;
    tri3_xyze(1, 1) = 0.064085748208095552614;
    tri3_xyze(2, 1) = 0.84911436253628247428;
    nids.push_back(7869);
    tri3_xyze(0, 2) = 0.073737192198314371883;
    tri3_xyze(1, 2) = 0.056580537591552231003;
    tri3_xyze(2, 2) = 0.84950784947877222564;
    nids.push_back(-1140);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.064085748208103893164;
    tri3_xyze(1, 0) = 0.064085748208095552614;
    tri3_xyze(2, 0) = 0.84911436253628247428;
    nids.push_back(7869);
    tri3_xyze(0, 1) = 0.078488691447165720438;
    tri3_xyze(1, 1) = 0.045315467135358493411;
    tri3_xyze(2, 1) = 0.84911436253628302939;
    nids.push_back(8119);
    tri3_xyze(0, 2) = 0.073737192198314371883;
    tri3_xyze(1, 2) = 0.056580537591552231003;
    tri3_xyze(2, 2) = 0.84950784947877222564;
    nids.push_back(-1140);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.078488691447165720438;
    tri3_xyze(1, 0) = 0.045315467135358493411;
    tri3_xyze(2, 0) = 0.84911436253628302939;
    nids.push_back(8119);
    tri3_xyze(0, 1) = 0.083883631126993332683;
    tri3_xyze(1, 1) = 0.048430237011768788313;
    tri3_xyze(2, 1) = 0.84990133642126208802;
    nids.push_back(8117);
    tri3_xyze(0, 2) = 0.073737192198314371883;
    tri3_xyze(1, 2) = 0.056580537591552231003;
    tri3_xyze(2, 2) = 0.84950784947877222564;
    nids.push_back(-1140);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.083883631126993332683;
    tri3_xyze(1, 0) = 0.048430237011768788313;
    tri3_xyze(2, 0) = 0.84990133642126208802;
    nids.push_back(8117);
    tri3_xyze(0, 1) = 0.068490698010994513489;
    tri3_xyze(1, 1) = 0.068490698010986075794;
    tri3_xyze(2, 1) = 0.84990133642126153291;
    nids.push_back(7867);
    tri3_xyze(0, 2) = 0.073737192198314371883;
    tri3_xyze(1, 2) = 0.056580537591552231003;
    tri3_xyze(2, 2) = 0.84950784947877222564;
    nids.push_back(-1140);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.064085748208103893164;
    tri3_xyze(1, 0) = 0.064085748208095552614;
    tri3_xyze(2, 0) = 0.84911436253628247428;
    nids.push_back(7869);
    tri3_xyze(0, 1) = 0.059785277507444195555;
    tri3_xyze(1, 1) = 0.059785277507435799493;
    tri3_xyze(2, 1) = 0.84755282581460589242;
    nids.push_back(7871);
    tri3_xyze(0, 2) = 0.068895357293669190479;
    tri3_xyze(1, 2) = 0.052865266997879460653;
    tri3_xyze(2, 2) = 0.84833359417544440539;
    nids.push_back(-1141);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.059785277507444195555;
    tri3_xyze(1, 0) = 0.059785277507435799493;
    tri3_xyze(2, 0) = 0.84755282581460589242;
    nids.push_back(7871);
    tri3_xyze(0, 1) = 0.073221712011962938882;
    tri3_xyze(1, 1) = 0.042274575140627983216;
    tri3_xyze(2, 1) = 0.84755282581460644753;
    nids.push_back(8121);
    tri3_xyze(0, 2) = 0.068895357293669190479;
    tri3_xyze(1, 2) = 0.052865266997879460653;
    tri3_xyze(2, 2) = 0.84833359417544440539;
    nids.push_back(-1141);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.073221712011962938882;
    tri3_xyze(1, 0) = 0.042274575140627983216;
    tri3_xyze(2, 0) = 0.84755282581460644753;
    nids.push_back(8121);
    tri3_xyze(0, 1) = 0.078488691447165720438;
    tri3_xyze(1, 1) = 0.045315467135358493411;
    tri3_xyze(2, 1) = 0.84911436253628302939;
    nids.push_back(8119);
    tri3_xyze(0, 2) = 0.068895357293669190479;
    tri3_xyze(1, 2) = 0.052865266997879460653;
    tri3_xyze(2, 2) = 0.84833359417544440539;
    nids.push_back(-1141);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.078488691447165720438;
    tri3_xyze(1, 0) = 0.045315467135358493411;
    tri3_xyze(2, 0) = 0.84911436253628302939;
    nids.push_back(8119);
    tri3_xyze(0, 1) = 0.064085748208103893164;
    tri3_xyze(1, 1) = 0.064085748208095552614;
    tri3_xyze(2, 1) = 0.84911436253628247428;
    nids.push_back(7869);
    tri3_xyze(0, 2) = 0.068895357293669190479;
    tri3_xyze(1, 2) = 0.052865266997879460653;
    tri3_xyze(2, 2) = 0.84833359417544440539;
    nids.push_back(-1141);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.059785277507444195555;
    tri3_xyze(1, 0) = 0.059785277507435799493;
    tri3_xyze(2, 0) = 0.84755282581460589242;
    nids.push_back(7871);
    tri3_xyze(0, 1) = 0.055657106900941374916;
    tri3_xyze(1, 1) = 0.055657106900932902527;
    tri3_xyze(2, 1) = 0.84524135262314925576;
    nids.push_back(7873);
    tri3_xyze(0, 2) = 0.064207463163441905185;
    tri3_xyze(1, 2) = 0.049268119314967828704;
    tri3_xyze(2, 2) = 0.84639708921887779614;
    nids.push_back(-1142);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.055657106900941374916;
    tri3_xyze(1, 0) = 0.055657106900932902527;
    tri3_xyze(2, 0) = 0.84524135262314925576;
    nids.push_back(7873);
    tri3_xyze(0, 1) = 0.068165756233419139143;
    tri3_xyze(1, 1) = 0.039355517710874629578;
    tri3_xyze(2, 1) = 0.84524135262314981087;
    nids.push_back(8123);
    tri3_xyze(0, 2) = 0.064207463163441905185;
    tri3_xyze(1, 2) = 0.049268119314967828704;
    tri3_xyze(2, 2) = 0.84639708921887779614;
    nids.push_back(-1142);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.073221712011962938882;
    tri3_xyze(1, 0) = 0.042274575140627983216;
    tri3_xyze(2, 0) = 0.84755282581460644753;
    nids.push_back(8121);
    tri3_xyze(0, 1) = 0.059785277507444195555;
    tri3_xyze(1, 1) = 0.059785277507435799493;
    tri3_xyze(2, 1) = 0.84755282581460589242;
    nids.push_back(7871);
    tri3_xyze(0, 2) = 0.064207463163441905185;
    tri3_xyze(1, 2) = 0.049268119314967828704;
    tri3_xyze(2, 2) = 0.84639708921887779614;
    nids.push_back(-1142);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.055657106900941374916;
    tri3_xyze(1, 0) = 0.055657106900932902527;
    tri3_xyze(2, 0) = 0.84524135262314925576;
    nids.push_back(7873);
    tri3_xyze(0, 1) = 0.051766340105109574787;
    tri3_xyze(1, 1) = 0.051766340105100998314;
    tri3_xyze(2, 1) = 0.84221639627494915281;
    nids.push_back(7875);
    tri3_xyze(0, 2) = 0.059747440698478385124;
    tri3_xyze(1, 2) = 0.045845823710608746282;
    tri3_xyze(2, 2) = 0.84372887444904953735;
    nids.push_back(-1143);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051766340105109574787;
    tri3_xyze(1, 0) = 0.051766340105100998314;
    tri3_xyze(2, 0) = 0.84221639627494915281;
    nids.push_back(7875);
    tri3_xyze(0, 1) = 0.063400559554443458588;
    tri3_xyze(1, 1) = 0.036604330125526468587;
    tri3_xyze(2, 1) = 0.84221639627494970792;
    nids.push_back(8125);
    tri3_xyze(0, 2) = 0.059747440698478385124;
    tri3_xyze(1, 2) = 0.045845823710608746282;
    tri3_xyze(2, 2) = 0.84372887444904953735;
    nids.push_back(-1143);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.068165756233419139143;
    tri3_xyze(1, 0) = 0.039355517710874629578;
    tri3_xyze(2, 0) = 0.84524135262314981087;
    nids.push_back(8123);
    tri3_xyze(0, 1) = 0.055657106900941374916;
    tri3_xyze(1, 1) = 0.055657106900932902527;
    tri3_xyze(2, 1) = 0.84524135262314925576;
    nids.push_back(7873);
    tri3_xyze(0, 2) = 0.059747440698478385124;
    tri3_xyze(1, 2) = 0.045845823710608746282;
    tri3_xyze(2, 2) = 0.84372887444904953735;
    nids.push_back(-1143);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051766340105109574787;
    tri3_xyze(1, 0) = 0.051766340105100998314;
    tri3_xyze(2, 0) = 0.84221639627494915281;
    nids.push_back(7875);
    tri3_xyze(0, 1) = 0.04817433683655062554;
    tri3_xyze(1, 1) = 0.048174336836541924167;
    tri3_xyze(2, 1) = 0.83852566213863799494;
    nids.push_back(7877);
    tri3_xyze(0, 2) = 0.055585627117339414527;
    tri3_xyze(1, 2) = 0.042652351830863362292;
    tri3_xyze(2, 2) = 0.84037102920679385143;
    nids.push_back(-1144);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.04817433683655062554;
    tri3_xyze(1, 0) = 0.048174336836541924167;
    tri3_xyze(2, 0) = 0.83852566213863799494;
    nids.push_back(7877);
    tri3_xyze(0, 1) = 0.05900127197325401307;
    tri3_xyze(1, 1) = 0.034064400256284051161;
    tri3_xyze(2, 1) = 0.83852566213863843902;
    nids.push_back(8127);
    tri3_xyze(0, 2) = 0.055585627117339414527;
    tri3_xyze(1, 2) = 0.042652351830863362292;
    tri3_xyze(2, 2) = 0.84037102920679385143;
    nids.push_back(-1144);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.063400559554443458588;
    tri3_xyze(1, 0) = 0.036604330125526468587;
    tri3_xyze(2, 0) = 0.84221639627494970792;
    nids.push_back(8125);
    tri3_xyze(0, 1) = 0.051766340105109574787;
    tri3_xyze(1, 1) = 0.051766340105100998314;
    tri3_xyze(2, 1) = 0.84221639627494915281;
    nids.push_back(7875);
    tri3_xyze(0, 2) = 0.055585627117339414527;
    tri3_xyze(1, 2) = 0.042652351830863362292;
    tri3_xyze(2, 2) = 0.84037102920679385143;
    nids.push_back(-1144);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.04817433683655062554;
    tri3_xyze(1, 0) = 0.048174336836541924167;
    tri3_xyze(2, 0) = 0.83852566213863799494;
    nids.push_back(7877);
    tri3_xyze(0, 1) = 0.044937745132568200634;
    tri3_xyze(1, 1) = 0.044937745132559360484;
    tri3_xyze(2, 1) = 0.83422735529628300188;
    nids.push_back(7879);
    tri3_xyze(0, 2) = 0.051787656706346454905;
    tri3_xyze(1, 2) = 0.039738066634962786183;
    tri3_xyze(2, 2) = 0.83637650871746072045;
    nids.push_back(-1145);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.044937745132568200634;
    tri3_xyze(1, 0) = 0.044937745132559360484;
    tri3_xyze(2, 0) = 0.83422735529628300188;
    nids.push_back(7879);
    tri3_xyze(0, 1) = 0.055037272883012987312;
    tri3_xyze(1, 1) = 0.0317757843144658228;
    tri3_xyze(2, 1) = 0.83422735529628344597;
    nids.push_back(8129);
    tri3_xyze(0, 2) = 0.051787656706346454905;
    tri3_xyze(1, 2) = 0.039738066634962786183;
    tri3_xyze(2, 2) = 0.83637650871746072045;
    nids.push_back(-1145);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.05900127197325401307;
    tri3_xyze(1, 0) = 0.034064400256284051161;
    tri3_xyze(2, 0) = 0.83852566213863843902;
    nids.push_back(8127);
    tri3_xyze(0, 1) = 0.04817433683655062554;
    tri3_xyze(1, 1) = 0.048174336836541924167;
    tri3_xyze(2, 1) = 0.83852566213863799494;
    nids.push_back(7877);
    tri3_xyze(0, 2) = 0.051787656706346454905;
    tri3_xyze(1, 2) = 0.039738066634962786183;
    tri3_xyze(2, 2) = 0.83637650871746072045;
    nids.push_back(-1145);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.044937745132568200634;
    tri3_xyze(1, 0) = 0.044937745132559360484;
    tri3_xyze(2, 0) = 0.83422735529628300188;
    nids.push_back(7879);
    tri3_xyze(0, 1) = 0.042107607977780528441;
    tri3_xyze(1, 1) = 0.042107607977771528696;
    tri3_xyze(2, 1) = 0.82938926261447243338;
    nids.push_back(7881);
    tri3_xyze(0, 2) = 0.048413425727678427446;
    tri3_xyze(1, 2) = 0.03714892814135600102;
    tri3_xyze(2, 2) = 0.83180830895537782865;
    nids.push_back(-1146);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042107607977780528441;
    tri3_xyze(1, 0) = 0.042107607977771528696;
    tri3_xyze(2, 0) = 0.82938926261447243338;
    nids.push_back(7881);
    tri3_xyze(0, 1) = 0.051571076917351993396;
    tri3_xyze(1, 1) = 0.029774575140627295572;
    tri3_xyze(2, 1) = 0.82938926261447276644;
    nids.push_back(8131);
    tri3_xyze(0, 2) = 0.048413425727678427446;
    tri3_xyze(1, 2) = 0.03714892814135600102;
    tri3_xyze(2, 2) = 0.83180830895537782865;
    nids.push_back(-1146);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051571076917351993396;
    tri3_xyze(1, 0) = 0.029774575140627295572;
    tri3_xyze(2, 0) = 0.82938926261447276644;
    nids.push_back(8131);
    tri3_xyze(0, 1) = 0.055037272883012987312;
    tri3_xyze(1, 1) = 0.0317757843144658228;
    tri3_xyze(2, 1) = 0.83422735529628344597;
    nids.push_back(8129);
    tri3_xyze(0, 2) = 0.048413425727678427446;
    tri3_xyze(1, 2) = 0.03714892814135600102;
    tri3_xyze(2, 2) = 0.83180830895537782865;
    nids.push_back(-1146);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.055037272883012987312;
    tri3_xyze(1, 0) = 0.0317757843144658228;
    tri3_xyze(2, 0) = 0.83422735529628344597;
    nids.push_back(8129);
    tri3_xyze(0, 1) = 0.044937745132568200634;
    tri3_xyze(1, 1) = 0.044937745132559360484;
    tri3_xyze(2, 1) = 0.83422735529628300188;
    nids.push_back(7879);
    tri3_xyze(0, 2) = 0.048413425727678427446;
    tri3_xyze(1, 2) = 0.03714892814135600102;
    tri3_xyze(2, 2) = 0.83180830895537782865;
    nids.push_back(-1146);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042107607977780528441;
    tri3_xyze(1, 0) = 0.042107607977771528696;
    tri3_xyze(2, 0) = 0.82938926261447243338;
    nids.push_back(7881);
    tri3_xyze(0, 1) = 0.039728558325760408854;
    tri3_xyze(1, 1) = 0.039728558325751242575;
    tri3_xyze(2, 1) = 0.82408768370493457045;
    nids.push_back(7883);
    tri3_xyze(0, 2) = 0.045516147819536827146;
    tri3_xyze(1, 2) = 0.034925768610763595345;
    tri3_xyze(2, 2) = 0.82673847315970372396;
    nids.push_back(-1147);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.039728558325760408854;
    tri3_xyze(1, 0) = 0.039728558325751242575;
    tri3_xyze(2, 0) = 0.82408768370493457045;
    nids.push_back(7883);
    tri3_xyze(0, 1) = 0.048657348057254377893;
    tri3_xyze(1, 1) = 0.028092332998904304131;
    tri3_xyze(2, 1) = 0.82408768370493490352;
    nids.push_back(8133);
    tri3_xyze(0, 2) = 0.045516147819536827146;
    tri3_xyze(1, 2) = 0.034925768610763595345;
    tri3_xyze(2, 2) = 0.82673847315970372396;
    nids.push_back(-1147);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048657348057254377893;
    tri3_xyze(1, 0) = 0.028092332998904304131;
    tri3_xyze(2, 0) = 0.82408768370493490352;
    nids.push_back(8133);
    tri3_xyze(0, 1) = 0.051571076917351993396;
    tri3_xyze(1, 1) = 0.029774575140627295572;
    tri3_xyze(2, 1) = 0.82938926261447276644;
    nids.push_back(8131);
    tri3_xyze(0, 2) = 0.045516147819536827146;
    tri3_xyze(1, 2) = 0.034925768610763595345;
    tri3_xyze(2, 2) = 0.82673847315970372396;
    nids.push_back(-1147);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051571076917351993396;
    tri3_xyze(1, 0) = 0.029774575140627295572;
    tri3_xyze(2, 0) = 0.82938926261447276644;
    nids.push_back(8131);
    tri3_xyze(0, 1) = 0.042107607977780528441;
    tri3_xyze(1, 1) = 0.042107607977771528696;
    tri3_xyze(2, 1) = 0.82938926261447243338;
    nids.push_back(7881);
    tri3_xyze(0, 2) = 0.045516147819536827146;
    tri3_xyze(1, 2) = 0.034925768610763595345;
    tri3_xyze(2, 2) = 0.82673847315970372396;
    nids.push_back(-1147);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.039728558325760408854;
    tri3_xyze(1, 0) = 0.039728558325751242575;
    tri3_xyze(2, 0) = 0.82408768370493457045;
    nids.push_back(7883);
    tri3_xyze(0, 1) = 0.037838115210695681678;
    tri3_xyze(1, 1) = 0.037838115210686328049;
    tri3_xyze(2, 1) = 0.8184062276340827724;
    nids.push_back(7885);
    tri3_xyze(0, 2) = 0.043141514785282969457;
    tri3_xyze(1, 2) = 0.033103648597034046608;
    tri3_xyze(2, 2) = 0.82124695566950878245;
    nids.push_back(-1148);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.037838115210695681678;
    tri3_xyze(1, 0) = 0.037838115210686328049;
    tri3_xyze(2, 0) = 0.8184062276340827724;
    nids.push_back(7885);
    tri3_xyze(0, 1) = 0.046342037547421395527;
    tri3_xyze(1, 1) = 0.026755587852794311676;
    tri3_xyze(2, 1) = 0.81840622763408310547;
    nids.push_back(8135);
    tri3_xyze(0, 2) = 0.043141514785282969457;
    tri3_xyze(1, 2) = 0.033103648597034046608;
    tri3_xyze(2, 2) = 0.82124695566950878245;
    nids.push_back(-1148);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.046342037547421395527;
    tri3_xyze(1, 0) = 0.026755587852794311676;
    tri3_xyze(2, 0) = 0.81840622763408310547;
    nids.push_back(8135);
    tri3_xyze(0, 1) = 0.048657348057254377893;
    tri3_xyze(1, 1) = 0.028092332998904304131;
    tri3_xyze(2, 1) = 0.82408768370493490352;
    nids.push_back(8133);
    tri3_xyze(0, 2) = 0.043141514785282969457;
    tri3_xyze(1, 2) = 0.033103648597034046608;
    tri3_xyze(2, 2) = 0.82124695566950878245;
    nids.push_back(-1148);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048657348057254377893;
    tri3_xyze(1, 0) = 0.028092332998904304131;
    tri3_xyze(2, 0) = 0.82408768370493490352;
    nids.push_back(8133);
    tri3_xyze(0, 1) = 0.039728558325760408854;
    tri3_xyze(1, 1) = 0.039728558325751242575;
    tri3_xyze(2, 1) = 0.82408768370493457045;
    nids.push_back(7883);
    tri3_xyze(0, 2) = 0.043141514785282969457;
    tri3_xyze(1, 2) = 0.033103648597034046608;
    tri3_xyze(2, 2) = 0.82124695566950878245;
    nids.push_back(-1148);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.037838115210695681678;
    tri3_xyze(1, 0) = 0.037838115210686328049;
    tri3_xyze(2, 0) = 0.8184062276340827724;
    nids.push_back(7885);
    tri3_xyze(0, 1) = 0.036466092049806886344;
    tri3_xyze(1, 1) = 0.036466092049797275976;
    tri3_xyze(2, 1) = 0.81243449435809167447;
    nids.push_back(7887);
    tri3_xyze(0, 2) = 0.041326976006404066155;
    tri3_xyze(1, 2) = 0.031711304021265644681;
    tri3_xyze(2, 2) = 0.81542036099608738997;
    nids.push_back(-1149);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036466092049806886344;
    tri3_xyze(1, 0) = 0.036466092049797275976;
    tri3_xyze(2, 0) = 0.81243449435809167447;
    nids.push_back(7887);
    tri3_xyze(0, 1) = 0.044661659217692294133;
    tri3_xyze(1, 1) = 0.025785420971784663025;
    tri3_xyze(2, 1) = 0.81243449435809200754;
    nids.push_back(8137);
    tri3_xyze(0, 2) = 0.041326976006404066155;
    tri3_xyze(1, 2) = 0.031711304021265644681;
    tri3_xyze(2, 2) = 0.81542036099608738997;
    nids.push_back(-1149);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.044661659217692294133;
    tri3_xyze(1, 0) = 0.025785420971784663025;
    tri3_xyze(2, 0) = 0.81243449435809200754;
    nids.push_back(8137);
    tri3_xyze(0, 1) = 0.046342037547421395527;
    tri3_xyze(1, 1) = 0.026755587852794311676;
    tri3_xyze(2, 1) = 0.81840622763408310547;
    nids.push_back(8135);
    tri3_xyze(0, 2) = 0.041326976006404066155;
    tri3_xyze(1, 2) = 0.031711304021265644681;
    tri3_xyze(2, 2) = 0.81542036099608738997;
    nids.push_back(-1149);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.046342037547421395527;
    tri3_xyze(1, 0) = 0.026755587852794311676;
    tri3_xyze(2, 0) = 0.81840622763408310547;
    nids.push_back(8135);
    tri3_xyze(0, 1) = 0.037838115210695681678;
    tri3_xyze(1, 1) = 0.037838115210686328049;
    tri3_xyze(2, 1) = 0.8184062276340827724;
    nids.push_back(7885);
    tri3_xyze(0, 2) = 0.041326976006404066155;
    tri3_xyze(1, 2) = 0.031711304021265644681;
    tri3_xyze(2, 2) = 0.81542036099608738997;
    nids.push_back(-1149);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036466092049806886344;
    tri3_xyze(1, 0) = 0.036466092049797275976;
    tri3_xyze(2, 0) = 0.81243449435809167447;
    nids.push_back(7887);
    tri3_xyze(0, 1) = 0.035634126467948050443;
    tri3_xyze(1, 1) = 0.035634126467938308236;
    tri3_xyze(2, 1) = 0.80626666167806415597;
    nids.push_back(7889);
    tri3_xyze(0, 2) = 0.040101147843395919257;
    tri3_xyze(1, 2) = 0.030770692989164627473;
    tri3_xyze(2, 2) = 0.80935057801807808175;
    nids.push_back(-1150);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.035634126467948050443;
    tri3_xyze(1, 0) = 0.035634126467938308236;
    tri3_xyze(2, 0) = 0.80626666167806415597;
    nids.push_back(7889);
    tri3_xyze(0, 1) = 0.043642713638136439169;
    tri3_xyze(1, 1) = 0.025197132467138269596;
    tri3_xyze(2, 1) = 0.80626666167806460006;
    nids.push_back(8139);
    tri3_xyze(0, 2) = 0.040101147843395919257;
    tri3_xyze(1, 2) = 0.030770692989164627473;
    tri3_xyze(2, 2) = 0.80935057801807808175;
    nids.push_back(-1150);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.043642713638136439169;
    tri3_xyze(1, 0) = 0.025197132467138269596;
    tri3_xyze(2, 0) = 0.80626666167806460006;
    nids.push_back(8139);
    tri3_xyze(0, 1) = 0.044661659217692294133;
    tri3_xyze(1, 1) = 0.025785420971784663025;
    tri3_xyze(2, 1) = 0.81243449435809200754;
    nids.push_back(8137);
    tri3_xyze(0, 2) = 0.040101147843395919257;
    tri3_xyze(1, 2) = 0.030770692989164627473;
    tri3_xyze(2, 2) = 0.80935057801807808175;
    nids.push_back(-1150);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.044661659217692294133;
    tri3_xyze(1, 0) = 0.025785420971784663025;
    tri3_xyze(2, 0) = 0.81243449435809200754;
    nids.push_back(8137);
    tri3_xyze(0, 1) = 0.036466092049806886344;
    tri3_xyze(1, 1) = 0.036466092049797275976;
    tri3_xyze(2, 1) = 0.81243449435809167447;
    nids.push_back(7887);
    tri3_xyze(0, 2) = 0.040101147843395919257;
    tri3_xyze(1, 2) = 0.030770692989164627473;
    tri3_xyze(2, 2) = 0.80935057801807808175;
    nids.push_back(-1150);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.043642713638136439169;
    tri3_xyze(1, 0) = 0.025197132467138269596;
    tri3_xyze(2, 0) = 0.80626666167806460006;
    nids.push_back(8139);
    tri3_xyze(0, 1) = 0.043301270189231896857;
    tri3_xyze(1, 1) = 0.025000000000000056899;
    tri3_xyze(2, 1) = 0.79999999999984927612;
    nids.push_back(8041);
    tri3_xyze(0, 2) = 0.04597934928467814486;
    tri3_xyze(1, 2) = 0.01904527006279937984;
    tri3_xyze(2, 2) = 0.80313333083895710462;
    nids.push_back(-1152);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.043301270189231896857;
    tri3_xyze(1, 0) = 0.025000000000000056899;
    tri3_xyze(2, 0) = 0.79999999999984927612;
    nids.push_back(8041);
    tri3_xyze(0, 1) = 0.048296291314463390965;
    tri3_xyze(1, 1) = 0.012940952255126120241;
    tri3_xyze(2, 1) = 0.79999999999984960919;
    nids.push_back(8291);
    tri3_xyze(0, 2) = 0.04597934928467814486;
    tri3_xyze(1, 2) = 0.01904527006279937984;
    tri3_xyze(2, 2) = 0.80313333083895710462;
    nids.push_back(-1152);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048296291314463390965;
    tri3_xyze(1, 0) = 0.012940952255126120241;
    tri3_xyze(2, 0) = 0.79999999999984960919;
    nids.push_back(8291);
    tri3_xyze(0, 1) = 0.04867712199688085245;
    tri3_xyze(1, 1) = 0.013042995528933058746;
    tri3_xyze(2, 1) = 0.80626666167806482211;
    nids.push_back(8389);
    tri3_xyze(0, 2) = 0.04597934928467814486;
    tri3_xyze(1, 2) = 0.01904527006279937984;
    tri3_xyze(2, 2) = 0.80313333083895710462;
    nids.push_back(-1152);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.04867712199688085245;
    tri3_xyze(1, 0) = 0.013042995528933058746;
    tri3_xyze(2, 0) = 0.80626666167806482211;
    nids.push_back(8389);
    tri3_xyze(0, 1) = 0.043642713638136439169;
    tri3_xyze(1, 1) = 0.025197132467138269596;
    tri3_xyze(2, 1) = 0.80626666167806460006;
    nids.push_back(8139);
    tri3_xyze(0, 2) = 0.04597934928467814486;
    tri3_xyze(1, 2) = 0.01904527006279937984;
    tri3_xyze(2, 2) = 0.80313333083895710462;
    nids.push_back(-1152);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.1050393245234885059;
    tri3_xyze(1, 0) = 0.060644482289128352759;
    tri3_xyze(2, 0) = 0.84524135262314914474;
    nids.push_back(8109);
    tri3_xyze(0, 1) = 0.09998336874494478943;
    tri3_xyze(1, 1) = 0.05772542485937524892;
    tri3_xyze(2, 1) = 0.84755282581460589242;
    nids.push_back(8111);
    tri3_xyze(0, 2) = 0.10842394850220110625;
    tri3_xyze(1, 2) = 0.044910669955651183471;
    tri3_xyze(2, 2) = 0.84639708921887790716;
    nids.push_back(-1186);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.09998336874494478943;
    tri3_xyze(1, 0) = 0.05772542485937524892;
    tri3_xyze(2, 0) = 0.84755282581460589242;
    nids.push_back(8111);
    tri3_xyze(0, 1) = 0.11151695741036604659;
    tri3_xyze(1, 1) = 0.029880878680482482468;
    tri3_xyze(2, 1) = 0.8475528258146067806;
    nids.push_back(8361);
    tri3_xyze(0, 2) = 0.10842394850220110625;
    tri3_xyze(1, 2) = 0.044910669955651183471;
    tri3_xyze(2, 2) = 0.84639708921887790716;
    nids.push_back(-1186);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.09998336874494478943;
    tri3_xyze(1, 0) = 0.05772542485937524892;
    tri3_xyze(2, 0) = 0.84755282581460589242;
    nids.push_back(8111);
    tri3_xyze(0, 1) = 0.094716389309741882974;
    tri3_xyze(1, 1) = 0.05468453286464477342;
    tri3_xyze(2, 1) = 0.84911436253628269633;
    nids.push_back(8113);
    tri3_xyze(0, 2) = 0.10296478016252173793;
    tri3_xyze(1, 2) = 0.042649408390078198672;
    tri3_xyze(2, 2) = 0.84833359417544473846;
    nids.push_back(-1187);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.094716389309741882974;
    tri3_xyze(1, 0) = 0.05468453286464477342;
    tri3_xyze(2, 0) = 0.84911436253628269633;
    nids.push_back(8113);
    tri3_xyze(0, 1) = 0.10564240518503423272;
    tri3_xyze(1, 1) = 0.028306797155810303757;
    tri3_xyze(2, 1) = 0.8491143625362835845;
    nids.push_back(8363);
    tri3_xyze(0, 2) = 0.10296478016252173793;
    tri3_xyze(1, 2) = 0.042649408390078198672;
    tri3_xyze(2, 2) = 0.84833359417544473846;
    nids.push_back(-1187);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.11151695741036604659;
    tri3_xyze(1, 0) = 0.029880878680482482468;
    tri3_xyze(2, 0) = 0.8475528258146067806;
    nids.push_back(8361);
    tri3_xyze(0, 1) = 0.09998336874494478943;
    tri3_xyze(1, 1) = 0.05772542485937524892;
    tri3_xyze(2, 1) = 0.84755282581460589242;
    nids.push_back(8111);
    tri3_xyze(0, 2) = 0.10296478016252173793;
    tri3_xyze(1, 2) = 0.042649408390078198672;
    tri3_xyze(2, 2) = 0.84833359417544473846;
    nids.push_back(-1187);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.094716389309741882974;
    tri3_xyze(1, 0) = 0.05468453286464477342;
    tri3_xyze(2, 0) = 0.84911436253628269633;
    nids.push_back(8113);
    tri3_xyze(0, 1) = 0.089321449629914353996;
    tri3_xyze(1, 1) = 0.051569762988234575662;
    tri3_xyze(2, 1) = 0.84990133642126197699;
    nids.push_back(8115);
    tri3_xyze(0, 2) = 0.097326343994145220195;
    tri3_xyze(1, 2) = 0.04031389165856171658;
    tri3_xyze(2, 2) = 0.84950784947877266973;
    nids.push_back(-1188);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.089321449629914353996;
    tri3_xyze(1, 0) = 0.051569762988234575662;
    tri3_xyze(2, 0) = 0.84990133642126197699;
    nids.push_back(8115);
    tri3_xyze(0, 1) = 0.099625131851890452728;
    tri3_xyze(1, 1) = 0.02669447362555721695;
    tri3_xyze(2, 1) = 0.84990133642126275415;
    nids.push_back(8365);
    tri3_xyze(0, 2) = 0.097326343994145220195;
    tri3_xyze(1, 2) = 0.04031389165856171658;
    tri3_xyze(2, 2) = 0.84950784947877266973;
    nids.push_back(-1188);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.10564240518503423272;
    tri3_xyze(1, 0) = 0.028306797155810303757;
    tri3_xyze(2, 0) = 0.8491143625362835845;
    nids.push_back(8363);
    tri3_xyze(0, 1) = 0.094716389309741882974;
    tri3_xyze(1, 1) = 0.05468453286464477342;
    tri3_xyze(2, 1) = 0.84911436253628269633;
    nids.push_back(8113);
    tri3_xyze(0, 2) = 0.097326343994145220195;
    tri3_xyze(1, 2) = 0.04031389165856171658;
    tri3_xyze(2, 2) = 0.84950784947877266973;
    nids.push_back(-1188);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.089321449629914353996;
    tri3_xyze(1, 0) = 0.051569762988234575662;
    tri3_xyze(2, 0) = 0.84990133642126197699;
    nids.push_back(8115);
    tri3_xyze(0, 1) = 0.083883631126993332683;
    tri3_xyze(1, 1) = 0.048430237011768788313;
    tri3_xyze(2, 1) = 0.84990133642126208802;
    nids.push_back(8117);
    tri3_xyze(0, 2) = 0.091597561503685337447;
    tri3_xyze(1, 2) = 0.037940952255127738391;
    tri3_xyze(2, 2) = 0.84990133642126242108;
    nids.push_back(-1189);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.099625131851890452728;
    tri3_xyze(1, 0) = 0.02669447362555721695;
    tri3_xyze(2, 0) = 0.84990133642126275415;
    nids.push_back(8365);
    tri3_xyze(0, 1) = 0.089321449629914353996;
    tri3_xyze(1, 1) = 0.051569762988234575662;
    tri3_xyze(2, 1) = 0.84990133642126197699;
    nids.push_back(8115);
    tri3_xyze(0, 2) = 0.091597561503685337447;
    tri3_xyze(1, 2) = 0.037940952255127738391;
    tri3_xyze(2, 2) = 0.84990133642126242108;
    nids.push_back(-1189);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051571076917351993396;
    tri3_xyze(1, 0) = 0.029774575140627295572;
    tri3_xyze(2, 0) = 0.82938926261447276644;
    nids.push_back(8131);
    tri3_xyze(0, 1) = 0.048657348057254377893;
    tri3_xyze(1, 1) = 0.028092332998904304131;
    tri3_xyze(2, 1) = 0.82408768370493490352;
    nids.push_back(8133);
    tri3_xyze(0, 2) = 0.053004676773391018463;
    tri3_xyze(1, 2) = 0.021955255988737560485;
    tri3_xyze(2, 2) = 0.82673847315970405703;
    nids.push_back(-1197);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048657348057254377893;
    tri3_xyze(1, 0) = 0.028092332998904304131;
    tri3_xyze(2, 0) = 0.82408768370493490352;
    nids.push_back(8133);
    tri3_xyze(0, 1) = 0.054270219928716861846;
    tri3_xyze(1, 1) = 0.014541661602957268312;
    tri3_xyze(2, 1) = 0.82408768370493534761;
    nids.push_back(8383);
    tri3_xyze(0, 2) = 0.053004676773391018463;
    tri3_xyze(1, 2) = 0.021955255988737560485;
    tri3_xyze(2, 2) = 0.82673847315970405703;
    nids.push_back(-1197);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048657348057254377893;
    tri3_xyze(1, 0) = 0.028092332998904304131;
    tri3_xyze(2, 0) = 0.82408768370493490352;
    nids.push_back(8133);
    tri3_xyze(0, 1) = 0.046342037547421395527;
    tri3_xyze(1, 1) = 0.026755587852794311676;
    tri3_xyze(2, 1) = 0.81840622763408310547;
    nids.push_back(8135);
    tri3_xyze(0, 2) = 0.050239358035630420074;
    tri3_xyze(1, 2) = 0.020809823463272470595;
    tri3_xyze(2, 2) = 0.82124695566950922654;
    nids.push_back(-1198);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.046342037547421395527;
    tri3_xyze(1, 0) = 0.026755587852794311676;
    tri3_xyze(2, 0) = 0.81840622763408310547;
    nids.push_back(8135);
    tri3_xyze(0, 1) = 0.051687826609129038091;
    tri3_xyze(1, 1) = 0.013849711398434003465;
    tri3_xyze(2, 1) = 0.81840622763408354956;
    nids.push_back(8385);
    tri3_xyze(0, 2) = 0.050239358035630420074;
    tri3_xyze(1, 2) = 0.020809823463272470595;
    tri3_xyze(2, 2) = 0.82124695566950922654;
    nids.push_back(-1198);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.054270219928716861846;
    tri3_xyze(1, 0) = 0.014541661602957268312;
    tri3_xyze(2, 0) = 0.82408768370493534761;
    nids.push_back(8383);
    tri3_xyze(0, 1) = 0.048657348057254377893;
    tri3_xyze(1, 1) = 0.028092332998904304131;
    tri3_xyze(2, 1) = 0.82408768370493490352;
    nids.push_back(8133);
    tri3_xyze(0, 2) = 0.050239358035630420074;
    tri3_xyze(1, 2) = 0.020809823463272470595;
    tri3_xyze(2, 2) = 0.82124695566950922654;
    nids.push_back(-1198);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.046342037547421395527;
    tri3_xyze(1, 0) = 0.026755587852794311676;
    tri3_xyze(2, 0) = 0.81840622763408310547;
    nids.push_back(8135);
    tri3_xyze(0, 1) = 0.044661659217692294133;
    tri3_xyze(1, 1) = 0.025785420971784663025;
    tri3_xyze(2, 1) = 0.81243449435809200754;
    nids.push_back(8137);
    tri3_xyze(0, 2) = 0.048126282872754250208;
    tri3_xyze(1, 2) = 0.019934559072495215692;
    tri3_xyze(2, 2) = 0.81542036099608772304;
    nids.push_back(-1199);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.044661659217692294133;
    tri3_xyze(1, 0) = 0.025785420971784663025;
    tri3_xyze(2, 0) = 0.81243449435809200754;
    nids.push_back(8137);
    tri3_xyze(0, 1) = 0.049813608116774273082;
    tri3_xyze(1, 1) = 0.013347516066967877665;
    tri3_xyze(2, 1) = 0.8124344943580923406;
    nids.push_back(8387);
    tri3_xyze(0, 2) = 0.048126282872754250208;
    tri3_xyze(1, 2) = 0.019934559072495215692;
    tri3_xyze(2, 2) = 0.81542036099608772304;
    nids.push_back(-1199);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.049813608116774273082;
    tri3_xyze(1, 0) = 0.013347516066967877665;
    tri3_xyze(2, 0) = 0.8124344943580923406;
    nids.push_back(8387);
    tri3_xyze(0, 1) = 0.051687826609129038091;
    tri3_xyze(1, 1) = 0.013849711398434003465;
    tri3_xyze(2, 1) = 0.81840622763408354956;
    nids.push_back(8385);
    tri3_xyze(0, 2) = 0.048126282872754250208;
    tri3_xyze(1, 2) = 0.019934559072495215692;
    tri3_xyze(2, 2) = 0.81542036099608772304;
    nids.push_back(-1199);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051687826609129038091;
    tri3_xyze(1, 0) = 0.013849711398434003465;
    tri3_xyze(2, 0) = 0.81840622763408354956;
    nids.push_back(8385);
    tri3_xyze(0, 1) = 0.046342037547421395527;
    tri3_xyze(1, 1) = 0.026755587852794311676;
    tri3_xyze(2, 1) = 0.81840622763408310547;
    nids.push_back(8135);
    tri3_xyze(0, 2) = 0.048126282872754250208;
    tri3_xyze(1, 2) = 0.019934559072495215692;
    tri3_xyze(2, 2) = 0.81542036099608772304;
    nids.push_back(-1199);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.044661659217692294133;
    tri3_xyze(1, 0) = 0.025785420971784663025;
    tri3_xyze(2, 0) = 0.81243449435809200754;
    nids.push_back(8137);
    tri3_xyze(0, 1) = 0.043642713638136439169;
    tri3_xyze(1, 1) = 0.025197132467138269596;
    tri3_xyze(2, 1) = 0.80626666167806460006;
    nids.push_back(8139);
    tri3_xyze(0, 2) = 0.046698775742370961239;
    tri3_xyze(1, 2) = 0.019343266258705967692;
    tri3_xyze(2, 2) = 0.80935057801807841482;
    nids.push_back(-1200);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.043642713638136439169;
    tri3_xyze(1, 0) = 0.025197132467138269596;
    tri3_xyze(2, 0) = 0.80626666167806460006;
    nids.push_back(8139);
    tri3_xyze(0, 1) = 0.04867712199688085245;
    tri3_xyze(1, 1) = 0.013042995528933058746;
    tri3_xyze(2, 1) = 0.80626666167806482211;
    nids.push_back(8389);
    tri3_xyze(0, 2) = 0.046698775742370961239;
    tri3_xyze(1, 2) = 0.019343266258705967692;
    tri3_xyze(2, 2) = 0.80935057801807841482;
    nids.push_back(-1200);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.04867712199688085245;
    tri3_xyze(1, 0) = 0.013042995528933058746;
    tri3_xyze(2, 0) = 0.80626666167806482211;
    nids.push_back(8389);
    tri3_xyze(0, 1) = 0.049813608116774273082;
    tri3_xyze(1, 1) = 0.013347516066967877665;
    tri3_xyze(2, 1) = 0.8124344943580923406;
    nids.push_back(8387);
    tri3_xyze(0, 2) = 0.046698775742370961239;
    tri3_xyze(1, 2) = 0.019343266258705967692;
    tri3_xyze(2, 2) = 0.80935057801807841482;
    nids.push_back(-1200);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.049813608116774273082;
    tri3_xyze(1, 0) = 0.013347516066967877665;
    tri3_xyze(2, 0) = 0.8124344943580923406;
    nids.push_back(8387);
    tri3_xyze(0, 1) = 0.044661659217692294133;
    tri3_xyze(1, 1) = 0.025785420971784663025;
    tri3_xyze(2, 1) = 0.81243449435809200754;
    nids.push_back(8137);
    tri3_xyze(0, 2) = 0.046698775742370961239;
    tri3_xyze(1, 2) = 0.019343266258705967692;
    tri3_xyze(2, 2) = 0.80935057801807841482;
    nids.push_back(-1200);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.063395669874484861239;
    tri3_xyze(1, 0) = 0.10980452120245559611;
    tri3_xyze(2, 0) = 0.84221639627494737645;
    nids.push_back(7607);
    tri3_xyze(0, 1) = 0.060644482289136734943;
    tri3_xyze(1, 1) = 0.10503932452347999882;
    tri3_xyze(2, 1) = 0.84524135262314770145;
    nids.push_back(7609);
    tri3_xyze(0, 2) = 0.074864854408057332491;
    tri3_xyze(1, 2) = 0.097565777798631625473;
    tri3_xyze(2, 2) = 0.84372887444904787202;
    nids.push_back(-1085);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.060644482289136734943;
    tri3_xyze(1, 0) = 0.10503932452347999882;
    tri3_xyze(2, 0) = 0.84524135262314770145;
    nids.push_back(7609);
    tri3_xyze(0, 1) = 0.085764249336388140299;
    tri3_xyze(1, 1) = 0.085764249336379549948;
    tri3_xyze(2, 1) = 0.84524135262314836758;
    nids.push_back(7859);
    tri3_xyze(0, 2) = 0.074864854408057332491;
    tri3_xyze(1, 2) = 0.097565777798631625473;
    tri3_xyze(2, 2) = 0.84372887444904787202;
    nids.push_back(-1085);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003276614502322774e-14;
    tri3_xyze(1, 0) = 0.053511175705588047424;
    tri3_xyze(2, 0) = 0.81840622763408232831;
    nids.push_back(7135);
    tri3_xyze(0, 1) = 1.0003157803325775479e-14;
    tri3_xyze(1, 1) = 0.051570841943568840327;
    tri3_xyze(2, 1) = 0.81243449435809111936;
    nids.push_back(7137);
    tri3_xyze(0, 2) = 0.0067993068663601882035;
    tri3_xyze(1, 2) = 0.051645863093760305262;
    tri3_xyze(2, 2) = 0.81542036099608672384;
    nids.push_back(-999);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003157803325775479e-14;
    tri3_xyze(1, 0) = 0.051570841943568840327;
    tri3_xyze(2, 0) = 0.81243449435809111936;
    nids.push_back(7137);
    tri3_xyze(0, 1) = 0.013347516066977403032;
    tri3_xyze(1, 1) = 0.04981360811676469047;
    tri3_xyze(2, 1) = 0.81243449435809123038;
    nids.push_back(7387);
    tri3_xyze(0, 2) = 0.0067993068663601882035;
    tri3_xyze(1, 2) = 0.051645863093760305262;
    tri3_xyze(2, 2) = 0.81542036099608672384;
    nids.push_back(-999);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013347516066977403032;
    tri3_xyze(1, 0) = 0.04981360811676469047;
    tri3_xyze(2, 0) = 0.81243449435809123038;
    nids.push_back(7387);
    tri3_xyze(0, 1) = 0.013849711398443341481;
    tri3_xyze(1, 1) = 0.051687826609119649768;
    tri3_xyze(2, 1) = 0.81840622763408232831;
    nids.push_back(7385);
    tri3_xyze(0, 2) = 0.0067993068663601882035;
    tri3_xyze(1, 2) = 0.051645863093760305262;
    tri3_xyze(2, 2) = 0.81542036099608672384;
    nids.push_back(-999);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013849711398443341481;
    tri3_xyze(1, 0) = 0.051687826609119649768;
    tri3_xyze(2, 0) = 0.81840622763408232831;
    nids.push_back(7385);
    tri3_xyze(0, 1) = 1.0003276614502322774e-14;
    tri3_xyze(1, 1) = 0.053511175705588047424;
    tri3_xyze(2, 1) = 0.81840622763408232831;
    nids.push_back(7135);
    tri3_xyze(0, 2) = 0.0067993068663601882035;
    tri3_xyze(1, 2) = 0.051645863093760305262;
    tri3_xyze(2, 2) = 0.81542036099608672384;
    nids.push_back(-999);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003157803325775479e-14;
    tri3_xyze(1, 0) = 0.051570841943568840327;
    tri3_xyze(2, 0) = 0.81243449435809111936;
    nids.push_back(7137);
    tri3_xyze(0, 1) = 1.0003085758762357486e-14;
    tri3_xyze(1, 1) = 0.050394264934276310208;
    tri3_xyze(2, 1) = 0.80626666167806371188;
    nids.push_back(7139);
    tri3_xyze(0, 2) = 0.0065976278989850582757;
    tri3_xyze(1, 2) = 0.050113959247870223934;
    tri3_xyze(2, 2) = 0.80935057801807741562;
    nids.push_back(-1000);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003085758762357486e-14;
    tri3_xyze(1, 0) = 0.050394264934276310208;
    tri3_xyze(2, 0) = 0.80626666167806371188;
    nids.push_back(7139);
    tri3_xyze(0, 1) = 0.013042995528942823505;
    tri3_xyze(1, 1) = 0.048677121996871054732;
    tri3_xyze(2, 1) = 0.80626666167806371188;
    nids.push_back(7389);
    tri3_xyze(0, 2) = 0.0065976278989850582757;
    tri3_xyze(1, 2) = 0.050113959247870223934;
    tri3_xyze(2, 2) = 0.80935057801807741562;
    nids.push_back(-1000);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013042995528942823505;
    tri3_xyze(1, 0) = 0.048677121996871054732;
    tri3_xyze(2, 0) = 0.80626666167806371188;
    nids.push_back(7389);
    tri3_xyze(0, 1) = 0.013347516066977403032;
    tri3_xyze(1, 1) = 0.04981360811676469047;
    tri3_xyze(2, 1) = 0.81243449435809123038;
    nids.push_back(7387);
    tri3_xyze(0, 2) = 0.0065976278989850582757;
    tri3_xyze(1, 2) = 0.050113959247870223934;
    tri3_xyze(2, 2) = 0.80935057801807741562;
    nids.push_back(-1000);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013347516066977403032;
    tri3_xyze(1, 0) = 0.04981360811676469047;
    tri3_xyze(2, 0) = 0.81243449435809123038;
    nids.push_back(7387);
    tri3_xyze(0, 1) = 1.0003157803325775479e-14;
    tri3_xyze(1, 1) = 0.051570841943568840327;
    tri3_xyze(2, 1) = 0.81243449435809111936;
    nids.push_back(7137);
    tri3_xyze(0, 2) = 0.0065976278989850582757;
    tri3_xyze(1, 2) = 0.050113959247870223934;
    tri3_xyze(2, 2) = 0.80935057801807741562;
    nids.push_back(-1000);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013042995528942823505;
    tri3_xyze(1, 0) = 0.048677121996871054732;
    tri3_xyze(2, 0) = 0.80626666167806371188;
    nids.push_back(7389);
    tri3_xyze(0, 1) = 0.012940952255136016838;
    tri3_xyze(1, 1) = 0.048296291314453426713;
    tri3_xyze(2, 1) = 0.79999999999984849897;
    nids.push_back(7291);
    tri3_xyze(0, 2) = 0.019045270062809201844;
    tri3_xyze(1, 2) = 0.045979349284668277753;
    tri3_xyze(2, 2) = 0.80313333083895621645;
    nids.push_back(-1002);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.012940952255136016838;
    tri3_xyze(1, 0) = 0.048296291314453426713;
    tri3_xyze(2, 0) = 0.79999999999984849897;
    nids.push_back(7291);
    tri3_xyze(0, 1) = 0.025000000000009965639;
    tri3_xyze(1, 1) = 0.043301270189221953422;
    tri3_xyze(2, 1) = 0.79999999999984872101;
    nids.push_back(7541);
    tri3_xyze(0, 2) = 0.019045270062809201844;
    tri3_xyze(1, 2) = 0.045979349284668277753;
    tri3_xyze(2, 2) = 0.80313333083895621645;
    nids.push_back(-1002);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025000000000009965639;
    tri3_xyze(1, 0) = 0.043301270189221953422;
    tri3_xyze(2, 0) = 0.79999999999984872101;
    nids.push_back(7541);
    tri3_xyze(0, 1) = 0.025197132467147994456;
    tri3_xyze(1, 1) = 0.043642713638126676146;
    tri3_xyze(2, 1) = 0.80626666167806393393;
    nids.push_back(7639);
    tri3_xyze(0, 2) = 0.019045270062809201844;
    tri3_xyze(1, 2) = 0.045979349284668277753;
    tri3_xyze(2, 2) = 0.80313333083895621645;
    nids.push_back(-1002);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025197132467147994456;
    tri3_xyze(1, 0) = 0.043642713638126676146;
    tri3_xyze(2, 0) = 0.80626666167806393393;
    nids.push_back(7639);
    tri3_xyze(0, 1) = 0.013042995528942823505;
    tri3_xyze(1, 1) = 0.048677121996871054732;
    tri3_xyze(2, 1) = 0.80626666167806371188;
    nids.push_back(7389);
    tri3_xyze(0, 2) = 0.019045270062809201844;
    tri3_xyze(1, 2) = 0.045979349284668277753;
    tri3_xyze(2, 2) = 0.80313333083895621645;
    nids.push_back(-1002);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.04867712199688085245;
    tri3_xyze(1, 0) = 0.013042995528933058746;
    tri3_xyze(2, 0) = 0.80626666167806482211;
    nids.push_back(8389);
    tri3_xyze(0, 1) = 0.048296291314463390965;
    tri3_xyze(1, 1) = 0.012940952255126120241;
    tri3_xyze(2, 1) = 0.79999999999984960919;
    nids.push_back(8291);
    tri3_xyze(0, 2) = 0.04934191956141009;
    tri3_xyze(1, 2) = 0.0064959869460148441864;
    tri3_xyze(2, 2) = 0.80313333083895743769;
    nids.push_back(-652);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048296291314463390965;
    tri3_xyze(1, 0) = 0.012940952255126120241;
    tri3_xyze(2, 0) = 0.79999999999984960919;
    nids.push_back(8291);
    tri3_xyze(0, 1) = 0.050000000000010001722;
    tri3_xyze(1, 1) = 0;
    tri3_xyze(2, 1) = 0.79999999999985005328;
    nids.push_back(5541);
    tri3_xyze(0, 2) = 0.04934191956141009;
    tri3_xyze(1, 2) = 0.0064959869460148441864;
    tri3_xyze(2, 2) = 0.80313333083895743769;
    nids.push_back(-652);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.050000000000010001722;
    tri3_xyze(1, 0) = 0;
    tri3_xyze(2, 0) = 0.79999999999985005328;
    nids.push_back(5541);
    tri3_xyze(0, 1) = 0.050394264934286114865;
    tri3_xyze(1, 1) = 1.9687298290813598652e-16;
    tri3_xyze(2, 1) = 0.80626666167806526619;
    nids.push_back(5639);
    tri3_xyze(0, 2) = 0.04934191956141009;
    tri3_xyze(1, 2) = 0.0064959869460148441864;
    tri3_xyze(2, 2) = 0.80313333083895743769;
    nids.push_back(-652);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.050394264934286114865;
    tri3_xyze(1, 0) = 1.9687298290813598652e-16;
    tri3_xyze(2, 0) = 0.80626666167806526619;
    nids.push_back(5639);
    tri3_xyze(0, 1) = 0.04867712199688085245;
    tri3_xyze(1, 1) = 0.013042995528933058746;
    tri3_xyze(2, 1) = 0.80626666167806482211;
    nids.push_back(8389);
    tri3_xyze(0, 2) = 0.04934191956141009;
    tri3_xyze(1, 2) = 0.0064959869460148441864;
    tri3_xyze(2, 2) = 0.80313333083895743769;
    nids.push_back(-652);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051687826609129038091;
    tri3_xyze(1, 0) = 0.013849711398434003465;
    tri3_xyze(2, 0) = 0.81840622763408354956;
    nids.push_back(8385);
    tri3_xyze(0, 1) = 0.049813608116774273082;
    tri3_xyze(1, 1) = 0.013347516066967877665;
    tri3_xyze(2, 1) = 0.8124344943580923406;
    nids.push_back(8387);
    tri3_xyze(0, 2) = 0.051645863093769797669;
    tri3_xyze(1, 2) = 0.0067993068663507131438;
    tri3_xyze(2, 2) = 0.81542036099608816713;
    nids.push_back(-699);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.049813608116774273082;
    tri3_xyze(1, 0) = 0.013347516066967877665;
    tri3_xyze(2, 0) = 0.8124344943580923406;
    nids.push_back(8387);
    tri3_xyze(0, 1) = 0.051570841943578457633;
    tri3_xyze(1, 1) = 3.9064116126959097725e-16;
    tri3_xyze(2, 1) = 0.81243449435809278469;
    nids.push_back(5637);
    tri3_xyze(0, 2) = 0.051645863093769797669;
    tri3_xyze(1, 2) = 0.0067993068663507131438;
    tri3_xyze(2, 2) = 0.81542036099608816713;
    nids.push_back(-699);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.049813608116774273082;
    tri3_xyze(1, 0) = 0.013347516066967877665;
    tri3_xyze(2, 0) = 0.8124344943580923406;
    nids.push_back(8387);
    tri3_xyze(0, 1) = 0.04867712199688085245;
    tri3_xyze(1, 1) = 0.013042995528933058746;
    tri3_xyze(2, 1) = 0.80626666167806482211;
    nids.push_back(8389);
    tri3_xyze(0, 2) = 0.050113959247879924508;
    tri3_xyze(1, 2) = 0.0065976278989753819881;
    tri3_xyze(2, 2) = 0.80935057801807874789;
    nids.push_back(-700);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.04867712199688085245;
    tri3_xyze(1, 0) = 0.013042995528933058746;
    tri3_xyze(2, 0) = 0.80626666167806482211;
    nids.push_back(8389);
    tri3_xyze(0, 1) = 0.050394264934286114865;
    tri3_xyze(1, 1) = 1.9687298290813598652e-16;
    tri3_xyze(2, 1) = 0.80626666167806526619;
    nids.push_back(5639);
    tri3_xyze(0, 2) = 0.050113959247879924508;
    tri3_xyze(1, 2) = 0.0065976278989753819881;
    tri3_xyze(2, 2) = 0.80935057801807874789;
    nids.push_back(-700);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051570841943578457633;
    tri3_xyze(1, 0) = 3.9064116126959097725e-16;
    tri3_xyze(2, 0) = 0.81243449435809278469;
    nids.push_back(5637);
    tri3_xyze(0, 1) = 0.049813608116774273082;
    tri3_xyze(1, 1) = 0.013347516066967877665;
    tri3_xyze(2, 1) = 0.8124344943580923406;
    nids.push_back(8387);
    tri3_xyze(0, 2) = 0.050113959247879924508;
    tri3_xyze(1, 2) = 0.0065976278989753819881;
    tri3_xyze(2, 2) = 0.80935057801807874789;
    nids.push_back(-700);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.048174336836530634587;
    tri3_xyze(1, 0) = 0.048174336836541813145;
    tri3_xyze(2, 0) = 0.83852566213863799494;
    nids.push_back(6377);
    tri3_xyze(0, 1) = -0.051766340105089618528;
    tri3_xyze(1, 1) = 0.051766340105101005253;
    tri3_xyze(2, 1) = 0.84221639627494915281;
    nids.push_back(6375);
    tri3_xyze(0, 2) = -0.05558562711731944439;
    tri3_xyze(1, 2) = 0.042652351830863306781;
    tri3_xyze(2, 2) = 0.84037102920679385143;
    nids.push_back(-794);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.068490698010974626619;
    tri3_xyze(1, 0) = 0.068490698010986020283;
    tri3_xyze(2, 0) = 0.84990133642126153291;
    nids.push_back(6367);
    tri3_xyze(0, 1) = -0.04843023701175713791;
    tri3_xyze(1, 1) = 0.083883631126985089277;
    tri3_xyze(2, 1) = 0.84990133642126097779;
    nids.push_back(6617);
    tri3_xyze(0, 2) = -0.06035533905931733667;
    tri3_xyze(1, 2) = 0.078656609248550862779;
    tri3_xyze(2, 2) = 0.84990133642126108882;
    nids.push_back(-839);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.04843023701175713791;
    tri3_xyze(1, 0) = 0.083883631126985089277;
    tri3_xyze(2, 0) = 0.84990133642126097779;
    nids.push_back(6617);
    tri3_xyze(0, 1) = -0.051569762988222737909;
    tri3_xyze(1, 1) = 0.089321449629905819156;
    tri3_xyze(2, 1) = 0.84990133642126086677;
    nids.push_back(6615);
    tri3_xyze(0, 2) = -0.06035533905931733667;
    tri3_xyze(1, 2) = 0.078656609248550862779;
    tri3_xyze(2, 2) = 0.84990133642126108882;
    nids.push_back(-839);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.064085748208083950783;
    tri3_xyze(1, 0) = 0.064085748208095622003;
    tri3_xyze(2, 0) = 0.84911436253628247428;
    nids.push_back(6369);
    tri3_xyze(0, 1) = -0.045315467135346912397;
    tri3_xyze(1, 1) = 0.078488691447157338255;
    tri3_xyze(2, 1) = 0.84911436253628203019;
    nids.push_back(6619);
    tri3_xyze(0, 2) = -0.056580537591540649989;
    tri3_xyze(1, 2) = 0.073737192198306017454;
    tri3_xyze(2, 2) = 0.84950784947877178155;
    nids.push_back(-840);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.045315467135346912397;
    tri3_xyze(1, 0) = 0.078488691447157338255;
    tri3_xyze(2, 0) = 0.84911436253628203019;
    nids.push_back(6619);
    tri3_xyze(0, 1) = -0.04843023701175713791;
    tri3_xyze(1, 1) = 0.083883631126985089277;
    tri3_xyze(2, 1) = 0.84990133642126097779;
    nids.push_back(6617);
    tri3_xyze(0, 2) = -0.056580537591540649989;
    tri3_xyze(1, 2) = 0.073737192198306017454;
    tri3_xyze(2, 2) = 0.84950784947877178155;
    nids.push_back(-840);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.04843023701175713791;
    tri3_xyze(1, 0) = 0.083883631126985089277;
    tri3_xyze(2, 0) = 0.84990133642126097779;
    nids.push_back(6617);
    tri3_xyze(0, 1) = -0.068490698010974626619;
    tri3_xyze(1, 1) = 0.068490698010986020283;
    tri3_xyze(2, 1) = 0.84990133642126153291;
    nids.push_back(6367);
    tri3_xyze(0, 2) = -0.056580537591540649989;
    tri3_xyze(1, 2) = 0.073737192198306017454;
    tri3_xyze(2, 2) = 0.84950784947877178155;
    nids.push_back(-840);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.059785277507424225418;
    tri3_xyze(1, 0) = 0.059785277507435716227;
    tri3_xyze(2, 0) = 0.84755282581460589242;
    nids.push_back(6371);
    tri3_xyze(0, 1) = -0.042274575140616263424;
    tri3_xyze(1, 1) = 0.073221712011954515065;
    tri3_xyze(2, 1) = 0.84755282581460544833;
    nids.push_back(6621);
    tri3_xyze(0, 2) = -0.052865266997867838006;
    tri3_xyze(1, 2) = 0.068895357293660808295;
    tri3_xyze(2, 2) = 0.8483335941754439613;
    nids.push_back(-841);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.042274575140616263424;
    tri3_xyze(1, 0) = 0.073221712011954515065;
    tri3_xyze(2, 0) = 0.84755282581460544833;
    nids.push_back(6621);
    tri3_xyze(0, 1) = -0.045315467135346912397;
    tri3_xyze(1, 1) = 0.078488691447157338255;
    tri3_xyze(2, 1) = 0.84911436253628203019;
    nids.push_back(6619);
    tri3_xyze(0, 2) = -0.052865266997867838006;
    tri3_xyze(1, 2) = 0.068895357293660808295;
    tri3_xyze(2, 2) = 0.8483335941754439613;
    nids.push_back(-841);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.045315467135346912397;
    tri3_xyze(1, 0) = 0.078488691447157338255;
    tri3_xyze(2, 0) = 0.84911436253628203019;
    nids.push_back(6619);
    tri3_xyze(0, 1) = -0.064085748208083950783;
    tri3_xyze(1, 1) = 0.064085748208095622003;
    tri3_xyze(2, 1) = 0.84911436253628247428;
    nids.push_back(6369);
    tri3_xyze(0, 2) = -0.052865266997867838006;
    tri3_xyze(1, 2) = 0.068895357293660808295;
    tri3_xyze(2, 2) = 0.8483335941754439613;
    nids.push_back(-841);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.059785277507424225418;
    tri3_xyze(1, 0) = 0.059785277507435716227;
    tri3_xyze(2, 0) = 0.84755282581460589242;
    nids.push_back(6371);
    tri3_xyze(0, 1) = -0.055657106900921397841;
    tri3_xyze(1, 1) = 0.05565710690093281926;
    tri3_xyze(2, 1) = 0.84524135262314925576;
    nids.push_back(6373);
    tri3_xyze(0, 2) = -0.049268119314956268506;
    tri3_xyze(1, 2) = 0.064207463163433398101;
    tri3_xyze(2, 2) = 0.84639708921887735205;
    nids.push_back(-842);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.055657106900921397841;
    tri3_xyze(1, 0) = 0.05565710690093281926;
    tri3_xyze(2, 0) = 0.84524135262314925576;
    nids.push_back(6373);
    tri3_xyze(0, 1) = -0.039355517710863187342;
    tri3_xyze(1, 1) = 0.068165756233410507159;
    tri3_xyze(2, 1) = 0.8452413526231489227;
    nids.push_back(6623);
    tri3_xyze(0, 2) = -0.049268119314956268506;
    tri3_xyze(1, 2) = 0.064207463163433398101;
    tri3_xyze(2, 2) = 0.84639708921887735205;
    nids.push_back(-842);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.039355517710863187342;
    tri3_xyze(1, 0) = 0.068165756233410507159;
    tri3_xyze(2, 0) = 0.8452413526231489227;
    nids.push_back(6623);
    tri3_xyze(0, 1) = -0.042274575140616263424;
    tri3_xyze(1, 1) = 0.073221712011954515065;
    tri3_xyze(2, 1) = 0.84755282581460544833;
    nids.push_back(6621);
    tri3_xyze(0, 2) = -0.049268119314956268506;
    tri3_xyze(1, 2) = 0.064207463163433398101;
    tri3_xyze(2, 2) = 0.84639708921887735205;
    nids.push_back(-842);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.042274575140616263424;
    tri3_xyze(1, 0) = 0.073221712011954515065;
    tri3_xyze(2, 0) = 0.84755282581460544833;
    nids.push_back(6621);
    tri3_xyze(0, 1) = -0.059785277507424225418;
    tri3_xyze(1, 1) = 0.059785277507435716227;
    tri3_xyze(2, 1) = 0.84755282581460589242;
    nids.push_back(6371);
    tri3_xyze(0, 2) = -0.049268119314956268506;
    tri3_xyze(1, 2) = 0.064207463163433398101;
    tri3_xyze(2, 2) = 0.84639708921887735205;
    nids.push_back(-842);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.055657106900921397841;
    tri3_xyze(1, 0) = 0.05565710690093281926;
    tri3_xyze(2, 0) = 0.84524135262314925576;
    nids.push_back(6373);
    tri3_xyze(0, 1) = -0.051766340105089618528;
    tri3_xyze(1, 1) = 0.051766340105101005253;
    tri3_xyze(2, 1) = 0.84221639627494915281;
    nids.push_back(6375);
    tri3_xyze(0, 2) = -0.045845823710597304046;
    tri3_xyze(1, 2) = 0.05974744069846981559;
    tri3_xyze(2, 2) = 0.84372887444904909326;
    nids.push_back(-843);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.051766340105089618528;
    tri3_xyze(1, 0) = 0.051766340105101005253;
    tri3_xyze(2, 0) = 0.84221639627494915281;
    nids.push_back(6375);
    tri3_xyze(0, 1) = -0.036604330125515026351;
    tri3_xyze(1, 1) = 0.063400559554434951504;
    tri3_xyze(2, 1) = 0.84221639627494881974;
    nids.push_back(6625);
    tri3_xyze(0, 2) = -0.045845823710597304046;
    tri3_xyze(1, 2) = 0.05974744069846981559;
    tri3_xyze(2, 2) = 0.84372887444904909326;
    nids.push_back(-843);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.036604330125515026351;
    tri3_xyze(1, 0) = 0.063400559554434951504;
    tri3_xyze(2, 0) = 0.84221639627494881974;
    nids.push_back(6625);
    tri3_xyze(0, 1) = -0.039355517710863187342;
    tri3_xyze(1, 1) = 0.068165756233410507159;
    tri3_xyze(2, 1) = 0.8452413526231489227;
    nids.push_back(6623);
    tri3_xyze(0, 2) = -0.045845823710597304046;
    tri3_xyze(1, 2) = 0.05974744069846981559;
    tri3_xyze(2, 2) = 0.84372887444904909326;
    nids.push_back(-843);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.039355517710863187342;
    tri3_xyze(1, 0) = 0.068165756233410507159;
    tri3_xyze(2, 0) = 0.8452413526231489227;
    nids.push_back(6623);
    tri3_xyze(0, 1) = -0.055657106900921397841;
    tri3_xyze(1, 1) = 0.05565710690093281926;
    tri3_xyze(2, 1) = 0.84524135262314925576;
    nids.push_back(6373);
    tri3_xyze(0, 2) = -0.045845823710597304046;
    tri3_xyze(1, 2) = 0.05974744069846981559;
    tri3_xyze(2, 2) = 0.84372887444904909326;
    nids.push_back(-843);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.051766340105089618528;
    tri3_xyze(1, 0) = 0.051766340105101005253;
    tri3_xyze(2, 0) = 0.84221639627494915281;
    nids.push_back(6375);
    tri3_xyze(0, 1) = -0.048174336836530634587;
    tri3_xyze(1, 1) = 0.048174336836541813145;
    tri3_xyze(2, 1) = 0.83852566213863799494;
    nids.push_back(6377);
    tri3_xyze(0, 2) = -0.0426523518308520172;
    tri3_xyze(1, 2) = 0.055585627117330754787;
    tri3_xyze(2, 2) = 0.84037102920679351836;
    nids.push_back(-844);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.048174336836530634587;
    tri3_xyze(1, 0) = 0.048174336836541813145;
    tri3_xyze(2, 0) = 0.83852566213863799494;
    nids.push_back(6377);
    tri3_xyze(0, 1) = -0.034064400256272782397;
    tri3_xyze(1, 1) = 0.05900127197324522843;
    tri3_xyze(2, 1) = 0.83852566213863766187;
    nids.push_back(6627);
    tri3_xyze(0, 2) = -0.0426523518308520172;
    tri3_xyze(1, 2) = 0.055585627117330754787;
    tri3_xyze(2, 2) = 0.84037102920679351836;
    nids.push_back(-844);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.034064400256272782397;
    tri3_xyze(1, 0) = 0.05900127197324522843;
    tri3_xyze(2, 0) = 0.83852566213863766187;
    nids.push_back(6627);
    tri3_xyze(0, 1) = -0.036604330125515026351;
    tri3_xyze(1, 1) = 0.063400559554434951504;
    tri3_xyze(2, 1) = 0.84221639627494881974;
    nids.push_back(6625);
    tri3_xyze(0, 2) = -0.0426523518308520172;
    tri3_xyze(1, 2) = 0.055585627117330754787;
    tri3_xyze(2, 2) = 0.84037102920679351836;
    nids.push_back(-844);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.036604330125515026351;
    tri3_xyze(1, 0) = 0.063400559554434951504;
    tri3_xyze(2, 0) = 0.84221639627494881974;
    nids.push_back(6625);
    tri3_xyze(0, 1) = -0.051766340105089618528;
    tri3_xyze(1, 1) = 0.051766340105101005253;
    tri3_xyze(2, 1) = 0.84221639627494915281;
    nids.push_back(6375);
    tri3_xyze(0, 2) = -0.0426523518308520172;
    tri3_xyze(1, 2) = 0.055585627117330754787;
    tri3_xyze(2, 2) = 0.84037102920679351836;
    nids.push_back(-844);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.048174336836530634587;
    tri3_xyze(1, 0) = 0.048174336836541813145;
    tri3_xyze(2, 0) = 0.83852566213863799494;
    nids.push_back(6377);
    tri3_xyze(0, 1) = -0.044937745132548299887;
    tri3_xyze(1, 1) = 0.044937745132559284156;
    tri3_xyze(2, 1) = 0.83422735529628300188;
    nids.push_back(6379);
    tri3_xyze(0, 2) = -0.039738066634951593747;
    tri3_xyze(1, 2) = 0.051787656706337614754;
    tri3_xyze(2, 2) = 0.83637650871746038739;
    nids.push_back(-845);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.044937745132548299887;
    tri3_xyze(1, 0) = 0.044937745132559284156;
    tri3_xyze(2, 0) = 0.83422735529628300188;
    nids.push_back(6379);
    tri3_xyze(0, 1) = -0.031775784314454665058;
    tri3_xyze(1, 1) = 0.055037272883004147161;
    tri3_xyze(2, 1) = 0.83422735529628266882;
    nids.push_back(6629);
    tri3_xyze(0, 2) = -0.039738066634951593747;
    tri3_xyze(1, 2) = 0.051787656706337614754;
    tri3_xyze(2, 2) = 0.83637650871746038739;
    nids.push_back(-845);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.031775784314454665058;
    tri3_xyze(1, 0) = 0.055037272883004147161;
    tri3_xyze(2, 0) = 0.83422735529628266882;
    nids.push_back(6629);
    tri3_xyze(0, 1) = -0.034064400256272782397;
    tri3_xyze(1, 1) = 0.05900127197324522843;
    tri3_xyze(2, 1) = 0.83852566213863766187;
    nids.push_back(6627);
    tri3_xyze(0, 2) = -0.039738066634951593747;
    tri3_xyze(1, 2) = 0.051787656706337614754;
    tri3_xyze(2, 2) = 0.83637650871746038739;
    nids.push_back(-845);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.034064400256272782397;
    tri3_xyze(1, 0) = 0.05900127197324522843;
    tri3_xyze(2, 0) = 0.83852566213863766187;
    nids.push_back(6627);
    tri3_xyze(0, 1) = -0.048174336836530634587;
    tri3_xyze(1, 1) = 0.048174336836541813145;
    tri3_xyze(2, 1) = 0.83852566213863799494;
    nids.push_back(6377);
    tri3_xyze(0, 2) = -0.039738066634951593747;
    tri3_xyze(1, 2) = 0.051787656706337614754;
    tri3_xyze(2, 2) = 0.83637650871746038739;
    nids.push_back(-845);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.042107607977760558304;
    tri3_xyze(1, 0) = 0.042107607977771480123;
    tri3_xyze(2, 0) = 0.82938926261447243338;
    nids.push_back(6381);
    tri3_xyze(0, 1) = -0.029774575140616293956;
    tri3_xyze(1, 1) = 0.051571076917342986712;
    tri3_xyze(2, 1) = 0.82938926261447210031;
    nids.push_back(6631);
    tri3_xyze(0, 2) = -0.037148928141344954301;
    tri3_xyze(1, 2) = 0.048413425727669476273;
    tri3_xyze(2, 2) = 0.83180830895537760661;
    nids.push_back(-846);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.029774575140616293956;
    tri3_xyze(1, 0) = 0.051571076917342986712;
    tri3_xyze(2, 0) = 0.82938926261447210031;
    nids.push_back(6631);
    tri3_xyze(0, 1) = -0.031775784314454665058;
    tri3_xyze(1, 1) = 0.055037272883004147161;
    tri3_xyze(2, 1) = 0.83422735529628266882;
    nids.push_back(6629);
    tri3_xyze(0, 2) = -0.037148928141344954301;
    tri3_xyze(1, 2) = 0.048413425727669476273;
    tri3_xyze(2, 2) = 0.83180830895537760661;
    nids.push_back(-846);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.031775784314454665058;
    tri3_xyze(1, 0) = 0.055037272883004147161;
    tri3_xyze(2, 0) = 0.83422735529628266882;
    nids.push_back(6629);
    tri3_xyze(0, 1) = -0.044937745132548299887;
    tri3_xyze(1, 1) = 0.044937745132559284156;
    tri3_xyze(2, 1) = 0.83422735529628300188;
    nids.push_back(6379);
    tri3_xyze(0, 2) = -0.037148928141344954301;
    tri3_xyze(1, 2) = 0.048413425727669476273;
    tri3_xyze(2, 2) = 0.83180830895537760661;
    nids.push_back(-846);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.028092332998893382312;
    tri3_xyze(1, 0) = 0.048657348057245253248;
    tri3_xyze(2, 0) = 0.82408768370493423738;
    nids.push_back(6633);
    tri3_xyze(0, 1) = -0.029774575140616293956;
    tri3_xyze(1, 1) = 0.051571076917342986712;
    tri3_xyze(2, 1) = 0.82938926261447210031;
    nids.push_back(6631);
    tri3_xyze(0, 2) = -0.034925768610752673526;
    tri3_xyze(1, 2) = 0.045516147819527730256;
    tri3_xyze(2, 2) = 0.82673847315970339089;
    nids.push_back(-847);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.029774575140616293956;
    tri3_xyze(1, 0) = 0.051571076917342986712;
    tri3_xyze(2, 0) = 0.82938926261447210031;
    nids.push_back(6631);
    tri3_xyze(0, 1) = -0.042107607977760558304;
    tri3_xyze(1, 1) = 0.042107607977771480123;
    tri3_xyze(2, 1) = 0.82938926261447243338;
    nids.push_back(6381);
    tri3_xyze(0, 2) = -0.034925768610752673526;
    tri3_xyze(1, 2) = 0.045516147819527730256;
    tri3_xyze(2, 2) = 0.82673847315970339089;
    nids.push_back(-847);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.057725424859363591579;
    tri3_xyze(1, 0) = 0.099983368744936226835;
    tri3_xyze(2, 0) = 0.84755282581460467117;
    nids.push_back(6611);
    tri3_xyze(0, 1) = -0.054684532864633067506;
    tri3_xyze(1, 1) = 0.094716389309733514668;
    tri3_xyze(2, 1) = 0.84911436253628147508;
    nids.push_back(6613);
    tri3_xyze(0, 2) = -0.042649408390066499697;
    tri3_xyze(1, 2) = 0.10296478016251330023;
    tri3_xyze(2, 2) = 0.84833359417544285108;
    nids.push_back(-887);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.054684532864633067506;
    tri3_xyze(1, 0) = 0.094716389309733514668;
    tri3_xyze(2, 0) = 0.84911436253628147508;
    nids.push_back(6613);
    tri3_xyze(0, 1) = -0.028306797155798663762;
    tri3_xyze(1, 1) = 0.10564240518502582278;
    tri3_xyze(2, 1) = 0.84911436253628114201;
    nids.push_back(6863);
    tri3_xyze(0, 2) = -0.042649408390066499697;
    tri3_xyze(1, 2) = 0.10296478016251330023;
    tri3_xyze(2, 2) = 0.84833359417544285108;
    nids.push_back(-887);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.054684532864633067506;
    tri3_xyze(1, 0) = 0.094716389309733514668;
    tri3_xyze(2, 0) = 0.84911436253628147508;
    nids.push_back(6613);
    tri3_xyze(0, 1) = -0.051569762988222737909;
    tri3_xyze(1, 1) = 0.089321449629905819156;
    tri3_xyze(2, 1) = 0.84990133642126086677;
    nids.push_back(6615);
    tri3_xyze(0, 2) = -0.04031389165854998291;
    tri3_xyze(1, 2) = 0.0973263439941367825;
    tri3_xyze(2, 2) = 0.84950784947877100439;
    nids.push_back(-888);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.051569762988222737909;
    tri3_xyze(1, 0) = 0.089321449629905819156;
    tri3_xyze(2, 0) = 0.84990133642126086677;
    nids.push_back(6615);
    tri3_xyze(0, 1) = -0.026694473625545441647;
    tri3_xyze(1, 1) = 0.099625131851881987277;
    tri3_xyze(2, 1) = 0.8499013364212605337;
    nids.push_back(6865);
    tri3_xyze(0, 2) = -0.04031389165854998291;
    tri3_xyze(1, 2) = 0.0973263439941367825;
    tri3_xyze(2, 2) = 0.84950784947877100439;
    nids.push_back(-888);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.026694473625545441647;
    tri3_xyze(1, 0) = 0.099625131851881987277;
    tri3_xyze(2, 0) = 0.8499013364212605337;
    nids.push_back(6865);
    tri3_xyze(0, 1) = -0.028306797155798663762;
    tri3_xyze(1, 1) = 0.10564240518502582278;
    tri3_xyze(2, 1) = 0.84911436253628114201;
    nids.push_back(6863);
    tri3_xyze(0, 2) = -0.04031389165854998291;
    tri3_xyze(1, 2) = 0.0973263439941367825;
    tri3_xyze(2, 2) = 0.84950784947877100439;
    nids.push_back(-888);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.028306797155798663762;
    tri3_xyze(1, 0) = 0.10564240518502582278;
    tri3_xyze(2, 0) = 0.84911436253628114201;
    nids.push_back(6863);
    tri3_xyze(0, 1) = -0.054684532864633067506;
    tri3_xyze(1, 1) = 0.094716389309733514668;
    tri3_xyze(2, 1) = 0.84911436253628147508;
    nids.push_back(6613);
    tri3_xyze(0, 2) = -0.04031389165854998291;
    tri3_xyze(1, 2) = 0.0973263439941367825;
    tri3_xyze(2, 2) = 0.84950784947877100439;
    nids.push_back(-888);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.051569762988222737909;
    tri3_xyze(1, 0) = 0.089321449629905819156;
    tri3_xyze(2, 0) = 0.84990133642126086677;
    nids.push_back(6615);
    tri3_xyze(0, 1) = -0.04843023701175713791;
    tri3_xyze(1, 1) = 0.083883631126985089277;
    tri3_xyze(2, 1) = 0.84990133642126097779;
    nids.push_back(6617);
    tri3_xyze(0, 2) = -0.037940952255115990843;
    tri3_xyze(1, 2) = 0.091597561503676927508;
    tri3_xyze(2, 2) = 0.84990133642126075575;
    nids.push_back(-889);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.04843023701175713791;
    tri3_xyze(1, 0) = 0.083883631126985089277;
    tri3_xyze(2, 0) = 0.84990133642126097779;
    nids.push_back(6617);
    tri3_xyze(0, 1) = -0.025069335394938645906;
    tri3_xyze(1, 1) = 0.093560033405934842077;
    tri3_xyze(2, 1) = 0.84990133642126064473;
    nids.push_back(6867);
    tri3_xyze(0, 2) = -0.037940952255115990843;
    tri3_xyze(1, 2) = 0.091597561503676927508;
    tri3_xyze(2, 2) = 0.84990133642126075575;
    nids.push_back(-889);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.025069335394938645906;
    tri3_xyze(1, 0) = 0.093560033405934842077;
    tri3_xyze(2, 0) = 0.84990133642126064473;
    nids.push_back(6867);
    tri3_xyze(0, 1) = -0.026694473625545441647;
    tri3_xyze(1, 1) = 0.099625131851881987277;
    tri3_xyze(2, 1) = 0.8499013364212605337;
    nids.push_back(6865);
    tri3_xyze(0, 2) = -0.037940952255115990843;
    tri3_xyze(1, 2) = 0.091597561503676927508;
    tri3_xyze(2, 2) = 0.84990133642126075575;
    nids.push_back(-889);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.026694473625545441647;
    tri3_xyze(1, 0) = 0.099625131851881987277;
    tri3_xyze(2, 0) = 0.8499013364212605337;
    nids.push_back(6865);
    tri3_xyze(0, 1) = -0.051569762988222737909;
    tri3_xyze(1, 1) = 0.089321449629905819156;
    tri3_xyze(2, 1) = 0.84990133642126086677;
    nids.push_back(6615);
    tri3_xyze(0, 2) = -0.037940952255115990843;
    tri3_xyze(1, 2) = 0.091597561503676927508;
    tri3_xyze(2, 2) = 0.84990133642126075575;
    nids.push_back(-889);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.04843023701175713791;
    tri3_xyze(1, 0) = 0.083883631126985089277;
    tri3_xyze(2, 0) = 0.84990133642126097779;
    nids.push_back(6617);
    tri3_xyze(0, 1) = -0.045315467135346912397;
    tri3_xyze(1, 1) = 0.078488691447157338255;
    tri3_xyze(2, 1) = 0.84911436253628203019;
    nids.push_back(6619);
    tri3_xyze(0, 2) = -0.035568012851682047348;
    tri3_xyze(1, 2) = 0.08586877901321704476;
    tri3_xyze(2, 2) = 0.84950784947877133746;
    nids.push_back(-890);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.045315467135346912397;
    tri3_xyze(1, 0) = 0.078488691447157338255;
    tri3_xyze(2, 0) = 0.84911436253628203019;
    nids.push_back(6619);
    tri3_xyze(0, 1) = -0.02345701186468549318;
    tri3_xyze(1, 1) = 0.087542760072790937187;
    tri3_xyze(2, 1) = 0.84911436253628169712;
    nids.push_back(6869);
    tri3_xyze(0, 2) = -0.035568012851682047348;
    tri3_xyze(1, 2) = 0.08586877901321704476;
    tri3_xyze(2, 2) = 0.84950784947877133746;
    nids.push_back(-890);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.02345701186468549318;
    tri3_xyze(1, 0) = 0.087542760072790937187;
    tri3_xyze(2, 0) = 0.84911436253628169712;
    nids.push_back(6869);
    tri3_xyze(0, 1) = -0.025069335394938645906;
    tri3_xyze(1, 1) = 0.093560033405934842077;
    tri3_xyze(2, 1) = 0.84990133642126064473;
    nids.push_back(6867);
    tri3_xyze(0, 2) = -0.035568012851682047348;
    tri3_xyze(1, 2) = 0.08586877901321704476;
    tri3_xyze(2, 2) = 0.84950784947877133746;
    nids.push_back(-890);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.025069335394938645906;
    tri3_xyze(1, 0) = 0.093560033405934842077;
    tri3_xyze(2, 0) = 0.84990133642126064473;
    nids.push_back(6867);
    tri3_xyze(0, 1) = -0.04843023701175713791;
    tri3_xyze(1, 1) = 0.083883631126985089277;
    tri3_xyze(2, 1) = 0.84990133642126097779;
    nids.push_back(6617);
    tri3_xyze(0, 2) = -0.035568012851682047348;
    tri3_xyze(1, 2) = 0.08586877901321704476;
    tri3_xyze(2, 2) = 0.84950784947877133746;
    nids.push_back(-890);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.045315467135346912397;
    tri3_xyze(1, 0) = 0.078488691447157338255;
    tri3_xyze(2, 0) = 0.84911436253628203019;
    nids.push_back(6619);
    tri3_xyze(0, 1) = -0.042274575140616263424;
    tri3_xyze(1, 1) = 0.073221712011954515065;
    tri3_xyze(2, 1) = 0.84755282581460544833;
    nids.push_back(6621);
    tri3_xyze(0, 2) = -0.033232496120165509745;
    tri3_xyze(1, 2) = 0.080230342844840457639;
    tri3_xyze(2, 2) = 0.84833359417544362824;
    nids.push_back(-891);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.042274575140616263424;
    tri3_xyze(1, 0) = 0.073221712011954515065;
    tri3_xyze(2, 0) = 0.84755282581460544833;
    nids.push_back(6621);
    tri3_xyze(0, 1) = -0.02188293034001336998;
    tri3_xyze(1, 1) = 0.081668207847459053927;
    tri3_xyze(2, 1) = 0.84755282581460522628;
    nids.push_back(6871);
    tri3_xyze(0, 2) = -0.033232496120165509745;
    tri3_xyze(1, 2) = 0.080230342844840457639;
    tri3_xyze(2, 2) = 0.84833359417544362824;
    nids.push_back(-891);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.02188293034001336998;
    tri3_xyze(1, 0) = 0.081668207847459053927;
    tri3_xyze(2, 0) = 0.84755282581460522628;
    nids.push_back(6871);
    tri3_xyze(0, 1) = -0.02345701186468549318;
    tri3_xyze(1, 1) = 0.087542760072790937187;
    tri3_xyze(2, 1) = 0.84911436253628169712;
    nids.push_back(6869);
    tri3_xyze(0, 2) = -0.033232496120165509745;
    tri3_xyze(1, 2) = 0.080230342844840457639;
    tri3_xyze(2, 2) = 0.84833359417544362824;
    nids.push_back(-891);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.02345701186468549318;
    tri3_xyze(1, 0) = 0.087542760072790937187;
    tri3_xyze(2, 0) = 0.84911436253628169712;
    nids.push_back(6869);
    tri3_xyze(0, 1) = -0.045315467135346912397;
    tri3_xyze(1, 1) = 0.078488691447157338255;
    tri3_xyze(2, 1) = 0.84911436253628203019;
    nids.push_back(6619);
    tri3_xyze(0, 2) = -0.033232496120165509745;
    tri3_xyze(1, 2) = 0.080230342844840457639;
    tri3_xyze(2, 2) = 0.84833359417544362824;
    nids.push_back(-891);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.042274575140616263424;
    tri3_xyze(1, 0) = 0.073221712011954515065;
    tri3_xyze(2, 0) = 0.84755282581460544833;
    nids.push_back(6621);
    tri3_xyze(0, 1) = -0.039355517710863187342;
    tri3_xyze(1, 1) = 0.068165756233410507159;
    tri3_xyze(2, 1) = 0.8452413526231489227;
    nids.push_back(6623);
    tri3_xyze(0, 2) = -0.030971234554592469435;
    tri3_xyze(1, 2) = 0.074771174505160992174;
    tri3_xyze(2, 2) = 0.84639708921887701898;
    nids.push_back(-892);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.039355517710863187342;
    tri3_xyze(1, 0) = 0.068165756233410507159;
    tri3_xyze(2, 0) = 0.8452413526231489227;
    nids.push_back(6623);
    tri3_xyze(0, 1) = -0.020371915026877046584;
    tri3_xyze(1, 1) = 0.076029021927819878668;
    tri3_xyze(2, 1) = 0.84524135262314858963;
    nids.push_back(6873);
    tri3_xyze(0, 2) = -0.030971234554592469435;
    tri3_xyze(1, 2) = 0.074771174505160992174;
    tri3_xyze(2, 2) = 0.84639708921887701898;
    nids.push_back(-892);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.020371915026877046584;
    tri3_xyze(1, 0) = 0.076029021927819878668;
    tri3_xyze(2, 0) = 0.84524135262314858963;
    nids.push_back(6873);
    tri3_xyze(0, 1) = -0.02188293034001336998;
    tri3_xyze(1, 1) = 0.081668207847459053927;
    tri3_xyze(2, 1) = 0.84755282581460522628;
    nids.push_back(6871);
    tri3_xyze(0, 2) = -0.030971234554592469435;
    tri3_xyze(1, 2) = 0.074771174505160992174;
    tri3_xyze(2, 2) = 0.84639708921887701898;
    nids.push_back(-892);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.02188293034001336998;
    tri3_xyze(1, 0) = 0.081668207847459053927;
    tri3_xyze(2, 0) = 0.84755282581460522628;
    nids.push_back(6871);
    tri3_xyze(0, 1) = -0.042274575140616263424;
    tri3_xyze(1, 1) = 0.073221712011954515065;
    tri3_xyze(2, 1) = 0.84755282581460544833;
    nids.push_back(6621);
    tri3_xyze(0, 2) = -0.030971234554592469435;
    tri3_xyze(1, 2) = 0.074771174505160992174;
    tri3_xyze(2, 2) = 0.84639708921887701898;
    nids.push_back(-892);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.039355517710863187342;
    tri3_xyze(1, 0) = 0.068165756233410507159;
    tri3_xyze(2, 0) = 0.8452413526231489227;
    nids.push_back(6623);
    tri3_xyze(0, 1) = -0.036604330125515026351;
    tri3_xyze(1, 1) = 0.063400559554434951504;
    tri3_xyze(2, 1) = 0.84221639627494881974;
    nids.push_back(6625);
    tri3_xyze(0, 2) = -0.028819889600664228002;
    tri3_xyze(1, 2) = 0.069577368340044493111;
    tri3_xyze(2, 2) = 0.8437288744490487602;
    nids.push_back(-893);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.036604330125515026351;
    tri3_xyze(1, 0) = 0.063400559554434951504;
    tri3_xyze(2, 0) = 0.84221639627494881974;
    nids.push_back(6625);
    tri3_xyze(0, 1) = -0.018947795539401655202;
    tri3_xyze(1, 1) = 0.070714135644512676748;
    tri3_xyze(2, 1) = 0.8422163962749485977;
    nids.push_back(6875);
    tri3_xyze(0, 2) = -0.028819889600664228002;
    tri3_xyze(1, 2) = 0.069577368340044493111;
    tri3_xyze(2, 2) = 0.8437288744490487602;
    nids.push_back(-893);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.018947795539401655202;
    tri3_xyze(1, 0) = 0.070714135644512676748;
    tri3_xyze(2, 0) = 0.8422163962749485977;
    nids.push_back(6875);
    tri3_xyze(0, 1) = -0.020371915026877046584;
    tri3_xyze(1, 1) = 0.076029021927819878668;
    tri3_xyze(2, 1) = 0.84524135262314858963;
    nids.push_back(6873);
    tri3_xyze(0, 2) = -0.028819889600664228002;
    tri3_xyze(1, 2) = 0.069577368340044493111;
    tri3_xyze(2, 2) = 0.8437288744490487602;
    nids.push_back(-893);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.020371915026877046584;
    tri3_xyze(1, 0) = 0.076029021927819878668;
    tri3_xyze(2, 0) = 0.84524135262314858963;
    nids.push_back(6873);
    tri3_xyze(0, 1) = -0.039355517710863187342;
    tri3_xyze(1, 1) = 0.068165756233410507159;
    tri3_xyze(2, 1) = 0.8452413526231489227;
    nids.push_back(6623);
    tri3_xyze(0, 2) = -0.028819889600664228002;
    tri3_xyze(1, 2) = 0.069577368340044493111;
    tri3_xyze(2, 2) = 0.8437288744490487602;
    nids.push_back(-893);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.036604330125515026351;
    tri3_xyze(1, 0) = 0.063400559554434951504;
    tri3_xyze(2, 0) = 0.84221639627494881974;
    nids.push_back(6625);
    tri3_xyze(0, 1) = -0.034064400256272782397;
    tri3_xyze(1, 1) = 0.05900127197324522843;
    tri3_xyze(2, 1) = 0.83852566213863766187;
    nids.push_back(6627);
    tri3_xyze(0, 2) = -0.026812389253455446575;
    tri3_xyze(1, 2) = 0.064730833775344262526;
    tri3_xyze(2, 2) = 0.84037102920679307427;
    nids.push_back(-894);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.034064400256272782397;
    tri3_xyze(1, 0) = 0.05900127197324522843;
    tri3_xyze(2, 0) = 0.83852566213863766187;
    nids.push_back(6627);
    tri3_xyze(0, 1) = -0.017633031092632322351;
    tri3_xyze(1, 1) = 0.0658073679291842073;
    tri3_xyze(2, 1) = 0.83852566213863743982;
    nids.push_back(6877);
    tri3_xyze(0, 2) = -0.026812389253455446575;
    tri3_xyze(1, 2) = 0.064730833775344262526;
    tri3_xyze(2, 2) = 0.84037102920679307427;
    nids.push_back(-894);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.017633031092632322351;
    tri3_xyze(1, 0) = 0.0658073679291842073;
    tri3_xyze(2, 0) = 0.83852566213863743982;
    nids.push_back(6877);
    tri3_xyze(0, 1) = -0.018947795539401655202;
    tri3_xyze(1, 1) = 0.070714135644512676748;
    tri3_xyze(2, 1) = 0.8422163962749485977;
    nids.push_back(6875);
    tri3_xyze(0, 2) = -0.026812389253455446575;
    tri3_xyze(1, 2) = 0.064730833775344262526;
    tri3_xyze(2, 2) = 0.84037102920679307427;
    nids.push_back(-894);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.018947795539401655202;
    tri3_xyze(1, 0) = 0.070714135644512676748;
    tri3_xyze(2, 0) = 0.8422163962749485977;
    nids.push_back(6875);
    tri3_xyze(0, 1) = -0.036604330125515026351;
    tri3_xyze(1, 1) = 0.063400559554434951504;
    tri3_xyze(2, 1) = 0.84221639627494881974;
    nids.push_back(6625);
    tri3_xyze(0, 2) = -0.026812389253455446575;
    tri3_xyze(1, 2) = 0.064730833775344262526;
    tri3_xyze(2, 2) = 0.84037102920679307427;
    nids.push_back(-894);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.034064400256272782397;
    tri3_xyze(1, 0) = 0.05900127197324522843;
    tri3_xyze(2, 0) = 0.83852566213863766187;
    nids.push_back(6627);
    tri3_xyze(0, 1) = -0.031775784314454665058;
    tri3_xyze(1, 1) = 0.055037272883004147161;
    tri3_xyze(2, 1) = 0.83422735529628266882;
    nids.push_back(6629);
    tri3_xyze(0, 2) = -0.024980392992664145896;
    tri3_xyze(1, 2) = 0.060308003556324936423;
    tri3_xyze(2, 2) = 0.83637650871746005432;
    nids.push_back(-895);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.031775784314454665058;
    tri3_xyze(1, 0) = 0.055037272883004147161;
    tri3_xyze(2, 0) = 0.83422735529628266882;
    nids.push_back(6629);
    tri3_xyze(0, 1) = -0.016448356307296806839;
    tri3_xyze(1, 1) = 0.061386101439866169738;
    tri3_xyze(2, 1) = 0.83422735529628255779;
    nids.push_back(6879);
    tri3_xyze(0, 2) = -0.024980392992664145896;
    tri3_xyze(1, 2) = 0.060308003556324936423;
    tri3_xyze(2, 2) = 0.83637650871746005432;
    nids.push_back(-895);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.016448356307296806839;
    tri3_xyze(1, 0) = 0.061386101439866169738;
    tri3_xyze(2, 0) = 0.83422735529628255779;
    nids.push_back(6879);
    tri3_xyze(0, 1) = -0.017633031092632322351;
    tri3_xyze(1, 1) = 0.0658073679291842073;
    tri3_xyze(2, 1) = 0.83852566213863743982;
    nids.push_back(6877);
    tri3_xyze(0, 2) = -0.024980392992664145896;
    tri3_xyze(1, 2) = 0.060308003556324936423;
    tri3_xyze(2, 2) = 0.83637650871746005432;
    nids.push_back(-895);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.017633031092632322351;
    tri3_xyze(1, 0) = 0.0658073679291842073;
    tri3_xyze(2, 0) = 0.83852566213863743982;
    nids.push_back(6877);
    tri3_xyze(0, 1) = -0.034064400256272782397;
    tri3_xyze(1, 1) = 0.05900127197324522843;
    tri3_xyze(2, 1) = 0.83852566213863766187;
    nids.push_back(6627);
    tri3_xyze(0, 2) = -0.024980392992664145896;
    tri3_xyze(1, 2) = 0.060308003556324936423;
    tri3_xyze(2, 2) = 0.83637650871746005432;
    nids.push_back(-895);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.031775784314454665058;
    tri3_xyze(1, 0) = 0.055037272883004147161;
    tri3_xyze(2, 0) = 0.83422735529628266882;
    nids.push_back(6629);
    tri3_xyze(0, 1) = -0.029774575140616293956;
    tri3_xyze(1, 1) = 0.051571076917342986712;
    tri3_xyze(2, 1) = 0.82938926261447210031;
    nids.push_back(6631);
    tri3_xyze(0, 2) = -0.023352792493704529336;
    tri3_xyze(1, 2) = 0.056378628357611280941;
    tri3_xyze(2, 2) = 0.83180830895537727354;
    nids.push_back(-896);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.029774575140616293956;
    tri3_xyze(1, 0) = 0.051571076917342986712;
    tri3_xyze(2, 0) = 0.82938926261447210031;
    nids.push_back(6631);
    tri3_xyze(0, 1) = -0.015412454212450349758;
    tri3_xyze(1, 1) = 0.057520062190231820154;
    tri3_xyze(2, 1) = 0.82938926261447198929;
    nids.push_back(6881);
    tri3_xyze(0, 2) = -0.023352792493704529336;
    tri3_xyze(1, 2) = 0.056378628357611280941;
    tri3_xyze(2, 2) = 0.83180830895537727354;
    nids.push_back(-896);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.015412454212450349758;
    tri3_xyze(1, 0) = 0.057520062190231820154;
    tri3_xyze(2, 0) = 0.82938926261447198929;
    nids.push_back(6881);
    tri3_xyze(0, 1) = -0.016448356307296806839;
    tri3_xyze(1, 1) = 0.061386101439866169738;
    tri3_xyze(2, 1) = 0.83422735529628255779;
    nids.push_back(6879);
    tri3_xyze(0, 2) = -0.023352792493704529336;
    tri3_xyze(1, 2) = 0.056378628357611280941;
    tri3_xyze(2, 2) = 0.83180830895537727354;
    nids.push_back(-896);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.016448356307296806839;
    tri3_xyze(1, 0) = 0.061386101439866169738;
    tri3_xyze(2, 0) = 0.83422735529628255779;
    nids.push_back(6879);
    tri3_xyze(0, 1) = -0.031775784314454665058;
    tri3_xyze(1, 1) = 0.055037272883004147161;
    tri3_xyze(2, 1) = 0.83422735529628266882;
    nids.push_back(6629);
    tri3_xyze(0, 2) = -0.023352792493704529336;
    tri3_xyze(1, 2) = 0.056378628357611280941;
    tri3_xyze(2, 2) = 0.83180830895537727354;
    nids.push_back(-896);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.029774575140616293956;
    tri3_xyze(1, 0) = 0.051571076917342986712;
    tri3_xyze(2, 0) = 0.82938926261447210031;
    nids.push_back(6631);
    tri3_xyze(0, 1) = -0.028092332998893382312;
    tri3_xyze(1, 1) = 0.048657348057245253248;
    tri3_xyze(2, 1) = 0.82408768370493423738;
    nids.push_back(6633);
    tri3_xyze(0, 2) = -0.021955255988726610911;
    tri3_xyze(1, 2) = 0.053004676773381921573;
    tri3_xyze(2, 2) = 0.82673847315970305782;
    nids.push_back(-897);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.028092332998893382312;
    tri3_xyze(1, 0) = 0.048657348057245253248;
    tri3_xyze(2, 0) = 0.82408768370493423738;
    nids.push_back(6633);
    tri3_xyze(0, 1) = -0.014541661602946407209;
    tri3_xyze(1, 1) = 0.054270219928707633117;
    tri3_xyze(2, 1) = 0.82408768370493412636;
    nids.push_back(6883);
    tri3_xyze(0, 2) = -0.021955255988726610911;
    tri3_xyze(1, 2) = 0.053004676773381921573;
    tri3_xyze(2, 2) = 0.82673847315970305782;
    nids.push_back(-897);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.014541661602946407209;
    tri3_xyze(1, 0) = 0.054270219928707633117;
    tri3_xyze(2, 0) = 0.82408768370493412636;
    nids.push_back(6883);
    tri3_xyze(0, 1) = -0.015412454212450349758;
    tri3_xyze(1, 1) = 0.057520062190231820154;
    tri3_xyze(2, 1) = 0.82938926261447198929;
    nids.push_back(6881);
    tri3_xyze(0, 2) = -0.021955255988726610911;
    tri3_xyze(1, 2) = 0.053004676773381921573;
    tri3_xyze(2, 2) = 0.82673847315970305782;
    nids.push_back(-897);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.015412454212450349758;
    tri3_xyze(1, 0) = 0.057520062190231820154;
    tri3_xyze(2, 0) = 0.82938926261447198929;
    nids.push_back(6881);
    tri3_xyze(0, 1) = -0.029774575140616293956;
    tri3_xyze(1, 1) = 0.051571076917342986712;
    tri3_xyze(2, 1) = 0.82938926261447210031;
    nids.push_back(6631);
    tri3_xyze(0, 2) = -0.021955255988726610911;
    tri3_xyze(1, 2) = 0.053004676773381921573;
    tri3_xyze(2, 2) = 0.82673847315970305782;
    nids.push_back(-897);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.028092332998893382312;
    tri3_xyze(1, 0) = 0.048657348057245253248;
    tri3_xyze(2, 0) = 0.82408768370493423738;
    nids.push_back(6633);
    tri3_xyze(0, 1) = -0.02675558785278368823;
    tri3_xyze(1, 1) = 0.046342037547411993326;
    tri3_xyze(2, 1) = 0.81840622763408255036;
    nids.push_back(6635);
    tri3_xyze(0, 2) = -0.020809823463261691023;
    tri3_xyze(1, 2) = 0.050239358035621128895;
    tri3_xyze(2, 2) = 0.82124695566950822734;
    nids.push_back(-898);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.02675558785278368823;
    tri3_xyze(1, 0) = 0.046342037547411993326;
    tri3_xyze(2, 0) = 0.81840622763408255036;
    nids.push_back(6635);
    tri3_xyze(0, 1) = -0.013849711398423282874;
    tri3_xyze(1, 1) = 0.051687826609119663646;
    tri3_xyze(2, 1) = 0.81840622763408232831;
    nids.push_back(6885);
    tri3_xyze(0, 2) = -0.020809823463261691023;
    tri3_xyze(1, 2) = 0.050239358035621128895;
    tri3_xyze(2, 2) = 0.82124695566950822734;
    nids.push_back(-898);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.013849711398423282874;
    tri3_xyze(1, 0) = 0.051687826609119663646;
    tri3_xyze(2, 0) = 0.81840622763408232831;
    nids.push_back(6885);
    tri3_xyze(0, 1) = -0.014541661602946407209;
    tri3_xyze(1, 1) = 0.054270219928707633117;
    tri3_xyze(2, 1) = 0.82408768370493412636;
    nids.push_back(6883);
    tri3_xyze(0, 2) = -0.020809823463261691023;
    tri3_xyze(1, 2) = 0.050239358035621128895;
    tri3_xyze(2, 2) = 0.82124695566950822734;
    nids.push_back(-898);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.014541661602946407209;
    tri3_xyze(1, 0) = 0.054270219928707633117;
    tri3_xyze(2, 0) = 0.82408768370493412636;
    nids.push_back(6883);
    tri3_xyze(0, 1) = -0.028092332998893382312;
    tri3_xyze(1, 1) = 0.048657348057245253248;
    tri3_xyze(2, 1) = 0.82408768370493423738;
    nids.push_back(6633);
    tri3_xyze(0, 2) = -0.020809823463261691023;
    tri3_xyze(1, 2) = 0.050239358035621128895;
    tri3_xyze(2, 2) = 0.82124695566950822734;
    nids.push_back(-898);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.013347516066957392997;
    tri3_xyze(1, 0) = 0.04981360811676469047;
    tri3_xyze(2, 0) = 0.81243449435809123038;
    nids.push_back(6887);
    tri3_xyze(0, 1) = -0.013849711398423282874;
    tri3_xyze(1, 1) = 0.051687826609119663646;
    tri3_xyze(2, 1) = 0.81840622763408232831;
    nids.push_back(6885);
    tri3_xyze(0, 2) = -0.019934559072484637349;
    tri3_xyze(1, 2) = 0.048126282872744771679;
    tri3_xyze(2, 2) = 0.81542036099608683486;
    nids.push_back(-899);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.013849711398423282874;
    tri3_xyze(1, 0) = 0.051687826609119663646;
    tri3_xyze(2, 0) = 0.81840622763408232831;
    nids.push_back(6885);
    tri3_xyze(0, 1) = -0.02675558785278368823;
    tri3_xyze(1, 1) = 0.046342037547411993326;
    tri3_xyze(2, 1) = 0.81840622763408255036;
    nids.push_back(6635);
    tri3_xyze(0, 2) = -0.019934559072484637349;
    tri3_xyze(1, 2) = 0.048126282872744771679;
    tri3_xyze(2, 2) = 0.81542036099608683486;
    nids.push_back(-899);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.013042995528922771836;
    tri3_xyze(1, 0) = 0.048677121996870686971;
    tri3_xyze(2, 0) = 0.79373333832163328605;
    nids.push_back(6792);
    tri3_xyze(0, 1) = 1.0003085758762357486e-14;
    tri3_xyze(1, 1) = 0.050394264934275928569;
    tri3_xyze(2, 1) = 0.79373333832163317503;
    nids.push_back(7042);
    tri3_xyze(0, 2) = -0.0064959869460047021256;
    tri3_xyze(1, 2) = 0.049341919561400007788;
    tri3_xyze(2, 2) = 0.79686666916074078149;
    nids.push_back(-901);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003085758762357486e-14;
    tri3_xyze(1, 0) = 0.050394264934275928569;
    tri3_xyze(2, 0) = 0.79373333832163317503;
    nids.push_back(7042);
    tri3_xyze(0, 1) = 9.9908151490063951676e-15;
    tri3_xyze(1, 1) = 0.050000000000000002776;
    tri3_xyze(2, 1) = 0.79999999999984849897;
    nids.push_back(7041);
    tri3_xyze(0, 2) = -0.0064959869460047021256;
    tri3_xyze(1, 2) = 0.049341919561400007788;
    tri3_xyze(2, 2) = 0.79686666916074078149;
    nids.push_back(-901);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.9908151490063951676e-15;
    tri3_xyze(1, 0) = 0.050000000000000002776;
    tri3_xyze(2, 0) = 0.79999999999984849897;
    nids.push_back(7041);
    tri3_xyze(0, 1) = -0.012940952255116031089;
    tri3_xyze(1, 1) = 0.048296291314453412835;
    tri3_xyze(2, 1) = 0.79999999999984849897;
    nids.push_back(6791);
    tri3_xyze(0, 2) = -0.0064959869460047021256;
    tri3_xyze(1, 2) = 0.049341919561400007788;
    tri3_xyze(2, 2) = 0.79686666916074078149;
    nids.push_back(-901);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.012940952255116031089;
    tri3_xyze(1, 0) = 0.048296291314453412835;
    tri3_xyze(2, 0) = 0.79999999999984849897;
    nids.push_back(6791);
    tri3_xyze(0, 1) = 9.9908151490063951676e-15;
    tri3_xyze(1, 1) = 0.050000000000000002776;
    tri3_xyze(2, 1) = 0.79999999999984849897;
    nids.push_back(7041);
    tri3_xyze(0, 2) = -0.0064959869460047021256;
    tri3_xyze(1, 2) = 0.049341919561400202077;
    tri3_xyze(2, 2) = 0.80313333083895610542;
    nids.push_back(-902);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.9908151490063951676e-15;
    tri3_xyze(1, 0) = 0.050000000000000002776;
    tri3_xyze(2, 0) = 0.79999999999984849897;
    nids.push_back(7041);
    tri3_xyze(0, 1) = 1.0003085758762357486e-14;
    tri3_xyze(1, 1) = 0.050394264934276310208;
    tri3_xyze(2, 1) = 0.80626666167806371188;
    nids.push_back(7139);
    tri3_xyze(0, 2) = -0.0064959869460047021256;
    tri3_xyze(1, 2) = 0.049341919561400202077;
    tri3_xyze(2, 2) = 0.80313333083895610542;
    nids.push_back(-902);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003085758762357486e-14;
    tri3_xyze(1, 0) = 0.050394264934276310208;
    tri3_xyze(2, 0) = 0.80626666167806371188;
    nids.push_back(7139);
    tri3_xyze(0, 1) = -0.013042995528922771836;
    tri3_xyze(1, 1) = 0.048677121996871075549;
    tri3_xyze(2, 1) = 0.80626666167806371188;
    nids.push_back(6889);
    tri3_xyze(0, 2) = -0.0064959869460047021256;
    tri3_xyze(1, 2) = 0.049341919561400202077;
    tri3_xyze(2, 2) = 0.80313333083895610542;
    nids.push_back(-902);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.013347516066957389527;
    tri3_xyze(1, 0) = 0.049813608116763899436;
    tri3_xyze(2, 0) = 0.78756550564160587857;
    nids.push_back(6795);
    tri3_xyze(0, 1) = 1.0003157803325775479e-14;
    tri3_xyze(1, 1) = 0.051570841943568049293;
    tri3_xyze(2, 1) = 0.78756550564160576755;
    nids.push_back(7045);
    tri3_xyze(0, 2) = -0.0065976278989650395668;
    tri3_xyze(1, 2) = 0.050113959247869641067;
    tri3_xyze(2, 2) = 0.79064942198161958231;
    nids.push_back(-903);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003157803325775479e-14;
    tri3_xyze(1, 0) = 0.051570841943568049293;
    tri3_xyze(2, 0) = 0.78756550564160576755;
    nids.push_back(7045);
    tri3_xyze(0, 1) = 1.0003085758762357486e-14;
    tri3_xyze(1, 1) = 0.050394264934275928569;
    tri3_xyze(2, 1) = 0.79373333832163317503;
    nids.push_back(7042);
    tri3_xyze(0, 2) = -0.0065976278989650395668;
    tri3_xyze(1, 2) = 0.050113959247869641067;
    tri3_xyze(2, 2) = 0.79064942198161958231;
    nids.push_back(-903);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003085758762357486e-14;
    tri3_xyze(1, 0) = 0.050394264934275928569;
    tri3_xyze(2, 0) = 0.79373333832163317503;
    nids.push_back(7042);
    tri3_xyze(0, 1) = -0.013042995528922771836;
    tri3_xyze(1, 1) = 0.048677121996870686971;
    tri3_xyze(2, 1) = 0.79373333832163328605;
    nids.push_back(6792);
    tri3_xyze(0, 2) = -0.0065976278989650395668;
    tri3_xyze(1, 2) = 0.050113959247869641067;
    tri3_xyze(2, 2) = 0.79064942198161958231;
    nids.push_back(-903);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.013849711398423343589;
    tri3_xyze(1, 0) = 0.051687826609118546484;
    tri3_xyze(2, 0) = 0.78159377236561433655;
    nids.push_back(6797);
    tri3_xyze(0, 1) = 9.9557491430278385319e-15;
    tri3_xyze(1, 1) = 0.053511175705586895568;
    tri3_xyze(2, 1) = 0.78159377236561444757;
    nids.push_back(7047);
    tri3_xyze(0, 2) = -0.0067993068663401929133;
    tri3_xyze(1, 2) = 0.051645863093759347695;
    tri3_xyze(2, 2) = 0.78457963900361016307;
    nids.push_back(-904);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.9557491430278385319e-15;
    tri3_xyze(1, 0) = 0.053511175705586895568;
    tri3_xyze(2, 0) = 0.78159377236561444757;
    nids.push_back(7047);
    tri3_xyze(0, 1) = 1.0003157803325775479e-14;
    tri3_xyze(1, 1) = 0.051570841943568049293;
    tri3_xyze(2, 1) = 0.78756550564160576755;
    nids.push_back(7045);
    tri3_xyze(0, 2) = -0.0067993068663401929133;
    tri3_xyze(1, 2) = 0.051645863093759347695;
    tri3_xyze(2, 2) = 0.78457963900361016307;
    nids.push_back(-904);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003157803325775479e-14;
    tri3_xyze(1, 0) = 0.051570841943568049293;
    tri3_xyze(2, 0) = 0.78756550564160576755;
    nids.push_back(7045);
    tri3_xyze(0, 1) = -0.013347516066957389527;
    tri3_xyze(1, 1) = 0.049813608116763899436;
    tri3_xyze(2, 1) = 0.78756550564160587857;
    nids.push_back(6795);
    tri3_xyze(0, 2) = -0.0067993068663401929133;
    tri3_xyze(1, 2) = 0.051645863093759347695;
    tri3_xyze(2, 2) = 0.78457963900361016307;
    nids.push_back(-904);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.014541661602946424556;
    tri3_xyze(1, 0) = 0.054270219928706189827;
    tri3_xyze(2, 0) = 0.7759123162947625385;
    nids.push_back(6799);
    tri3_xyze(0, 1) = 1.0053342326426079629e-14;
    tri3_xyze(1, 1) = 0.056184665997806082505;
    tri3_xyze(2, 1) = 0.77591231629476242748;
    nids.push_back(7049);
    tri3_xyze(0, 2) = -0.0070978432503374395274;
    tri3_xyze(1, 2) = 0.053913472060304425126;
    tri3_xyze(2, 2) = 0.77875304433018843753;
    nids.push_back(-905);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0053342326426079629e-14;
    tri3_xyze(1, 0) = 0.056184665997806082505;
    tri3_xyze(2, 0) = 0.77591231629476242748;
    nids.push_back(7049);
    tri3_xyze(0, 1) = 9.9557491430278385319e-15;
    tri3_xyze(1, 1) = 0.053511175705586895568;
    tri3_xyze(2, 1) = 0.78159377236561444757;
    nids.push_back(7047);
    tri3_xyze(0, 2) = -0.0070978432503374395274;
    tri3_xyze(1, 2) = 0.053913472060304425126;
    tri3_xyze(2, 2) = 0.77875304433018843753;
    nids.push_back(-905);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.9557491430278385319e-15;
    tri3_xyze(1, 0) = 0.053511175705586895568;
    tri3_xyze(2, 0) = 0.78159377236561444757;
    nids.push_back(7047);
    tri3_xyze(0, 1) = -0.013849711398423343589;
    tri3_xyze(1, 1) = 0.051687826609118546484;
    tri3_xyze(2, 1) = 0.78159377236561433655;
    nids.push_back(6797);
    tri3_xyze(0, 2) = -0.0070978432503374395274;
    tri3_xyze(1, 2) = 0.053913472060304425126;
    tri3_xyze(2, 2) = 0.77875304433018843753;
    nids.push_back(-905);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.015412454212450197102;
    tri3_xyze(1, 0) = 0.057520062190230016042;
    tri3_xyze(2, 0) = 0.77061073738522456456;
    nids.push_back(6801);
    tri3_xyze(0, 1) = 1.0056536604005482796e-14;
    tri3_xyze(1, 1) = 0.059549150281251718442;
    tri3_xyze(2, 1) = 0.77061073738522445353;
    nids.push_back(7051);
    tri3_xyze(0, 2) = -0.0074885289538441273186;
    tri3_xyze(1, 2) = 0.056881024599498501704;
    tri3_xyze(2, 2) = 0.77326152683999349602;
    nids.push_back(-906);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0056536604005482796e-14;
    tri3_xyze(1, 0) = 0.059549150281251718442;
    tri3_xyze(2, 0) = 0.77061073738522445353;
    nids.push_back(7051);
    tri3_xyze(0, 1) = 1.0053342326426079629e-14;
    tri3_xyze(1, 1) = 0.056184665997806082505;
    tri3_xyze(2, 1) = 0.77591231629476242748;
    nids.push_back(7049);
    tri3_xyze(0, 2) = -0.0074885289538441273186;
    tri3_xyze(1, 2) = 0.056881024599498501704;
    tri3_xyze(2, 2) = 0.77326152683999349602;
    nids.push_back(-906);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0053342326426079629e-14;
    tri3_xyze(1, 0) = 0.056184665997806082505;
    tri3_xyze(2, 0) = 0.77591231629476242748;
    nids.push_back(7049);
    tri3_xyze(0, 1) = -0.014541661602946424556;
    tri3_xyze(1, 1) = 0.054270219928706189827;
    tri3_xyze(2, 1) = 0.7759123162947625385;
    nids.push_back(6799);
    tri3_xyze(0, 2) = -0.0074885289538441273186;
    tri3_xyze(1, 2) = 0.056881024599498501704;
    tri3_xyze(2, 2) = 0.77326152683999349602;
    nids.push_back(-906);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.0164483563072969248;
    tri3_xyze(1, 0) = 0.061386101439864018681;
    tri3_xyze(2, 0) = 0.76577264470341366298;
    nids.push_back(6803);
    tri3_xyze(0, 1) = 1.0060336543049417096e-14;
    tri3_xyze(1, 1) = 0.063551568628928342686;
    tri3_xyze(2, 1) = 0.76577264470341366298;
    nids.push_back(7053);
    tri3_xyze(0, 2) = -0.0079652026299317519459;
    tri3_xyze(1, 2) = 0.060501720635068527432;
    tri3_xyze(2, 2) = 0.76819169104431916928;
    nids.push_back(-907);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0060336543049417096e-14;
    tri3_xyze(1, 0) = 0.063551568628928342686;
    tri3_xyze(2, 0) = 0.76577264470341366298;
    nids.push_back(7053);
    tri3_xyze(0, 1) = 1.0056536604005482796e-14;
    tri3_xyze(1, 1) = 0.059549150281251718442;
    tri3_xyze(2, 1) = 0.77061073738522445353;
    nids.push_back(7051);
    tri3_xyze(0, 2) = -0.0079652026299317519459;
    tri3_xyze(1, 2) = 0.060501720635068527432;
    tri3_xyze(2, 2) = 0.76819169104431916928;
    nids.push_back(-907);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0056536604005482796e-14;
    tri3_xyze(1, 0) = 0.059549150281251718442;
    tri3_xyze(2, 0) = 0.77061073738522445353;
    nids.push_back(7051);
    tri3_xyze(0, 1) = -0.015412454212450197102;
    tri3_xyze(1, 1) = 0.057520062190230016042;
    tri3_xyze(2, 1) = 0.77061073738522456456;
    nids.push_back(6801);
    tri3_xyze(0, 2) = -0.0079652026299317519459;
    tri3_xyze(1, 2) = 0.060501720635068527432;
    tri3_xyze(2, 2) = 0.76819169104431916928;
    nids.push_back(-907);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.017633031092632350106;
    tri3_xyze(1, 0) = 0.065807367929181903587;
    tri3_xyze(2, 0) = 0.76147433786105833686;
    nids.push_back(6805);
    tri3_xyze(0, 1) = 1.0004171685873872409e-14;
    tri3_xyze(1, 1) = 0.068128800512564396952;
    tri3_xyze(2, 1) = 0.76147433786105844788;
    nids.push_back(7055);
    tri3_xyze(0, 2) = -0.0085203468499773019063;
    tri3_xyze(1, 2) = 0.064718459627634658538;
    tri3_xyze(2, 2) = 0.76362349128223594441;
    nids.push_back(-908);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004171685873872409e-14;
    tri3_xyze(1, 0) = 0.068128800512564396952;
    tri3_xyze(2, 0) = 0.76147433786105844788;
    nids.push_back(7055);
    tri3_xyze(0, 1) = 1.0060336543049417096e-14;
    tri3_xyze(1, 1) = 0.063551568628928342686;
    tri3_xyze(2, 1) = 0.76577264470341366298;
    nids.push_back(7053);
    tri3_xyze(0, 2) = -0.0085203468499773019063;
    tri3_xyze(1, 2) = 0.064718459627634658538;
    tri3_xyze(2, 2) = 0.76362349128223594441;
    nids.push_back(-908);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0060336543049417096e-14;
    tri3_xyze(1, 0) = 0.063551568628928342686;
    tri3_xyze(2, 0) = 0.76577264470341366298;
    nids.push_back(7053);
    tri3_xyze(0, 1) = -0.0164483563072969248;
    tri3_xyze(1, 1) = 0.061386101439864018681;
    tri3_xyze(2, 1) = 0.76577264470341366298;
    nids.push_back(6803);
    tri3_xyze(0, 2) = -0.0085203468499773019063;
    tri3_xyze(1, 2) = 0.064718459627634658538;
    tri3_xyze(2, 2) = 0.76362349128223594441;
    nids.push_back(-908);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.018947795539401686427;
    tri3_xyze(1, 0) = 0.070714135644510123235;
    tri3_xyze(2, 0) = 0.75778360372474695694;
    nids.push_back(6807);
    tri3_xyze(0, 1) = 1.0004482737572315978e-14;
    tri3_xyze(1, 1) = 0.073208660251048829348;
    tri3_xyze(2, 1) = 0.75778360372474695694;
    nids.push_back(7057);
    tri3_xyze(0, 2) = -0.0091452066580035053234;
    tri3_xyze(1, 2) = 0.06946474108432632022;
    tri3_xyze(2, 2) = 0.75962897079290270241;
    nids.push_back(-909);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004482737572315978e-14;
    tri3_xyze(1, 0) = 0.073208660251048829348;
    tri3_xyze(2, 0) = 0.75778360372474695694;
    nids.push_back(7057);
    tri3_xyze(0, 1) = 1.0004171685873872409e-14;
    tri3_xyze(1, 1) = 0.068128800512564396952;
    tri3_xyze(2, 1) = 0.76147433786105844788;
    nids.push_back(7055);
    tri3_xyze(0, 2) = -0.0091452066580035053234;
    tri3_xyze(1, 2) = 0.06946474108432632022;
    tri3_xyze(2, 2) = 0.75962897079290270241;
    nids.push_back(-909);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004171685873872409e-14;
    tri3_xyze(1, 0) = 0.068128800512564396952;
    tri3_xyze(2, 0) = 0.76147433786105844788;
    nids.push_back(7055);
    tri3_xyze(0, 1) = -0.017633031092632350106;
    tri3_xyze(1, 1) = 0.065807367929181903587;
    tri3_xyze(2, 1) = 0.76147433786105833686;
    nids.push_back(6805);
    tri3_xyze(0, 2) = -0.0091452066580035053234;
    tri3_xyze(1, 2) = 0.06946474108432632022;
    tri3_xyze(2, 2) = 0.75962897079290270241;
    nids.push_back(-909);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.020371915026877178423;
    tri3_xyze(1, 0) = 0.076029021927816992088;
    tri3_xyze(2, 0) = 0.75475864737654663195;
    nids.push_back(6809);
    tri3_xyze(0, 1) = 1.0004819660879340093e-14;
    tri3_xyze(1, 1) = 0.078711035421744998675;
    tri3_xyze(2, 1) = 0.75475864737654663195;
    nids.push_back(7059);
    tri3_xyze(0, 2) = -0.0098299276415647132699;
    tri3_xyze(1, 2) = 0.074665713311280232367;
    tri3_xyze(2, 2) = 0.75627112555064679444;
    nids.push_back(-910);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004819660879340093e-14;
    tri3_xyze(1, 0) = 0.078711035421744998675;
    tri3_xyze(2, 0) = 0.75475864737654663195;
    nids.push_back(7059);
    tri3_xyze(0, 1) = 1.0004482737572315978e-14;
    tri3_xyze(1, 1) = 0.073208660251048829348;
    tri3_xyze(2, 1) = 0.75778360372474695694;
    nids.push_back(7057);
    tri3_xyze(0, 2) = -0.0098299276415647132699;
    tri3_xyze(1, 2) = 0.074665713311280232367;
    tri3_xyze(2, 2) = 0.75627112555064679444;
    nids.push_back(-910);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004482737572315978e-14;
    tri3_xyze(1, 0) = 0.073208660251048829348;
    tri3_xyze(2, 0) = 0.75778360372474695694;
    nids.push_back(7057);
    tri3_xyze(0, 1) = -0.018947795539401686427;
    tri3_xyze(1, 1) = 0.070714135644510123235;
    tri3_xyze(2, 1) = 0.75778360372474695694;
    nids.push_back(6807);
    tri3_xyze(0, 2) = -0.0098299276415647132699;
    tri3_xyze(1, 2) = 0.074665713311280232367;
    tri3_xyze(2, 2) = 0.75627112555064679444;
    nids.push_back(-910);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.02188293034001330406;
    tri3_xyze(1, 0) = 0.081668207847456070203;
    tri3_xyze(2, 0) = 0.75244717418508977325;
    nids.push_back(6811);
    tri3_xyze(0, 1) = 1.0005177142313128246e-14;
    tri3_xyze(1, 1) = 0.084549150281251053696;
    tri3_xyze(2, 1) = 0.75244717418508977325;
    nids.push_back(7061);
    tri3_xyze(0, 2) = -0.010563711341717617678;
    tri3_xyze(1, 2) = 0.080239353869567278665;
    tri3_xyze(2, 2) = 0.7536029107808182026;
    nids.push_back(-911);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005177142313128246e-14;
    tri3_xyze(1, 0) = 0.084549150281251053696;
    tri3_xyze(2, 0) = 0.75244717418508977325;
    nids.push_back(7061);
    tri3_xyze(0, 1) = 1.0004819660879340093e-14;
    tri3_xyze(1, 1) = 0.078711035421744998675;
    tri3_xyze(2, 1) = 0.75475864737654663195;
    nids.push_back(7059);
    tri3_xyze(0, 2) = -0.010563711341717617678;
    tri3_xyze(1, 2) = 0.080239353869567278665;
    tri3_xyze(2, 2) = 0.7536029107808182026;
    nids.push_back(-911);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004819660879340093e-14;
    tri3_xyze(1, 0) = 0.078711035421744998675;
    tri3_xyze(2, 0) = 0.75475864737654663195;
    nids.push_back(7059);
    tri3_xyze(0, 1) = -0.020371915026877178423;
    tri3_xyze(1, 1) = 0.076029021927816992088;
    tri3_xyze(2, 1) = 0.75475864737654663195;
    nids.push_back(6809);
    tri3_xyze(0, 2) = -0.010563711341717617678;
    tri3_xyze(1, 2) = 0.080239353869567278665;
    tri3_xyze(2, 2) = 0.7536029107808182026;
    nids.push_back(-911);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.023457011864685423791;
    tri3_xyze(1, 0) = 0.087542760072787884074;
    tri3_xyze(2, 0) = 0.75088563746341285832;
    nids.push_back(6813);
    tri3_xyze(0, 1) = 1.0005549544177918912e-14;
    tri3_xyze(1, 1) = 0.090630934270712198986;
    tri3_xyze(2, 1) = 0.75088563746341274729;
    nids.push_back(7063);
    tri3_xyze(0, 2) = -0.011334985551169678153;
    tri3_xyze(1, 2) = 0.08609776311805179827;
    tri3_xyze(2, 2) = 0.75166640582425126027;
    nids.push_back(-912);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005549544177918912e-14;
    tri3_xyze(1, 0) = 0.090630934270712198986;
    tri3_xyze(2, 0) = 0.75088563746341274729;
    nids.push_back(7063);
    tri3_xyze(0, 1) = 1.0005177142313128246e-14;
    tri3_xyze(1, 1) = 0.084549150281251053696;
    tri3_xyze(2, 1) = 0.75244717418508977325;
    nids.push_back(7061);
    tri3_xyze(0, 2) = -0.011334985551169678153;
    tri3_xyze(1, 2) = 0.08609776311805179827;
    tri3_xyze(2, 2) = 0.75166640582425126027;
    nids.push_back(-912);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005177142313128246e-14;
    tri3_xyze(1, 0) = 0.084549150281251053696;
    tri3_xyze(2, 0) = 0.75244717418508977325;
    nids.push_back(7061);
    tri3_xyze(0, 1) = -0.02188293034001330406;
    tri3_xyze(1, 1) = 0.081668207847456070203;
    tri3_xyze(2, 1) = 0.75244717418508977325;
    nids.push_back(6811);
    tri3_xyze(0, 2) = -0.011334985551169678153;
    tri3_xyze(1, 2) = 0.08609776311805179827;
    tri3_xyze(2, 2) = 0.75166640582425126027;
    nids.push_back(-912);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.025069335394938625089;
    tri3_xyze(1, 0) = 0.09356003340593162243;
    tri3_xyze(2, 0) = 0.75009866357843346663;
    nids.push_back(6815);
    tri3_xyze(0, 1) = 1.0005930993473840787e-14;
    tri3_xyze(1, 1) = 0.096860474023532747156;
    tri3_xyze(2, 1) = 0.75009866357843346663;
    nids.push_back(7065);
    tri3_xyze(0, 2) = -0.012131586814901009277;
    tri3_xyze(1, 2) = 0.092148550443241106223;
    tri3_xyze(2, 2) = 0.75049215052092321798;
    nids.push_back(-913);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005930993473840787e-14;
    tri3_xyze(1, 0) = 0.096860474023532747156;
    tri3_xyze(2, 0) = 0.75009866357843346663;
    nids.push_back(7065);
    tri3_xyze(0, 1) = 1.0005549544177918912e-14;
    tri3_xyze(1, 1) = 0.090630934270712198986;
    tri3_xyze(2, 1) = 0.75088563746341274729;
    nids.push_back(7063);
    tri3_xyze(0, 2) = -0.012131586814901009277;
    tri3_xyze(1, 2) = 0.092148550443241106223;
    tri3_xyze(2, 2) = 0.75049215052092321798;
    nids.push_back(-913);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005549544177918912e-14;
    tri3_xyze(1, 0) = 0.090630934270712198986;
    tri3_xyze(2, 0) = 0.75088563746341274729;
    nids.push_back(7063);
    tri3_xyze(0, 1) = -0.023457011864685423791;
    tri3_xyze(1, 1) = 0.087542760072787884074;
    tri3_xyze(2, 1) = 0.75088563746341285832;
    nids.push_back(6813);
    tri3_xyze(0, 2) = -0.012131586814901009277;
    tri3_xyze(1, 2) = 0.092148550443241106223;
    tri3_xyze(2, 2) = 0.75049215052092321798;
    nids.push_back(-913);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0006315474517633335e-14;
    tri3_xyze(1, 0) = 0.10313952597646400267;
    tri3_xyze(2, 0) = 0.75009866357843324458;
    nids.push_back(7067);
    tri3_xyze(0, 1) = 1.0005930993473840787e-14;
    tri3_xyze(1, 1) = 0.096860474023532747156;
    tri3_xyze(2, 1) = 0.75009866357843346663;
    nids.push_back(7065);
    tri3_xyze(0, 2) = -0.012940952255116039762;
    tri3_xyze(1, 2) = 0.098296291314451819665;
    tri3_xyze(2, 2) = 0.75009866357843346663;
    nids.push_back(-914);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005930993473840787e-14;
    tri3_xyze(1, 0) = 0.096860474023532747156;
    tri3_xyze(2, 0) = 0.75009866357843346663;
    nids.push_back(7065);
    tri3_xyze(0, 1) = -0.025069335394938625089;
    tri3_xyze(1, 1) = 0.09356003340593162243;
    tri3_xyze(2, 1) = 0.75009866357843346663;
    nids.push_back(6815);
    tri3_xyze(0, 2) = -0.012940952255116039762;
    tri3_xyze(1, 2) = 0.098296291314451819665;
    tri3_xyze(2, 2) = 0.75009866357843346663;
    nids.push_back(-914);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.038720813491561319186;
    tri3_xyze(1, 0) = 0.14450804326094301988;
    tri3_xyze(2, 0) = 0.80626666167806060326;
    nids.push_back(6843);
    tri3_xyze(0, 1) = 9.7434075385316158324e-15;
    tri3_xyze(1, 1) = 0.14960573506572408253;
    tri3_xyze(2, 1) = 0.80626666167806060326;
    nids.push_back(7093);
    tri3_xyze(0, 2) = -0.019385917564227441584;
    tri3_xyze(1, 2) = 0.14725066306750683731;
    tri3_xyze(2, 2) = 0.8031333308389529968;
    nids.push_back(-927);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.7434075385316158324e-15;
    tri3_xyze(1, 0) = 0.14960573506572408253;
    tri3_xyze(2, 0) = 0.80626666167806060326;
    nids.push_back(7093);
    tri3_xyze(0, 1) = 9.8759580880385857084e-15;
    tri3_xyze(1, 1) = 0.1500000000000000222;
    tri3_xyze(2, 1) = 0.79999999999984527932;
    nids.push_back(7091);
    tri3_xyze(0, 2) = -0.019385917564227441584;
    tri3_xyze(1, 2) = 0.14725066306750683731;
    tri3_xyze(2, 2) = 0.8031333308389529968;
    nids.push_back(-927);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.038416292953526692822;
    tri3_xyze(1, 0) = 0.14337155714104979354;
    tri3_xyze(2, 0) = 0.8124344943580883438;
    nids.push_back(6845);
    tri3_xyze(0, 1) = 9.7454255145857347186e-15;
    tri3_xyze(1, 1) = 0.14842915805643194793;
    tri3_xyze(2, 1) = 0.81243449435808812176;
    nids.push_back(7095);
    tri3_xyze(0, 2) = -0.019284276611267130164;
    tri3_xyze(1, 2) = 0.14647862338103720403;
    tri3_xyze(2, 2) = 0.80935057801807452904;
    nids.push_back(-928);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.7454255145857347186e-15;
    tri3_xyze(1, 0) = 0.14842915805643194793;
    tri3_xyze(2, 0) = 0.81243449435808812176;
    nids.push_back(7095);
    tri3_xyze(0, 1) = 9.7434075385316158324e-15;
    tri3_xyze(1, 1) = 0.14960573506572408253;
    tri3_xyze(2, 1) = 0.80626666167806060326;
    nids.push_back(7093);
    tri3_xyze(0, 2) = -0.019284276611267130164;
    tri3_xyze(1, 2) = 0.14647862338103720403;
    tri3_xyze(2, 2) = 0.80935057801807452904;
    nids.push_back(-928);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.7434075385316158324e-15;
    tri3_xyze(1, 0) = 0.14960573506572408253;
    tri3_xyze(2, 0) = 0.80626666167806060326;
    nids.push_back(7093);
    tri3_xyze(0, 1) = -0.038720813491561319186;
    tri3_xyze(1, 1) = 0.14450804326094301988;
    tri3_xyze(2, 1) = 0.80626666167806060326;
    nids.push_back(6843);
    tri3_xyze(0, 2) = -0.019284276611267130164;
    tri3_xyze(1, 2) = 0.14647862338103720403;
    tri3_xyze(2, 2) = 0.80935057801807452904;
    nids.push_back(-928);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.037914097622060743964;
    tri3_xyze(1, 0) = 0.14149733864869518118;
    tri3_xyze(2, 0) = 0.81840622763407966378;
    nids.push_back(6847);
    tri3_xyze(0, 1) = 1.0008969853489149771e-14;
    tri3_xyze(1, 1) = 0.14648882429441315023;
    tri3_xyze(2, 1) = 0.81840622763407933071;
    nids.push_back(7097);
    tri3_xyze(0, 2) = -0.019082597643891922173;
    tri3_xyze(1, 2) = 0.14494671953514751128;
    tri3_xyze(2, 2) = 0.81542036099608394828;
    nids.push_back(-929);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008969853489149771e-14;
    tri3_xyze(1, 0) = 0.14648882429441315023;
    tri3_xyze(2, 0) = 0.81840622763407933071;
    nids.push_back(7097);
    tri3_xyze(0, 1) = 9.7454255145857347186e-15;
    tri3_xyze(1, 1) = 0.14842915805643194793;
    tri3_xyze(2, 1) = 0.81243449435808812176;
    nids.push_back(7095);
    tri3_xyze(0, 2) = -0.019082597643891922173;
    tri3_xyze(1, 2) = 0.14494671953514751128;
    tri3_xyze(2, 2) = 0.81542036099608394828;
    nids.push_back(-929);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.7454255145857347186e-15;
    tri3_xyze(1, 0) = 0.14842915805643194793;
    tri3_xyze(2, 0) = 0.81243449435808812176;
    nids.push_back(7095);
    tri3_xyze(0, 1) = -0.038416292953526692822;
    tri3_xyze(1, 1) = 0.14337155714104979354;
    tri3_xyze(2, 1) = 0.8124344943580883438;
    nids.push_back(6845);
    tri3_xyze(0, 2) = -0.019082597643891922173;
    tri3_xyze(1, 2) = 0.14494671953514751128;
    tri3_xyze(2, 2) = 0.81542036099608394828;
    nids.push_back(-929);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.037222147417537689018;
    tri3_xyze(1, 0) = 0.13891494532910761417;
    tri3_xyze(2, 0) = 0.8240876837049313508;
    nids.push_back(6849);
    tri3_xyze(0, 1) = 1.0008806149422705407e-14;
    tri3_xyze(1, 1) = 0.14381533400219395635;
    tri3_xyze(2, 1) = 0.82408768370493123978;
    nids.push_back(7099);
    tri3_xyze(0, 2) = -0.018784061259894603568;
    tri3_xyze(1, 2) = 0.14267911056860246855;
    tri3_xyze(2, 2) = 0.82124695566950545178;
    nids.push_back(-930);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008806149422705407e-14;
    tri3_xyze(1, 0) = 0.14381533400219395635;
    tri3_xyze(2, 0) = 0.82408768370493123978;
    nids.push_back(7099);
    tri3_xyze(0, 1) = 1.0008969853489149771e-14;
    tri3_xyze(1, 1) = 0.14648882429441315023;
    tri3_xyze(2, 1) = 0.81840622763407933071;
    nids.push_back(7097);
    tri3_xyze(0, 2) = -0.018784061259894603568;
    tri3_xyze(1, 2) = 0.14267911056860246855;
    tri3_xyze(2, 2) = 0.82124695566950545178;
    nids.push_back(-930);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008969853489149771e-14;
    tri3_xyze(1, 0) = 0.14648882429441315023;
    tri3_xyze(2, 0) = 0.81840622763407933071;
    nids.push_back(7097);
    tri3_xyze(0, 1) = -0.037914097622060743964;
    tri3_xyze(1, 1) = 0.14149733864869518118;
    tri3_xyze(2, 1) = 0.81840622763407966378;
    nids.push_back(6847);
    tri3_xyze(0, 2) = -0.018784061259894603568;
    tri3_xyze(1, 2) = 0.14267911056860246855;
    tri3_xyze(2, 2) = 0.82124695566950545178;
    nids.push_back(-930);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.036351354808033670141;
    tri3_xyze(1, 0) = 0.13566510306758378102;
    tri3_xyze(2, 0) = 0.8293892626144695468;
    nids.push_back(6851);
    tri3_xyze(0, 1) = 1.000860013417727975e-14;
    tri3_xyze(1, 1) = 0.14045084971874827184;
    tri3_xyze(2, 1) = 0.82938926261446932475;
    nids.push_back(7101);
    tri3_xyze(0, 2) = -0.018393375556387835112;
    tri3_xyze(1, 2) = 0.13971155802940840585;
    tri3_xyze(2, 2) = 0.82673847315970039329;
    nids.push_back(-931);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.000860013417727975e-14;
    tri3_xyze(1, 0) = 0.14045084971874827184;
    tri3_xyze(2, 0) = 0.82938926261446932475;
    nids.push_back(7101);
    tri3_xyze(0, 1) = 1.0008806149422705407e-14;
    tri3_xyze(1, 1) = 0.14381533400219395635;
    tri3_xyze(2, 1) = 0.82408768370493123978;
    nids.push_back(7099);
    tri3_xyze(0, 2) = -0.018393375556387835112;
    tri3_xyze(1, 2) = 0.13971155802940840585;
    tri3_xyze(2, 2) = 0.82673847315970039329;
    nids.push_back(-931);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008806149422705407e-14;
    tri3_xyze(1, 0) = 0.14381533400219395635;
    tri3_xyze(2, 0) = 0.82408768370493123978;
    nids.push_back(7099);
    tri3_xyze(0, 1) = -0.037222147417537689018;
    tri3_xyze(1, 1) = 0.13891494532910761417;
    tri3_xyze(2, 1) = 0.8240876837049313508;
    nids.push_back(6849);
    tri3_xyze(0, 2) = -0.018393375556387835112;
    tri3_xyze(1, 2) = 0.13971155802940840585;
    tri3_xyze(2, 2) = 0.82673847315970039329;
    nids.push_back(-931);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.035315452713187284184;
    tri3_xyze(1, 0) = 0.13179906381794970205;
    tri3_xyze(2, 0) = 0.83422735529628033735;
    nids.push_back(6853);
    tri3_xyze(0, 1) = 1.0008355056736363408e-14;
    tri3_xyze(1, 1) = 0.13644843137107168229;
    tri3_xyze(2, 1) = 0.8342273552962801153;
    nids.push_back(7103);
    tri3_xyze(0, 2) = -0.017916701880300235639;
    tri3_xyze(1, 2) = 0.13609086199383835236;
    tri3_xyze(2, 2) = 0.83180830895537483105;
    nids.push_back(-932);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008355056736363408e-14;
    tri3_xyze(1, 0) = 0.13644843137107168229;
    tri3_xyze(2, 0) = 0.8342273552962801153;
    nids.push_back(7103);
    tri3_xyze(0, 1) = 1.000860013417727975e-14;
    tri3_xyze(1, 1) = 0.14045084971874827184;
    tri3_xyze(2, 1) = 0.82938926261446932475;
    nids.push_back(7101);
    tri3_xyze(0, 2) = -0.017916701880300235639;
    tri3_xyze(1, 2) = 0.13609086199383835236;
    tri3_xyze(2, 2) = 0.83180830895537483105;
    nids.push_back(-932);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.000860013417727975e-14;
    tri3_xyze(1, 0) = 0.14045084971874827184;
    tri3_xyze(2, 0) = 0.82938926261446932475;
    nids.push_back(7101);
    tri3_xyze(0, 1) = -0.036351354808033670141;
    tri3_xyze(1, 1) = 0.13566510306758378102;
    tri3_xyze(2, 1) = 0.8293892626144695468;
    nids.push_back(6851);
    tri3_xyze(0, 2) = -0.017916701880300235639;
    tri3_xyze(1, 2) = 0.13609086199383835236;
    tri3_xyze(2, 2) = 0.83180830895537483105;
    nids.push_back(-932);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.034130777927851768672;
    tri3_xyze(1, 0) = 0.12737779732863194204;
    tri3_xyze(2, 0) = 0.83852566213863555245;
    nids.push_back(6855);
    tri3_xyze(0, 1) = 1.0008074782117600136e-14;
    tri3_xyze(1, 1) = 0.13187119948743578068;
    tri3_xyze(2, 1) = 0.8385256621386353304;
    nids.push_back(7105);
    tri3_xyze(0, 2) = -0.017361557660254760271;
    tri3_xyze(1, 2) = 0.13187412300127226983;
    tri3_xyze(2, 2) = 0.83637650871745783387;
    nids.push_back(-933);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008074782117600136e-14;
    tri3_xyze(1, 0) = 0.13187119948743578068;
    tri3_xyze(2, 0) = 0.8385256621386353304;
    nids.push_back(7105);
    tri3_xyze(0, 1) = 1.0008355056736363408e-14;
    tri3_xyze(1, 1) = 0.13644843137107168229;
    tri3_xyze(2, 1) = 0.8342273552962801153;
    nids.push_back(7103);
    tri3_xyze(0, 2) = -0.017361557660254760271;
    tri3_xyze(1, 2) = 0.13187412300127226983;
    tri3_xyze(2, 2) = 0.83637650871745783387;
    nids.push_back(-933);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008355056736363408e-14;
    tri3_xyze(1, 0) = 0.13644843137107168229;
    tri3_xyze(2, 0) = 0.8342273552962801153;
    nids.push_back(7103);
    tri3_xyze(0, 1) = -0.035315452713187284184;
    tri3_xyze(1, 1) = 0.13179906381794970205;
    tri3_xyze(2, 1) = 0.83422735529628033735;
    nids.push_back(6853);
    tri3_xyze(0, 2) = -0.017361557660254760271;
    tri3_xyze(1, 2) = 0.13187412300127226983;
    tri3_xyze(2, 2) = 0.83637650871745783387;
    nids.push_back(-933);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.03281601348108244276;
    tri3_xyze(1, 0) = 0.12247102961330372239;
    tri3_xyze(2, 0) = 0.84221639627494693237;
    nids.push_back(6857);
    tri3_xyze(0, 1) = 1.0007763730419158145e-14;
    tri3_xyze(1, 1) = 0.12679133974895126502;
    tri3_xyze(2, 1) = 0.84221639627494671032;
    nids.push_back(7107);
    tri3_xyze(0, 2) = -0.016736697852228549915;
    tri3_xyze(1, 2) = 0.12712784154458067754;
    tri3_xyze(2, 2) = 0.84037102920679107587;
    nids.push_back(-934);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0007763730419158145e-14;
    tri3_xyze(1, 0) = 0.12679133974895126502;
    tri3_xyze(2, 0) = 0.84221639627494671032;
    nids.push_back(7107);
    tri3_xyze(0, 1) = 1.0008074782117600136e-14;
    tri3_xyze(1, 1) = 0.13187119948743578068;
    tri3_xyze(2, 1) = 0.8385256621386353304;
    nids.push_back(7105);
    tri3_xyze(0, 2) = -0.016736697852228549915;
    tri3_xyze(1, 2) = 0.12712784154458067754;
    tri3_xyze(2, 2) = 0.84037102920679107587;
    nids.push_back(-934);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008074782117600136e-14;
    tri3_xyze(1, 0) = 0.13187119948743578068;
    tri3_xyze(2, 0) = 0.8385256621386353304;
    nids.push_back(7105);
    tri3_xyze(0, 1) = -0.034130777927851768672;
    tri3_xyze(1, 1) = 0.12737779732863194204;
    tri3_xyze(2, 1) = 0.83852566213863555245;
    nids.push_back(6855);
    tri3_xyze(0, 2) = -0.016736697852228549915;
    tri3_xyze(1, 2) = 0.12712784154458067754;
    tri3_xyze(2, 2) = 0.84037102920679107587;
    nids.push_back(-934);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.031391893993607040969;
    tri3_xyze(1, 0) = 0.11715614332999667313;
    tri3_xyze(2, 0) = 0.84524135262314736838;
    nids.push_back(6859);
    tri3_xyze(0, 1) = 1.0222879288884487869e-14;
    tri3_xyze(1, 1) = 0.12128896457825505406;
    tri3_xyze(2, 1) = 0.84524135262314725736;
    nids.push_back(7109);
    tri3_xyze(0, 2) = -0.016051976868667314213;
    tri3_xyze(1, 2) = 0.12192686931762666824;
    tri3_xyze(2, 2) = 0.84372887444904709486;
    nids.push_back(-935);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0222879288884487869e-14;
    tri3_xyze(1, 0) = 0.12128896457825505406;
    tri3_xyze(2, 0) = 0.84524135262314725736;
    nids.push_back(7109);
    tri3_xyze(0, 1) = 1.0007763730419158145e-14;
    tri3_xyze(1, 1) = 0.12679133974895126502;
    tri3_xyze(2, 1) = 0.84221639627494671032;
    nids.push_back(7107);
    tri3_xyze(0, 2) = -0.016051976868667314213;
    tri3_xyze(1, 2) = 0.12192686931762666824;
    tri3_xyze(2, 2) = 0.84372887444904709486;
    nids.push_back(-935);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0007763730419158145e-14;
    tri3_xyze(1, 0) = 0.12679133974895126502;
    tri3_xyze(2, 0) = 0.84221639627494671032;
    nids.push_back(7107);
    tri3_xyze(0, 1) = -0.03281601348108244276;
    tri3_xyze(1, 1) = 0.12247102961330372239;
    tri3_xyze(2, 1) = 0.84221639627494693237;
    nids.push_back(6857);
    tri3_xyze(0, 2) = -0.016051976868667314213;
    tri3_xyze(1, 2) = 0.12192686931762666824;
    tri3_xyze(2, 2) = 0.84372887444904709486;
    nids.push_back(-935);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.029880878680470679409;
    tri3_xyze(1, 0) = 0.11151695741035765053;
    tri3_xyze(2, 0) = 0.84755282581460422708;
    nids.push_back(6861);
    tri3_xyze(0, 1) = 1.0109610278934578421e-14;
    tri3_xyze(1, 1) = 0.11545084971874891577;
    tri3_xyze(2, 1) = 0.84755282581460400504;
    nids.push_back(7111);
    tri3_xyze(0, 2) = -0.015318193168514347355;
    tri3_xyze(1, 2) = 0.11635322875933958031;
    tri3_xyze(2, 2) = 0.84639708921887568671;
    nids.push_back(-936);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0109610278934578421e-14;
    tri3_xyze(1, 0) = 0.11545084971874891577;
    tri3_xyze(2, 0) = 0.84755282581460400504;
    nids.push_back(7111);
    tri3_xyze(0, 1) = 1.0222879288884487869e-14;
    tri3_xyze(1, 1) = 0.12128896457825505406;
    tri3_xyze(2, 1) = 0.84524135262314725736;
    nids.push_back(7109);
    tri3_xyze(0, 2) = -0.015318193168514347355;
    tri3_xyze(1, 2) = 0.11635322875933958031;
    tri3_xyze(2, 2) = 0.84639708921887568671;
    nids.push_back(-936);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0222879288884487869e-14;
    tri3_xyze(1, 0) = 0.12128896457825505406;
    tri3_xyze(2, 0) = 0.84524135262314725736;
    nids.push_back(7109);
    tri3_xyze(0, 1) = -0.031391893993607040969;
    tri3_xyze(1, 1) = 0.11715614332999667313;
    tri3_xyze(2, 1) = 0.84524135262314736838;
    nids.push_back(6859);
    tri3_xyze(0, 2) = -0.015318193168514347355;
    tri3_xyze(1, 2) = 0.11635322875933958031;
    tri3_xyze(2, 2) = 0.84639708921887568671;
    nids.push_back(-936);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.028306797155798663762;
    tri3_xyze(1, 0) = 0.10564240518502582278;
    tri3_xyze(2, 0) = 0.84911436253628114201;
    nids.push_back(6863);
    tri3_xyze(0, 1) = 9.8124184358865223691e-15;
    tri3_xyze(1, 1) = 0.10936906572928781212;
    tri3_xyze(2, 1) = 0.84911436253628103099;
    nids.push_back(7113);
    tri3_xyze(0, 2) = -0.014546918959062354534;
    tri3_xyze(1, 2) = 0.11049481951085506071;
    tri3_xyze(2, 2) = 0.84833359417544262904;
    nids.push_back(-937);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.8124184358865223691e-15;
    tri3_xyze(1, 0) = 0.10936906572928781212;
    tri3_xyze(2, 0) = 0.84911436253628103099;
    nids.push_back(7113);
    tri3_xyze(0, 1) = 1.0109610278934578421e-14;
    tri3_xyze(1, 1) = 0.11545084971874891577;
    tri3_xyze(2, 1) = 0.84755282581460400504;
    nids.push_back(7111);
    tri3_xyze(0, 2) = -0.014546918959062354534;
    tri3_xyze(1, 2) = 0.11049481951085506071;
    tri3_xyze(2, 2) = 0.84833359417544262904;
    nids.push_back(-937);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0109610278934578421e-14;
    tri3_xyze(1, 0) = 0.11545084971874891577;
    tri3_xyze(2, 0) = 0.84755282581460400504;
    nids.push_back(7111);
    tri3_xyze(0, 1) = -0.029880878680470679409;
    tri3_xyze(1, 1) = 0.11151695741035765053;
    tri3_xyze(2, 1) = 0.84755282581460422708;
    nids.push_back(6861);
    tri3_xyze(0, 2) = -0.014546918959062354534;
    tri3_xyze(1, 2) = 0.11049481951085506071;
    tri3_xyze(2, 2) = 0.84833359417544262904;
    nids.push_back(-937);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.028306797155798663762;
    tri3_xyze(1, 0) = 0.10564240518502582278;
    tri3_xyze(2, 0) = 0.84911436253628114201;
    nids.push_back(6863);
    tri3_xyze(0, 1) = -0.026694473625545441647;
    tri3_xyze(1, 1) = 0.099625131851881987277;
    tri3_xyze(2, 1) = 0.8499013364212605337;
    nids.push_back(6865);
    tri3_xyze(0, 2) = -0.013750317695331094534;
    tri3_xyze(1, 2) = 0.10444403218566569724;
    tri3_xyze(2, 2) = 0.84950784947877089337;
    nids.push_back(-938);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.026694473625545441647;
    tri3_xyze(1, 0) = 0.099625131851881987277;
    tri3_xyze(2, 0) = 0.8499013364212605337;
    nids.push_back(6865);
    tri3_xyze(0, 1) = 9.9147091733272353276e-15;
    tri3_xyze(1, 1) = 0.10313952597646713905;
    tri3_xyze(2, 1) = 0.84990133642126042268;
    nids.push_back(7115);
    tri3_xyze(0, 2) = -0.013750317695331094534;
    tri3_xyze(1, 2) = 0.10444403218566569724;
    tri3_xyze(2, 2) = 0.84950784947877089337;
    nids.push_back(-938);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.9147091733272353276e-15;
    tri3_xyze(1, 0) = 0.10313952597646713905;
    tri3_xyze(2, 0) = 0.84990133642126042268;
    nids.push_back(7115);
    tri3_xyze(0, 1) = 9.8124184358865223691e-15;
    tri3_xyze(1, 1) = 0.10936906572928781212;
    tri3_xyze(2, 1) = 0.84911436253628103099;
    nids.push_back(7113);
    tri3_xyze(0, 2) = -0.013750317695331094534;
    tri3_xyze(1, 2) = 0.10444403218566569724;
    tri3_xyze(2, 2) = 0.84950784947877089337;
    nids.push_back(-938);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.8124184358865223691e-15;
    tri3_xyze(1, 0) = 0.10936906572928781212;
    tri3_xyze(2, 0) = 0.84911436253628103099;
    nids.push_back(7113);
    tri3_xyze(0, 1) = -0.028306797155798663762;
    tri3_xyze(1, 1) = 0.10564240518502582278;
    tri3_xyze(2, 1) = 0.84911436253628114201;
    nids.push_back(6863);
    tri3_xyze(0, 2) = -0.013750317695331094534;
    tri3_xyze(1, 2) = 0.10444403218566569724;
    tri3_xyze(2, 2) = 0.84950784947877089337;
    nids.push_back(-938);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.026694473625545441647;
    tri3_xyze(1, 0) = 0.099625131851881987277;
    tri3_xyze(2, 0) = 0.8499013364212605337;
    nids.push_back(6865);
    tri3_xyze(0, 1) = -0.025069335394938645906;
    tri3_xyze(1, 1) = 0.093560033405934842077;
    tri3_xyze(2, 1) = 0.84990133642126064473;
    nids.push_back(6867);
    tri3_xyze(0, 2) = -0.012940952255116041497;
    tri3_xyze(1, 2) = 0.098296291314454983801;
    tri3_xyze(2, 2) = 0.84990133642126064473;
    nids.push_back(-939);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.025069335394938645906;
    tri3_xyze(1, 0) = 0.093560033405934842077;
    tri3_xyze(2, 0) = 0.84990133642126064473;
    nids.push_back(6867);
    tri3_xyze(0, 1) = 1.0005930993473840787e-14;
    tri3_xyze(1, 1) = 0.096860474023535966803;
    tri3_xyze(2, 1) = 0.84990133642126064473;
    nids.push_back(7117);
    tri3_xyze(0, 2) = -0.012940952255116041497;
    tri3_xyze(1, 2) = 0.098296291314454983801;
    tri3_xyze(2, 2) = 0.84990133642126064473;
    nids.push_back(-939);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005930993473840787e-14;
    tri3_xyze(1, 0) = 0.096860474023535966803;
    tri3_xyze(2, 0) = 0.84990133642126064473;
    nids.push_back(7117);
    tri3_xyze(0, 1) = 9.9147091733272353276e-15;
    tri3_xyze(1, 1) = 0.10313952597646713905;
    tri3_xyze(2, 1) = 0.84990133642126042268;
    nids.push_back(7115);
    tri3_xyze(0, 2) = -0.012940952255116041497;
    tri3_xyze(1, 2) = 0.098296291314454983801;
    tri3_xyze(2, 2) = 0.84990133642126064473;
    nids.push_back(-939);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.9147091733272353276e-15;
    tri3_xyze(1, 0) = 0.10313952597646713905;
    tri3_xyze(2, 0) = 0.84990133642126042268;
    nids.push_back(7115);
    tri3_xyze(0, 1) = -0.026694473625545441647;
    tri3_xyze(1, 1) = 0.099625131851881987277;
    tri3_xyze(2, 1) = 0.8499013364212605337;
    nids.push_back(6865);
    tri3_xyze(0, 2) = -0.012940952255116041497;
    tri3_xyze(1, 2) = 0.098296291314454983801;
    tri3_xyze(2, 2) = 0.84990133642126064473;
    nids.push_back(-939);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.025069335394938645906;
    tri3_xyze(1, 0) = 0.093560033405934842077;
    tri3_xyze(2, 0) = 0.84990133642126064473;
    nids.push_back(6867);
    tri3_xyze(0, 1) = -0.02345701186468549318;
    tri3_xyze(1, 1) = 0.087542760072790937187;
    tri3_xyze(2, 1) = 0.84911436253628169712;
    nids.push_back(6869);
    tri3_xyze(0, 2) = -0.012131586814901031829;
    tri3_xyze(1, 2) = 0.09214855044324425648;
    tri3_xyze(2, 2) = 0.84950784947877111541;
    nids.push_back(-940);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.02345701186468549318;
    tri3_xyze(1, 0) = 0.087542760072790937187;
    tri3_xyze(2, 0) = 0.84911436253628169712;
    nids.push_back(6869);
    tri3_xyze(0, 1) = 1.0005549544177918912e-14;
    tri3_xyze(1, 1) = 0.090630934270715279855;
    tri3_xyze(2, 1) = 0.8491143625362815861;
    nids.push_back(7119);
    tri3_xyze(0, 2) = -0.012131586814901031829;
    tri3_xyze(1, 2) = 0.09214855044324425648;
    tri3_xyze(2, 2) = 0.84950784947877111541;
    nids.push_back(-940);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005549544177918912e-14;
    tri3_xyze(1, 0) = 0.090630934270715279855;
    tri3_xyze(2, 0) = 0.8491143625362815861;
    nids.push_back(7119);
    tri3_xyze(0, 1) = 1.0005930993473840787e-14;
    tri3_xyze(1, 1) = 0.096860474023535966803;
    tri3_xyze(2, 1) = 0.84990133642126064473;
    nids.push_back(7117);
    tri3_xyze(0, 2) = -0.012131586814901031829;
    tri3_xyze(1, 2) = 0.09214855044324425648;
    tri3_xyze(2, 2) = 0.84950784947877111541;
    nids.push_back(-940);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005930993473840787e-14;
    tri3_xyze(1, 0) = 0.096860474023535966803;
    tri3_xyze(2, 0) = 0.84990133642126064473;
    nids.push_back(7117);
    tri3_xyze(0, 1) = -0.025069335394938645906;
    tri3_xyze(1, 1) = 0.093560033405934842077;
    tri3_xyze(2, 1) = 0.84990133642126064473;
    nids.push_back(6867);
    tri3_xyze(0, 2) = -0.012131586814901031829;
    tri3_xyze(1, 2) = 0.09214855044324425648;
    tri3_xyze(2, 2) = 0.84950784947877111541;
    nids.push_back(-940);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.02345701186468549318;
    tri3_xyze(1, 0) = 0.087542760072790937187;
    tri3_xyze(2, 0) = 0.84911436253628169712;
    nids.push_back(6869);
    tri3_xyze(0, 1) = -0.02188293034001336998;
    tri3_xyze(1, 1) = 0.081668207847459053927;
    tri3_xyze(2, 1) = 0.84755282581460522628;
    nids.push_back(6871);
    tri3_xyze(0, 2) = -0.011334985551169712847;
    tri3_xyze(1, 2) = 0.086097763118054851383;
    tri3_xyze(2, 2) = 0.84833359417544340619;
    nids.push_back(-941);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.02188293034001336998;
    tri3_xyze(1, 0) = 0.081668207847459053927;
    tri3_xyze(2, 0) = 0.84755282581460522628;
    nids.push_back(6871);
    tri3_xyze(0, 1) = 1.0005177142313128246e-14;
    tri3_xyze(1, 1) = 0.084549150281254134565;
    tri3_xyze(2, 1) = 0.84755282581460511526;
    nids.push_back(7121);
    tri3_xyze(0, 2) = -0.011334985551169712847;
    tri3_xyze(1, 2) = 0.086097763118054851383;
    tri3_xyze(2, 2) = 0.84833359417544340619;
    nids.push_back(-941);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005177142313128246e-14;
    tri3_xyze(1, 0) = 0.084549150281254134565;
    tri3_xyze(2, 0) = 0.84755282581460511526;
    nids.push_back(7121);
    tri3_xyze(0, 1) = 1.0005549544177918912e-14;
    tri3_xyze(1, 1) = 0.090630934270715279855;
    tri3_xyze(2, 1) = 0.8491143625362815861;
    nids.push_back(7119);
    tri3_xyze(0, 2) = -0.011334985551169712847;
    tri3_xyze(1, 2) = 0.086097763118054851383;
    tri3_xyze(2, 2) = 0.84833359417544340619;
    nids.push_back(-941);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005549544177918912e-14;
    tri3_xyze(1, 0) = 0.090630934270715279855;
    tri3_xyze(2, 0) = 0.8491143625362815861;
    nids.push_back(7119);
    tri3_xyze(0, 1) = -0.02345701186468549318;
    tri3_xyze(1, 1) = 0.087542760072790937187;
    tri3_xyze(2, 1) = 0.84911436253628169712;
    nids.push_back(6869);
    tri3_xyze(0, 2) = -0.011334985551169712847;
    tri3_xyze(1, 2) = 0.086097763118054851383;
    tri3_xyze(2, 2) = 0.84833359417544340619;
    nids.push_back(-941);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.02188293034001336998;
    tri3_xyze(1, 0) = 0.081668207847459053927;
    tri3_xyze(2, 0) = 0.84755282581460522628;
    nids.push_back(6871);
    tri3_xyze(0, 1) = -0.020371915026877046584;
    tri3_xyze(1, 1) = 0.076029021927819878668;
    tri3_xyze(2, 1) = 0.84524135262314858963;
    nids.push_back(6873);
    tri3_xyze(0, 2) = -0.010563711341717602066;
    tri3_xyze(1, 2) = 0.080239353869570234634;
    tri3_xyze(2, 2) = 0.84639708921887690796;
    nids.push_back(-942);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.020371915026877046584;
    tri3_xyze(1, 0) = 0.076029021927819878668;
    tri3_xyze(2, 0) = 0.84524135262314858963;
    nids.push_back(6873);
    tri3_xyze(0, 1) = 1.0004819660879340093e-14;
    tri3_xyze(1, 1) = 0.078711035421747829743;
    tri3_xyze(2, 1) = 0.84524135262314858963;
    nids.push_back(7123);
    tri3_xyze(0, 2) = -0.010563711341717602066;
    tri3_xyze(1, 2) = 0.080239353869570234634;
    tri3_xyze(2, 2) = 0.84639708921887690796;
    nids.push_back(-942);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004819660879340093e-14;
    tri3_xyze(1, 0) = 0.078711035421747829743;
    tri3_xyze(2, 0) = 0.84524135262314858963;
    nids.push_back(7123);
    tri3_xyze(0, 1) = 1.0005177142313128246e-14;
    tri3_xyze(1, 1) = 0.084549150281254134565;
    tri3_xyze(2, 1) = 0.84755282581460511526;
    nids.push_back(7121);
    tri3_xyze(0, 2) = -0.010563711341717602066;
    tri3_xyze(1, 2) = 0.080239353869570234634;
    tri3_xyze(2, 2) = 0.84639708921887690796;
    nids.push_back(-942);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005177142313128246e-14;
    tri3_xyze(1, 0) = 0.084549150281254134565;
    tri3_xyze(2, 0) = 0.84755282581460511526;
    nids.push_back(7121);
    tri3_xyze(0, 1) = -0.02188293034001336998;
    tri3_xyze(1, 1) = 0.081668207847459053927;
    tri3_xyze(2, 1) = 0.84755282581460522628;
    nids.push_back(6871);
    tri3_xyze(0, 2) = -0.010563711341717602066;
    tri3_xyze(1, 2) = 0.080239353869570234634;
    tri3_xyze(2, 2) = 0.84639708921887690796;
    nids.push_back(-942);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.020371915026877046584;
    tri3_xyze(1, 0) = 0.076029021927819878668;
    tri3_xyze(2, 0) = 0.84524135262314858963;
    nids.push_back(6873);
    tri3_xyze(0, 1) = -0.018947795539401655202;
    tri3_xyze(1, 1) = 0.070714135644512676748;
    tri3_xyze(2, 1) = 0.8422163962749485977;
    nids.push_back(6875);
    tri3_xyze(0, 2) = -0.009829927641564656024;
    tri3_xyze(1, 2) = 0.074665713311282966291;
    tri3_xyze(2, 2) = 0.84372887444904853815;
    nids.push_back(-943);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.018947795539401655202;
    tri3_xyze(1, 0) = 0.070714135644512676748;
    tri3_xyze(2, 0) = 0.8422163962749485977;
    nids.push_back(6875);
    tri3_xyze(0, 1) = 1.0069505089742456541e-14;
    tri3_xyze(1, 1) = 0.073208660251051493884;
    tri3_xyze(2, 1) = 0.84221639627494848668;
    nids.push_back(7125);
    tri3_xyze(0, 2) = -0.009829927641564656024;
    tri3_xyze(1, 2) = 0.074665713311282966291;
    tri3_xyze(2, 2) = 0.84372887444904853815;
    nids.push_back(-943);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0069505089742456541e-14;
    tri3_xyze(1, 0) = 0.073208660251051493884;
    tri3_xyze(2, 0) = 0.84221639627494848668;
    nids.push_back(7125);
    tri3_xyze(0, 1) = 1.0004819660879340093e-14;
    tri3_xyze(1, 1) = 0.078711035421747829743;
    tri3_xyze(2, 1) = 0.84524135262314858963;
    nids.push_back(7123);
    tri3_xyze(0, 2) = -0.009829927641564656024;
    tri3_xyze(1, 2) = 0.074665713311282966291;
    tri3_xyze(2, 2) = 0.84372887444904853815;
    nids.push_back(-943);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004819660879340093e-14;
    tri3_xyze(1, 0) = 0.078711035421747829743;
    tri3_xyze(2, 0) = 0.84524135262314858963;
    nids.push_back(7123);
    tri3_xyze(0, 1) = -0.020371915026877046584;
    tri3_xyze(1, 1) = 0.076029021927819878668;
    tri3_xyze(2, 1) = 0.84524135262314858963;
    nids.push_back(6873);
    tri3_xyze(0, 2) = -0.009829927641564656024;
    tri3_xyze(1, 2) = 0.074665713311282966291;
    tri3_xyze(2, 2) = 0.84372887444904853815;
    nids.push_back(-943);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.018947795539401655202;
    tri3_xyze(1, 0) = 0.070714135644512676748;
    tri3_xyze(2, 0) = 0.8422163962749485977;
    nids.push_back(6875);
    tri3_xyze(0, 1) = -0.017633031092632322351;
    tri3_xyze(1, 1) = 0.0658073679291842073;
    tri3_xyze(2, 1) = 0.83852566213863743982;
    nids.push_back(6877);
    tri3_xyze(0, 2) = -0.0091452066580034758331;
    tri3_xyze(1, 2) = 0.069464741084328790466;
    tri3_xyze(2, 2) = 0.84037102920679296325;
    nids.push_back(-944);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.017633031092632322351;
    tri3_xyze(1, 0) = 0.0658073679291842073;
    tri3_xyze(2, 0) = 0.83852566213863743982;
    nids.push_back(6877);
    tri3_xyze(0, 1) = 1.0004171685873872409e-14;
    tri3_xyze(1, 1) = 0.068128800512566770053;
    tri3_xyze(2, 1) = 0.83852566213863743982;
    nids.push_back(7127);
    tri3_xyze(0, 2) = -0.0091452066580034758331;
    tri3_xyze(1, 2) = 0.069464741084328790466;
    tri3_xyze(2, 2) = 0.84037102920679296325;
    nids.push_back(-944);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004171685873872409e-14;
    tri3_xyze(1, 0) = 0.068128800512566770053;
    tri3_xyze(2, 0) = 0.83852566213863743982;
    nids.push_back(7127);
    tri3_xyze(0, 1) = 1.0069505089742456541e-14;
    tri3_xyze(1, 1) = 0.073208660251051493884;
    tri3_xyze(2, 1) = 0.84221639627494848668;
    nids.push_back(7125);
    tri3_xyze(0, 2) = -0.0091452066580034758331;
    tri3_xyze(1, 2) = 0.069464741084328790466;
    tri3_xyze(2, 2) = 0.84037102920679296325;
    nids.push_back(-944);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0069505089742456541e-14;
    tri3_xyze(1, 0) = 0.073208660251051493884;
    tri3_xyze(2, 0) = 0.84221639627494848668;
    nids.push_back(7125);
    tri3_xyze(0, 1) = -0.018947795539401655202;
    tri3_xyze(1, 1) = 0.070714135644512676748;
    tri3_xyze(2, 1) = 0.8422163962749485977;
    nids.push_back(6875);
    tri3_xyze(0, 2) = -0.0091452066580034758331;
    tri3_xyze(1, 2) = 0.069464741084328790466;
    tri3_xyze(2, 2) = 0.84037102920679296325;
    nids.push_back(-944);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.017633031092632322351;
    tri3_xyze(1, 0) = 0.0658073679291842073;
    tri3_xyze(2, 0) = 0.83852566213863743982;
    nids.push_back(6877);
    tri3_xyze(0, 1) = -0.016448356307296806839;
    tri3_xyze(1, 1) = 0.061386101439866169738;
    tri3_xyze(2, 1) = 0.83422735529628255779;
    nids.push_back(6879);
    tri3_xyze(0, 2) = -0.0085203468499772654771;
    tri3_xyze(1, 2) = 0.064718459627636906739;
    tri3_xyze(2, 2) = 0.8363765087174599433;
    nids.push_back(-945);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.016448356307296806839;
    tri3_xyze(1, 0) = 0.061386101439866169738;
    tri3_xyze(2, 0) = 0.83422735529628255779;
    nids.push_back(6879);
    tri3_xyze(0, 1) = 1.0060336543049417096e-14;
    tri3_xyze(1, 1) = 0.063551568628930479865;
    tri3_xyze(2, 1) = 0.83422735529628244677;
    nids.push_back(7129);
    tri3_xyze(0, 2) = -0.0085203468499772654771;
    tri3_xyze(1, 2) = 0.064718459627636906739;
    tri3_xyze(2, 2) = 0.8363765087174599433;
    nids.push_back(-945);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0060336543049417096e-14;
    tri3_xyze(1, 0) = 0.063551568628930479865;
    tri3_xyze(2, 0) = 0.83422735529628244677;
    nids.push_back(7129);
    tri3_xyze(0, 1) = 1.0004171685873872409e-14;
    tri3_xyze(1, 1) = 0.068128800512566770053;
    tri3_xyze(2, 1) = 0.83852566213863743982;
    nids.push_back(7127);
    tri3_xyze(0, 2) = -0.0085203468499772654771;
    tri3_xyze(1, 2) = 0.064718459627636906739;
    tri3_xyze(2, 2) = 0.8363765087174599433;
    nids.push_back(-945);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004171685873872409e-14;
    tri3_xyze(1, 0) = 0.068128800512566770053;
    tri3_xyze(2, 0) = 0.83852566213863743982;
    nids.push_back(7127);
    tri3_xyze(0, 1) = -0.017633031092632322351;
    tri3_xyze(1, 1) = 0.0658073679291842073;
    tri3_xyze(2, 1) = 0.83852566213863743982;
    nids.push_back(6877);
    tri3_xyze(0, 2) = -0.0085203468499772654771;
    tri3_xyze(1, 2) = 0.064718459627636906739;
    tri3_xyze(2, 2) = 0.8363765087174599433;
    nids.push_back(-945);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.016448356307296806839;
    tri3_xyze(1, 0) = 0.061386101439866169738;
    tri3_xyze(2, 0) = 0.83422735529628255779;
    nids.push_back(6879);
    tri3_xyze(0, 1) = -0.015412454212450349758;
    tri3_xyze(1, 1) = 0.057520062190231820154;
    tri3_xyze(2, 1) = 0.82938926261447198929;
    nids.push_back(6881);
    tri3_xyze(0, 2) = -0.0079652026299317467417;
    tri3_xyze(1, 2) = 0.060501720635070505017;
    tri3_xyze(2, 2) = 0.83180830895537716252;
    nids.push_back(-946);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.015412454212450349758;
    tri3_xyze(1, 0) = 0.057520062190231820154;
    tri3_xyze(2, 0) = 0.82938926261447198929;
    nids.push_back(6881);
    tri3_xyze(0, 1) = 1.010942687419677122e-14;
    tri3_xyze(1, 1) = 0.059549150281253564188;
    tri3_xyze(2, 1) = 0.82938926261447187827;
    nids.push_back(7131);
    tri3_xyze(0, 2) = -0.0079652026299317467417;
    tri3_xyze(1, 2) = 0.060501720635070505017;
    tri3_xyze(2, 2) = 0.83180830895537716252;
    nids.push_back(-946);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.010942687419677122e-14;
    tri3_xyze(1, 0) = 0.059549150281253564188;
    tri3_xyze(2, 0) = 0.82938926261447187827;
    nids.push_back(7131);
    tri3_xyze(0, 1) = 1.0060336543049417096e-14;
    tri3_xyze(1, 1) = 0.063551568628930479865;
    tri3_xyze(2, 1) = 0.83422735529628244677;
    nids.push_back(7129);
    tri3_xyze(0, 2) = -0.0079652026299317467417;
    tri3_xyze(1, 2) = 0.060501720635070505017;
    tri3_xyze(2, 2) = 0.83180830895537716252;
    nids.push_back(-946);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0060336543049417096e-14;
    tri3_xyze(1, 0) = 0.063551568628930479865;
    tri3_xyze(2, 0) = 0.83422735529628244677;
    nids.push_back(7129);
    tri3_xyze(0, 1) = -0.016448356307296806839;
    tri3_xyze(1, 1) = 0.061386101439866169738;
    tri3_xyze(2, 1) = 0.83422735529628255779;
    nids.push_back(6879);
    tri3_xyze(0, 2) = -0.0079652026299317467417;
    tri3_xyze(1, 2) = 0.060501720635070505017;
    tri3_xyze(2, 2) = 0.83180830895537716252;
    nids.push_back(-946);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.015412454212450349758;
    tri3_xyze(1, 0) = 0.057520062190231820154;
    tri3_xyze(2, 0) = 0.82938926261447198929;
    nids.push_back(6881);
    tri3_xyze(0, 1) = -0.014541661602946407209;
    tri3_xyze(1, 1) = 0.054270219928707633117;
    tri3_xyze(2, 1) = 0.82408768370493412636;
    nids.push_back(6883);
    tri3_xyze(0, 2) = -0.0074885289538441611457;
    tri3_xyze(1, 2) = 0.056881024599500160099;
    tri3_xyze(2, 2) = 0.82673847315970305782;
    nids.push_back(-947);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.014541661602946407209;
    tri3_xyze(1, 0) = 0.054270219928707633117;
    tri3_xyze(2, 0) = 0.82408768370493412636;
    nids.push_back(6883);
    tri3_xyze(0, 1) = 1.0003440318568768715e-14;
    tri3_xyze(1, 1) = 0.056184665997807595184;
    tri3_xyze(2, 1) = 0.82408768370493412636;
    nids.push_back(7133);
    tri3_xyze(0, 2) = -0.0074885289538441611457;
    tri3_xyze(1, 2) = 0.056881024599500160099;
    tri3_xyze(2, 2) = 0.82673847315970305782;
    nids.push_back(-947);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003440318568768715e-14;
    tri3_xyze(1, 0) = 0.056184665997807595184;
    tri3_xyze(2, 0) = 0.82408768370493412636;
    nids.push_back(7133);
    tri3_xyze(0, 1) = 1.010942687419677122e-14;
    tri3_xyze(1, 1) = 0.059549150281253564188;
    tri3_xyze(2, 1) = 0.82938926261447187827;
    nids.push_back(7131);
    tri3_xyze(0, 2) = -0.0074885289538441611457;
    tri3_xyze(1, 2) = 0.056881024599500160099;
    tri3_xyze(2, 2) = 0.82673847315970305782;
    nids.push_back(-947);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.010942687419677122e-14;
    tri3_xyze(1, 0) = 0.059549150281253564188;
    tri3_xyze(2, 0) = 0.82938926261447187827;
    nids.push_back(7131);
    tri3_xyze(0, 1) = -0.015412454212450349758;
    tri3_xyze(1, 1) = 0.057520062190231820154;
    tri3_xyze(2, 1) = 0.82938926261447198929;
    nids.push_back(6881);
    tri3_xyze(0, 2) = -0.0074885289538441611457;
    tri3_xyze(1, 2) = 0.056881024599500160099;
    tri3_xyze(2, 2) = 0.82673847315970305782;
    nids.push_back(-947);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.014541661602946407209;
    tri3_xyze(1, 0) = 0.054270219928707633117;
    tri3_xyze(2, 0) = 0.82408768370493412636;
    nids.push_back(6883);
    tri3_xyze(0, 1) = -0.013849711398423282874;
    tri3_xyze(1, 1) = 0.051687826609119663646;
    tri3_xyze(2, 1) = 0.81840622763408232831;
    nids.push_back(6885);
    tri3_xyze(0, 2) = -0.0070978432503374213128;
    tri3_xyze(1, 2) = 0.053913472060305736577;
    tri3_xyze(2, 2) = 0.82124695566950822734;
    nids.push_back(-948);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.013849711398423282874;
    tri3_xyze(1, 0) = 0.051687826609119663646;
    tri3_xyze(2, 0) = 0.81840622763408232831;
    nids.push_back(6885);
    tri3_xyze(0, 1) = 1.0003276614502322774e-14;
    tri3_xyze(1, 1) = 0.053511175705588047424;
    tri3_xyze(2, 1) = 0.81840622763408232831;
    nids.push_back(7135);
    tri3_xyze(0, 2) = -0.0070978432503374213128;
    tri3_xyze(1, 2) = 0.053913472060305736577;
    tri3_xyze(2, 2) = 0.82124695566950822734;
    nids.push_back(-948);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003276614502322774e-14;
    tri3_xyze(1, 0) = 0.053511175705588047424;
    tri3_xyze(2, 0) = 0.81840622763408232831;
    nids.push_back(7135);
    tri3_xyze(0, 1) = 1.0003440318568768715e-14;
    tri3_xyze(1, 1) = 0.056184665997807595184;
    tri3_xyze(2, 1) = 0.82408768370493412636;
    nids.push_back(7133);
    tri3_xyze(0, 2) = -0.0070978432503374213128;
    tri3_xyze(1, 2) = 0.053913472060305736577;
    tri3_xyze(2, 2) = 0.82124695566950822734;
    nids.push_back(-948);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003440318568768715e-14;
    tri3_xyze(1, 0) = 0.056184665997807595184;
    tri3_xyze(2, 0) = 0.82408768370493412636;
    nids.push_back(7133);
    tri3_xyze(0, 1) = -0.014541661602946407209;
    tri3_xyze(1, 1) = 0.054270219928707633117;
    tri3_xyze(2, 1) = 0.82408768370493412636;
    nids.push_back(6883);
    tri3_xyze(0, 2) = -0.0070978432503374213128;
    tri3_xyze(1, 2) = 0.053913472060305736577;
    tri3_xyze(2, 2) = 0.82124695566950822734;
    nids.push_back(-948);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.013849711398423282874;
    tri3_xyze(1, 0) = 0.051687826609119663646;
    tri3_xyze(2, 0) = 0.81840622763408232831;
    nids.push_back(6885);
    tri3_xyze(0, 1) = -0.013347516066957392997;
    tri3_xyze(1, 1) = 0.04981360811676469047;
    tri3_xyze(2, 1) = 0.81243449435809123038;
    nids.push_back(6887);
    tri3_xyze(0, 2) = -0.0067993068663401677598;
    tri3_xyze(1, 2) = 0.051645863093760312201;
    tri3_xyze(2, 2) = 0.81542036099608672384;
    nids.push_back(-949);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.013347516066957392997;
    tri3_xyze(1, 0) = 0.04981360811676469047;
    tri3_xyze(2, 0) = 0.81243449435809123038;
    nids.push_back(6887);
    tri3_xyze(0, 1) = 1.0003157803325775479e-14;
    tri3_xyze(1, 1) = 0.051570841943568840327;
    tri3_xyze(2, 1) = 0.81243449435809111936;
    nids.push_back(7137);
    tri3_xyze(0, 2) = -0.0067993068663401677598;
    tri3_xyze(1, 2) = 0.051645863093760312201;
    tri3_xyze(2, 2) = 0.81542036099608672384;
    nids.push_back(-949);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003157803325775479e-14;
    tri3_xyze(1, 0) = 0.051570841943568840327;
    tri3_xyze(2, 0) = 0.81243449435809111936;
    nids.push_back(7137);
    tri3_xyze(0, 1) = 1.0003276614502322774e-14;
    tri3_xyze(1, 1) = 0.053511175705588047424;
    tri3_xyze(2, 1) = 0.81840622763408232831;
    nids.push_back(7135);
    tri3_xyze(0, 2) = -0.0067993068663401677598;
    tri3_xyze(1, 2) = 0.051645863093760312201;
    tri3_xyze(2, 2) = 0.81542036099608672384;
    nids.push_back(-949);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003276614502322774e-14;
    tri3_xyze(1, 0) = 0.053511175705588047424;
    tri3_xyze(2, 0) = 0.81840622763408232831;
    nids.push_back(7135);
    tri3_xyze(0, 1) = -0.013849711398423282874;
    tri3_xyze(1, 1) = 0.051687826609119663646;
    tri3_xyze(2, 1) = 0.81840622763408232831;
    nids.push_back(6885);
    tri3_xyze(0, 2) = -0.0067993068663401677598;
    tri3_xyze(1, 2) = 0.051645863093760312201;
    tri3_xyze(2, 2) = 0.81542036099608672384;
    nids.push_back(-949);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.013347516066957392997;
    tri3_xyze(1, 0) = 0.04981360811676469047;
    tri3_xyze(2, 0) = 0.81243449435809123038;
    nids.push_back(6887);
    tri3_xyze(0, 1) = -0.013042995528922771836;
    tri3_xyze(1, 1) = 0.048677121996871075549;
    tri3_xyze(2, 1) = 0.80626666167806371188;
    nids.push_back(6889);
    tri3_xyze(0, 2) = -0.0065976278989650395668;
    tri3_xyze(1, 2) = 0.050113959247870230873;
    tri3_xyze(2, 2) = 0.80935057801807741562;
    nids.push_back(-950);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.013042995528922771836;
    tri3_xyze(1, 0) = 0.048677121996871075549;
    tri3_xyze(2, 0) = 0.80626666167806371188;
    nids.push_back(6889);
    tri3_xyze(0, 1) = 1.0003085758762357486e-14;
    tri3_xyze(1, 1) = 0.050394264934276310208;
    tri3_xyze(2, 1) = 0.80626666167806371188;
    nids.push_back(7139);
    tri3_xyze(0, 2) = -0.0065976278989650395668;
    tri3_xyze(1, 2) = 0.050113959247870230873;
    tri3_xyze(2, 2) = 0.80935057801807741562;
    nids.push_back(-950);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003085758762357486e-14;
    tri3_xyze(1, 0) = 0.050394264934276310208;
    tri3_xyze(2, 0) = 0.80626666167806371188;
    nids.push_back(7139);
    tri3_xyze(0, 1) = 1.0003157803325775479e-14;
    tri3_xyze(1, 1) = 0.051570841943568840327;
    tri3_xyze(2, 1) = 0.81243449435809111936;
    nids.push_back(7137);
    tri3_xyze(0, 2) = -0.0065976278989650395668;
    tri3_xyze(1, 2) = 0.050113959247870230873;
    tri3_xyze(2, 2) = 0.80935057801807741562;
    nids.push_back(-950);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003157803325775479e-14;
    tri3_xyze(1, 0) = 0.051570841943568840327;
    tri3_xyze(2, 0) = 0.81243449435809111936;
    nids.push_back(7137);
    tri3_xyze(0, 1) = -0.013347516066957392997;
    tri3_xyze(1, 1) = 0.04981360811676469047;
    tri3_xyze(2, 1) = 0.81243449435809123038;
    nids.push_back(6887);
    tri3_xyze(0, 2) = -0.0065976278989650395668;
    tri3_xyze(1, 2) = 0.050113959247870230873;
    tri3_xyze(2, 2) = 0.80935057801807741562;
    nids.push_back(-950);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.9908151490063951676e-15;
    tri3_xyze(1, 0) = 0.050000000000000002776;
    tri3_xyze(2, 0) = 0.79999999999984849897;
    nids.push_back(7041);
    tri3_xyze(0, 1) = 1.0003085758762357486e-14;
    tri3_xyze(1, 1) = 0.050394264934275928569;
    tri3_xyze(2, 1) = 0.79373333832163317503;
    nids.push_back(7042);
    tri3_xyze(0, 2) = 0.0064959869460246974157;
    tri3_xyze(1, 2) = 0.049341919561400007788;
    tri3_xyze(2, 2) = 0.79686666916074089251;
    nids.push_back(-951);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003085758762357486e-14;
    tri3_xyze(1, 0) = 0.050394264934275928569;
    tri3_xyze(2, 0) = 0.79373333832163317503;
    nids.push_back(7042);
    tri3_xyze(0, 1) = 0.013042995528942778402;
    tri3_xyze(1, 1) = 0.048677121996870680032;
    tri3_xyze(2, 1) = 0.79373333832163328605;
    nids.push_back(7292);
    tri3_xyze(0, 2) = 0.0064959869460246974157;
    tri3_xyze(1, 2) = 0.049341919561400007788;
    tri3_xyze(2, 2) = 0.79686666916074089251;
    nids.push_back(-951);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.012940952255136016838;
    tri3_xyze(1, 0) = 0.048296291314453426713;
    tri3_xyze(2, 0) = 0.79999999999984849897;
    nids.push_back(7291);
    tri3_xyze(0, 1) = 9.9908151490063951676e-15;
    tri3_xyze(1, 1) = 0.050000000000000002776;
    tri3_xyze(2, 1) = 0.79999999999984849897;
    nids.push_back(7041);
    tri3_xyze(0, 2) = 0.0064959869460246974157;
    tri3_xyze(1, 2) = 0.049341919561400007788;
    tri3_xyze(2, 2) = 0.79686666916074089251;
    nids.push_back(-951);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003085758762357486e-14;
    tri3_xyze(1, 0) = 0.050394264934276310208;
    tri3_xyze(2, 0) = 0.80626666167806371188;
    nids.push_back(7139);
    tri3_xyze(0, 1) = 9.9908151490063951676e-15;
    tri3_xyze(1, 1) = 0.050000000000000002776;
    tri3_xyze(2, 1) = 0.79999999999984849897;
    nids.push_back(7041);
    tri3_xyze(0, 2) = 0.0064959869460247086914;
    tri3_xyze(1, 2) = 0.049341919561400195138;
    tri3_xyze(2, 2) = 0.80313333083895610542;
    nids.push_back(-952);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.9908151490063951676e-15;
    tri3_xyze(1, 0) = 0.050000000000000002776;
    tri3_xyze(2, 0) = 0.79999999999984849897;
    nids.push_back(7041);
    tri3_xyze(0, 1) = 0.012940952255136016838;
    tri3_xyze(1, 1) = 0.048296291314453426713;
    tri3_xyze(2, 1) = 0.79999999999984849897;
    nids.push_back(7291);
    tri3_xyze(0, 2) = 0.0064959869460247086914;
    tri3_xyze(1, 2) = 0.049341919561400195138;
    tri3_xyze(2, 2) = 0.80313333083895610542;
    nids.push_back(-952);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.012940952255136016838;
    tri3_xyze(1, 0) = 0.048296291314453426713;
    tri3_xyze(2, 0) = 0.79999999999984849897;
    nids.push_back(7291);
    tri3_xyze(0, 1) = 0.013042995528942823505;
    tri3_xyze(1, 1) = 0.048677121996871054732;
    tri3_xyze(2, 1) = 0.80626666167806371188;
    nids.push_back(7389);
    tri3_xyze(0, 2) = 0.0064959869460247086914;
    tri3_xyze(1, 2) = 0.049341919561400195138;
    tri3_xyze(2, 2) = 0.80313333083895610542;
    nids.push_back(-952);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013042995528942823505;
    tri3_xyze(1, 0) = 0.048677121996871054732;
    tri3_xyze(2, 0) = 0.80626666167806371188;
    nids.push_back(7389);
    tri3_xyze(0, 1) = 1.0003085758762357486e-14;
    tri3_xyze(1, 1) = 0.050394264934276310208;
    tri3_xyze(2, 1) = 0.80626666167806371188;
    nids.push_back(7139);
    tri3_xyze(0, 2) = 0.0064959869460247086914;
    tri3_xyze(1, 2) = 0.049341919561400195138;
    tri3_xyze(2, 2) = 0.80313333083895610542;
    nids.push_back(-952);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003085758762357486e-14;
    tri3_xyze(1, 0) = 0.050394264934275928569;
    tri3_xyze(2, 0) = 0.79373333832163317503;
    nids.push_back(7042);
    tri3_xyze(0, 1) = 1.0003157803325775479e-14;
    tri3_xyze(1, 1) = 0.051570841943568049293;
    tri3_xyze(2, 1) = 0.78756550564160576755;
    nids.push_back(7045);
    tri3_xyze(0, 2) = 0.0065976278989850478673;
    tri3_xyze(1, 2) = 0.050113959247869648006;
    tri3_xyze(2, 2) = 0.79064942198161947129;
    nids.push_back(-953);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003157803325775479e-14;
    tri3_xyze(1, 0) = 0.051570841943568049293;
    tri3_xyze(2, 0) = 0.78756550564160576755;
    nids.push_back(7045);
    tri3_xyze(0, 1) = 0.013347516066977408236;
    tri3_xyze(1, 1) = 0.04981360811676393413;
    tri3_xyze(2, 1) = 0.78756550564160565653;
    nids.push_back(7295);
    tri3_xyze(0, 2) = 0.0065976278989850478673;
    tri3_xyze(1, 2) = 0.050113959247869648006;
    tri3_xyze(2, 2) = 0.79064942198161947129;
    nids.push_back(-953);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013347516066977408236;
    tri3_xyze(1, 0) = 0.04981360811676393413;
    tri3_xyze(2, 0) = 0.78756550564160565653;
    nids.push_back(7295);
    tri3_xyze(0, 1) = 0.013042995528942778402;
    tri3_xyze(1, 1) = 0.048677121996870680032;
    tri3_xyze(2, 1) = 0.79373333832163328605;
    nids.push_back(7292);
    tri3_xyze(0, 2) = 0.0065976278989850478673;
    tri3_xyze(1, 2) = 0.050113959247869648006;
    tri3_xyze(2, 2) = 0.79064942198161947129;
    nids.push_back(-953);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013042995528942778402;
    tri3_xyze(1, 0) = 0.048677121996870680032;
    tri3_xyze(2, 0) = 0.79373333832163328605;
    nids.push_back(7292);
    tri3_xyze(0, 1) = 1.0003085758762357486e-14;
    tri3_xyze(1, 1) = 0.050394264934275928569;
    tri3_xyze(2, 1) = 0.79373333832163317503;
    nids.push_back(7042);
    tri3_xyze(0, 2) = 0.0065976278989850478673;
    tri3_xyze(1, 2) = 0.050113959247869648006;
    tri3_xyze(2, 2) = 0.79064942198161947129;
    nids.push_back(-953);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003157803325775479e-14;
    tri3_xyze(1, 0) = 0.051570841943568049293;
    tri3_xyze(2, 0) = 0.78756550564160576755;
    nids.push_back(7045);
    tri3_xyze(0, 1) = 9.9557491430278385319e-15;
    tri3_xyze(1, 1) = 0.053511175705586895568;
    tri3_xyze(2, 1) = 0.78159377236561444757;
    nids.push_back(7047);
    tri3_xyze(0, 2) = 0.0067993068663601777951;
    tri3_xyze(1, 2) = 0.051645863093759347695;
    tri3_xyze(2, 2) = 0.78457963900361016307;
    nids.push_back(-954);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.9557491430278385319e-15;
    tri3_xyze(1, 0) = 0.053511175705586895568;
    tri3_xyze(2, 0) = 0.78159377236561444757;
    nids.push_back(7047);
    tri3_xyze(0, 1) = 0.013849711398443344951;
    tri3_xyze(1, 1) = 0.051687826609118511789;
    tri3_xyze(2, 1) = 0.78159377236561444757;
    nids.push_back(7297);
    tri3_xyze(0, 2) = 0.0067993068663601777951;
    tri3_xyze(1, 2) = 0.051645863093759347695;
    tri3_xyze(2, 2) = 0.78457963900361016307;
    nids.push_back(-954);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013849711398443344951;
    tri3_xyze(1, 0) = 0.051687826609118511789;
    tri3_xyze(2, 0) = 0.78159377236561444757;
    nids.push_back(7297);
    tri3_xyze(0, 1) = 0.013347516066977408236;
    tri3_xyze(1, 1) = 0.04981360811676393413;
    tri3_xyze(2, 1) = 0.78756550564160565653;
    nids.push_back(7295);
    tri3_xyze(0, 2) = 0.0067993068663601777951;
    tri3_xyze(1, 2) = 0.051645863093759347695;
    tri3_xyze(2, 2) = 0.78457963900361016307;
    nids.push_back(-954);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013347516066977408236;
    tri3_xyze(1, 0) = 0.04981360811676393413;
    tri3_xyze(2, 0) = 0.78756550564160565653;
    nids.push_back(7295);
    tri3_xyze(0, 1) = 1.0003157803325775479e-14;
    tri3_xyze(1, 1) = 0.051570841943568049293;
    tri3_xyze(2, 1) = 0.78756550564160576755;
    nids.push_back(7045);
    tri3_xyze(0, 2) = 0.0067993068663601777951;
    tri3_xyze(1, 2) = 0.051645863093759347695;
    tri3_xyze(2, 2) = 0.78457963900361016307;
    nids.push_back(-954);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.9557491430278385319e-15;
    tri3_xyze(1, 0) = 0.053511175705586895568;
    tri3_xyze(2, 0) = 0.78159377236561444757;
    nids.push_back(7047);
    tri3_xyze(0, 1) = 1.0053342326426079629e-14;
    tri3_xyze(1, 1) = 0.056184665997806082505;
    tri3_xyze(2, 1) = 0.77591231629476242748;
    nids.push_back(7049);
    tri3_xyze(0, 2) = 0.0070978432503574556342;
    tri3_xyze(1, 2) = 0.05391347206030440431;
    tri3_xyze(2, 2) = 0.77875304433018854855;
    nids.push_back(-955);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0053342326426079629e-14;
    tri3_xyze(1, 0) = 0.056184665997806082505;
    tri3_xyze(2, 0) = 0.77591231629476242748;
    nids.push_back(7049);
    tri3_xyze(0, 1) = 0.014541661602966469286;
    tri3_xyze(1, 1) = 0.054270219928706113499;
    tri3_xyze(2, 1) = 0.7759123162947625385;
    nids.push_back(7299);
    tri3_xyze(0, 2) = 0.0070978432503574556342;
    tri3_xyze(1, 2) = 0.05391347206030440431;
    tri3_xyze(2, 2) = 0.77875304433018854855;
    nids.push_back(-955);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.014541661602966469286;
    tri3_xyze(1, 0) = 0.054270219928706113499;
    tri3_xyze(2, 0) = 0.7759123162947625385;
    nids.push_back(7299);
    tri3_xyze(0, 1) = 0.013849711398443344951;
    tri3_xyze(1, 1) = 0.051687826609118511789;
    tri3_xyze(2, 1) = 0.78159377236561444757;
    nids.push_back(7297);
    tri3_xyze(0, 2) = 0.0070978432503574556342;
    tri3_xyze(1, 2) = 0.05391347206030440431;
    tri3_xyze(2, 2) = 0.77875304433018854855;
    nids.push_back(-955);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013849711398443344951;
    tri3_xyze(1, 0) = 0.051687826609118511789;
    tri3_xyze(2, 0) = 0.78159377236561444757;
    nids.push_back(7297);
    tri3_xyze(0, 1) = 9.9557491430278385319e-15;
    tri3_xyze(1, 1) = 0.053511175705586895568;
    tri3_xyze(2, 1) = 0.78159377236561444757;
    nids.push_back(7047);
    tri3_xyze(0, 2) = 0.0070978432503574556342;
    tri3_xyze(1, 2) = 0.05391347206030440431;
    tri3_xyze(2, 2) = 0.77875304433018854855;
    nids.push_back(-955);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0053342326426079629e-14;
    tri3_xyze(1, 0) = 0.056184665997806082505;
    tri3_xyze(2, 0) = 0.77591231629476242748;
    nids.push_back(7049);
    tri3_xyze(0, 1) = 1.0056536604005482796e-14;
    tri3_xyze(1, 1) = 0.059549150281251718442;
    tri3_xyze(2, 1) = 0.77061073738522445353;
    nids.push_back(7051);
    tri3_xyze(0, 2) = 0.0074885289538642197532;
    tri3_xyze(1, 2) = 0.05688102459949846007;
    tri3_xyze(2, 2) = 0.773261526839993385;
    nids.push_back(-956);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0056536604005482796e-14;
    tri3_xyze(1, 0) = 0.059549150281251718442;
    tri3_xyze(2, 0) = 0.77061073738522445353;
    nids.push_back(7051);
    tri3_xyze(0, 1) = 0.015412454212470300813;
    tri3_xyze(1, 1) = 0.057520062190229932775;
    tri3_xyze(2, 1) = 0.77061073738522456456;
    nids.push_back(7301);
    tri3_xyze(0, 2) = 0.0074885289538642197532;
    tri3_xyze(1, 2) = 0.05688102459949846007;
    tri3_xyze(2, 2) = 0.773261526839993385;
    nids.push_back(-956);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.015412454212470300813;
    tri3_xyze(1, 0) = 0.057520062190229932775;
    tri3_xyze(2, 0) = 0.77061073738522456456;
    nids.push_back(7301);
    tri3_xyze(0, 1) = 0.014541661602966469286;
    tri3_xyze(1, 1) = 0.054270219928706113499;
    tri3_xyze(2, 1) = 0.7759123162947625385;
    nids.push_back(7299);
    tri3_xyze(0, 2) = 0.0074885289538642197532;
    tri3_xyze(1, 2) = 0.05688102459949846007;
    tri3_xyze(2, 2) = 0.773261526839993385;
    nids.push_back(-956);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.014541661602966469286;
    tri3_xyze(1, 0) = 0.054270219928706113499;
    tri3_xyze(2, 0) = 0.7759123162947625385;
    nids.push_back(7299);
    tri3_xyze(0, 1) = 1.0053342326426079629e-14;
    tri3_xyze(1, 1) = 0.056184665997806082505;
    tri3_xyze(2, 1) = 0.77591231629476242748;
    nids.push_back(7049);
    tri3_xyze(0, 2) = 0.0074885289538642197532;
    tri3_xyze(1, 2) = 0.05688102459949846007;
    tri3_xyze(2, 2) = 0.773261526839993385;
    nids.push_back(-956);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0056536604005482796e-14;
    tri3_xyze(1, 0) = 0.059549150281251718442;
    tri3_xyze(2, 0) = 0.77061073738522445353;
    nids.push_back(7051);
    tri3_xyze(0, 1) = 1.0060336543049417096e-14;
    tri3_xyze(1, 1) = 0.063551568628928342686;
    tri3_xyze(2, 1) = 0.76577264470341366298;
    nids.push_back(7053);
    tri3_xyze(0, 2) = 0.0079652026299518140229;
    tri3_xyze(1, 2) = 0.060501720635068520493;
    tri3_xyze(2, 2) = 0.76819169104431905826;
    nids.push_back(-957);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0060336543049417096e-14;
    tri3_xyze(1, 0) = 0.063551568628928342686;
    tri3_xyze(2, 0) = 0.76577264470341366298;
    nids.push_back(7053);
    tri3_xyze(0, 1) = 0.016448356307316835956;
    tri3_xyze(1, 1) = 0.061386101439864095008;
    tri3_xyze(2, 1) = 0.76577264470341366298;
    nids.push_back(7303);
    tri3_xyze(0, 2) = 0.0079652026299518140229;
    tri3_xyze(1, 2) = 0.060501720635068520493;
    tri3_xyze(2, 2) = 0.76819169104431905826;
    nids.push_back(-957);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.016448356307316835956;
    tri3_xyze(1, 0) = 0.061386101439864095008;
    tri3_xyze(2, 0) = 0.76577264470341366298;
    nids.push_back(7303);
    tri3_xyze(0, 1) = 0.015412454212470300813;
    tri3_xyze(1, 1) = 0.057520062190229932775;
    tri3_xyze(2, 1) = 0.77061073738522456456;
    nids.push_back(7301);
    tri3_xyze(0, 2) = 0.0079652026299518140229;
    tri3_xyze(1, 2) = 0.060501720635068520493;
    tri3_xyze(2, 2) = 0.76819169104431905826;
    nids.push_back(-957);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.015412454212470300813;
    tri3_xyze(1, 0) = 0.057520062190229932775;
    tri3_xyze(2, 0) = 0.77061073738522456456;
    nids.push_back(7301);
    tri3_xyze(0, 1) = 1.0056536604005482796e-14;
    tri3_xyze(1, 1) = 0.059549150281251718442;
    tri3_xyze(2, 1) = 0.77061073738522445353;
    nids.push_back(7051);
    tri3_xyze(0, 2) = 0.0079652026299518140229;
    tri3_xyze(1, 2) = 0.060501720635068520493;
    tri3_xyze(2, 2) = 0.76819169104431905826;
    nids.push_back(-957);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0060336543049417096e-14;
    tri3_xyze(1, 0) = 0.063551568628928342686;
    tri3_xyze(2, 0) = 0.76577264470341366298;
    nids.push_back(7053);
    tri3_xyze(0, 1) = 1.0004171685873872409e-14;
    tri3_xyze(1, 1) = 0.068128800512564396952;
    tri3_xyze(2, 1) = 0.76147433786105844788;
    nids.push_back(7055);
    tri3_xyze(0, 2) = 0.0085203468499973015332;
    tri3_xyze(1, 2) = 0.064718459627634630782;
    tri3_xyze(2, 2) = 0.76362349128223616646;
    nids.push_back(-958);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004171685873872409e-14;
    tri3_xyze(1, 0) = 0.068128800512564396952;
    tri3_xyze(2, 0) = 0.76147433786105844788;
    nids.push_back(7055);
    tri3_xyze(0, 1) = 0.017633031092652306365;
    tri3_xyze(1, 1) = 0.065807367929181695421;
    tri3_xyze(2, 1) = 0.76147433786105855891;
    nids.push_back(7305);
    tri3_xyze(0, 2) = 0.0085203468499973015332;
    tri3_xyze(1, 2) = 0.064718459627634630782;
    tri3_xyze(2, 2) = 0.76362349128223616646;
    nids.push_back(-958);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.017633031092652306365;
    tri3_xyze(1, 0) = 0.065807367929181695421;
    tri3_xyze(2, 0) = 0.76147433786105855891;
    nids.push_back(7305);
    tri3_xyze(0, 1) = 0.016448356307316835956;
    tri3_xyze(1, 1) = 0.061386101439864095008;
    tri3_xyze(2, 1) = 0.76577264470341366298;
    nids.push_back(7303);
    tri3_xyze(0, 2) = 0.0085203468499973015332;
    tri3_xyze(1, 2) = 0.064718459627634630782;
    tri3_xyze(2, 2) = 0.76362349128223616646;
    nids.push_back(-958);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.016448356307316835956;
    tri3_xyze(1, 0) = 0.061386101439864095008;
    tri3_xyze(2, 0) = 0.76577264470341366298;
    nids.push_back(7303);
    tri3_xyze(0, 1) = 1.0060336543049417096e-14;
    tri3_xyze(1, 1) = 0.063551568628928342686;
    tri3_xyze(2, 1) = 0.76577264470341366298;
    nids.push_back(7053);
    tri3_xyze(0, 2) = 0.0085203468499973015332;
    tri3_xyze(1, 2) = 0.064718459627634630782;
    tri3_xyze(2, 2) = 0.76362349128223616646;
    nids.push_back(-958);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004171685873872409e-14;
    tri3_xyze(1, 0) = 0.068128800512564396952;
    tri3_xyze(2, 0) = 0.76147433786105844788;
    nids.push_back(7055);
    tri3_xyze(0, 1) = 1.0004482737572315978e-14;
    tri3_xyze(1, 1) = 0.073208660251048829348;
    tri3_xyze(2, 1) = 0.75778360372474695694;
    nids.push_back(7057);
    tri3_xyze(0, 2) = 0.0091452066580234858684;
    tri3_xyze(1, 2) = 0.069464741084326264708;
    tri3_xyze(2, 2) = 0.75962897079290270241;
    nids.push_back(-959);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004482737572315978e-14;
    tri3_xyze(1, 0) = 0.073208660251048829348;
    tri3_xyze(2, 0) = 0.75778360372474695694;
    nids.push_back(7057);
    tri3_xyze(0, 1) = 0.018947795539421625338;
    tri3_xyze(1, 1) = 0.070714135644510123235;
    tri3_xyze(2, 1) = 0.75778360372474695694;
    nids.push_back(7307);
    tri3_xyze(0, 2) = 0.0091452066580234858684;
    tri3_xyze(1, 2) = 0.069464741084326264708;
    tri3_xyze(2, 2) = 0.75962897079290270241;
    nids.push_back(-959);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.018947795539421625338;
    tri3_xyze(1, 0) = 0.070714135644510123235;
    tri3_xyze(2, 0) = 0.75778360372474695694;
    nids.push_back(7307);
    tri3_xyze(0, 1) = 0.017633031092652306365;
    tri3_xyze(1, 1) = 0.065807367929181695421;
    tri3_xyze(2, 1) = 0.76147433786105855891;
    nids.push_back(7305);
    tri3_xyze(0, 2) = 0.0091452066580234858684;
    tri3_xyze(1, 2) = 0.069464741084326264708;
    tri3_xyze(2, 2) = 0.75962897079290270241;
    nids.push_back(-959);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.017633031092652306365;
    tri3_xyze(1, 0) = 0.065807367929181695421;
    tri3_xyze(2, 0) = 0.76147433786105855891;
    nids.push_back(7305);
    tri3_xyze(0, 1) = 1.0004171685873872409e-14;
    tri3_xyze(1, 1) = 0.068128800512564396952;
    tri3_xyze(2, 1) = 0.76147433786105844788;
    nids.push_back(7055);
    tri3_xyze(0, 2) = 0.0091452066580234858684;
    tri3_xyze(1, 2) = 0.069464741084326264708;
    tri3_xyze(2, 2) = 0.75962897079290270241;
    nids.push_back(-959);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004482737572315978e-14;
    tri3_xyze(1, 0) = 0.073208660251048829348;
    tri3_xyze(2, 0) = 0.75778360372474695694;
    nids.push_back(7057);
    tri3_xyze(0, 1) = 1.0004819660879340093e-14;
    tri3_xyze(1, 1) = 0.078711035421744998675;
    tri3_xyze(2, 1) = 0.75475864737654663195;
    nids.push_back(7059);
    tri3_xyze(0, 2) = 0.0098299276415846816718;
    tri3_xyze(1, 2) = 0.074665713311280287878;
    tri3_xyze(2, 2) = 0.75627112555064679444;
    nids.push_back(-960);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004819660879340093e-14;
    tri3_xyze(1, 0) = 0.078711035421744998675;
    tri3_xyze(2, 0) = 0.75475864737654663195;
    nids.push_back(7059);
    tri3_xyze(0, 1) = 0.020371915026897093048;
    tri3_xyze(1, 1) = 0.076029021927817172499;
    tri3_xyze(2, 1) = 0.75475864737654663195;
    nids.push_back(7309);
    tri3_xyze(0, 2) = 0.0098299276415846816718;
    tri3_xyze(1, 2) = 0.074665713311280287878;
    tri3_xyze(2, 2) = 0.75627112555064679444;
    nids.push_back(-960);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.020371915026897093048;
    tri3_xyze(1, 0) = 0.076029021927817172499;
    tri3_xyze(2, 0) = 0.75475864737654663195;
    nids.push_back(7309);
    tri3_xyze(0, 1) = 0.018947795539421625338;
    tri3_xyze(1, 1) = 0.070714135644510123235;
    tri3_xyze(2, 1) = 0.75778360372474695694;
    nids.push_back(7307);
    tri3_xyze(0, 2) = 0.0098299276415846816718;
    tri3_xyze(1, 2) = 0.074665713311280287878;
    tri3_xyze(2, 2) = 0.75627112555064679444;
    nids.push_back(-960);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.018947795539421625338;
    tri3_xyze(1, 0) = 0.070714135644510123235;
    tri3_xyze(2, 0) = 0.75778360372474695694;
    nids.push_back(7307);
    tri3_xyze(0, 1) = 1.0004482737572315978e-14;
    tri3_xyze(1, 1) = 0.073208660251048829348;
    tri3_xyze(2, 1) = 0.75778360372474695694;
    nids.push_back(7057);
    tri3_xyze(0, 2) = 0.0098299276415846816718;
    tri3_xyze(1, 2) = 0.074665713311280287878;
    tri3_xyze(2, 2) = 0.75627112555064679444;
    nids.push_back(-960);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004819660879340093e-14;
    tri3_xyze(1, 0) = 0.078711035421744998675;
    tri3_xyze(2, 0) = 0.75475864737654663195;
    nids.push_back(7059);
    tri3_xyze(0, 1) = 1.0005177142313128246e-14;
    tri3_xyze(1, 1) = 0.084549150281251053696;
    tri3_xyze(2, 1) = 0.75244717418508977325;
    nids.push_back(7061);
    tri3_xyze(0, 2) = 0.010563711341737603427;
    tri3_xyze(1, 2) = 0.080239353869567320299;
    tri3_xyze(2, 2) = 0.7536029107808182026;
    nids.push_back(-961);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005177142313128246e-14;
    tri3_xyze(1, 0) = 0.084549150281251053696;
    tri3_xyze(2, 0) = 0.75244717418508977325;
    nids.push_back(7061);
    tri3_xyze(0, 1) = 0.021882930340033308891;
    tri3_xyze(1, 1) = 0.081668207847456056325;
    tri3_xyze(2, 1) = 0.75244717418508977325;
    nids.push_back(7311);
    tri3_xyze(0, 2) = 0.010563711341737603427;
    tri3_xyze(1, 2) = 0.080239353869567320299;
    tri3_xyze(2, 2) = 0.7536029107808182026;
    nids.push_back(-961);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.021882930340033308891;
    tri3_xyze(1, 0) = 0.081668207847456056325;
    tri3_xyze(2, 0) = 0.75244717418508977325;
    nids.push_back(7311);
    tri3_xyze(0, 1) = 0.020371915026897093048;
    tri3_xyze(1, 1) = 0.076029021927817172499;
    tri3_xyze(2, 1) = 0.75475864737654663195;
    nids.push_back(7309);
    tri3_xyze(0, 2) = 0.010563711341737603427;
    tri3_xyze(1, 2) = 0.080239353869567320299;
    tri3_xyze(2, 2) = 0.7536029107808182026;
    nids.push_back(-961);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.020371915026897093048;
    tri3_xyze(1, 0) = 0.076029021927817172499;
    tri3_xyze(2, 0) = 0.75475864737654663195;
    nids.push_back(7309);
    tri3_xyze(0, 1) = 1.0004819660879340093e-14;
    tri3_xyze(1, 1) = 0.078711035421744998675;
    tri3_xyze(2, 1) = 0.75475864737654663195;
    nids.push_back(7059);
    tri3_xyze(0, 2) = 0.010563711341737603427;
    tri3_xyze(1, 2) = 0.080239353869567320299;
    tri3_xyze(2, 2) = 0.7536029107808182026;
    nids.push_back(-961);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005177142313128246e-14;
    tri3_xyze(1, 0) = 0.084549150281251053696;
    tri3_xyze(2, 0) = 0.75244717418508977325;
    nids.push_back(7061);
    tri3_xyze(0, 1) = 1.0005549544177918912e-14;
    tri3_xyze(1, 1) = 0.090630934270712198986;
    tri3_xyze(2, 1) = 0.75088563746341274729;
    nids.push_back(7063);
    tri3_xyze(0, 2) = 0.011334985551189688188;
    tri3_xyze(1, 2) = 0.08609776311805179827;
    tri3_xyze(2, 2) = 0.75166640582425126027;
    nids.push_back(-962);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005549544177918912e-14;
    tri3_xyze(1, 0) = 0.090630934270712198986;
    tri3_xyze(2, 0) = 0.75088563746341274729;
    nids.push_back(7063);
    tri3_xyze(0, 1) = 0.023457011864705432092;
    tri3_xyze(1, 1) = 0.087542760072787856318;
    tri3_xyze(2, 1) = 0.75088563746341285832;
    nids.push_back(7313);
    tri3_xyze(0, 2) = 0.011334985551189688188;
    tri3_xyze(1, 2) = 0.08609776311805179827;
    tri3_xyze(2, 2) = 0.75166640582425126027;
    nids.push_back(-962);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.023457011864705432092;
    tri3_xyze(1, 0) = 0.087542760072787856318;
    tri3_xyze(2, 0) = 0.75088563746341285832;
    nids.push_back(7313);
    tri3_xyze(0, 1) = 0.021882930340033308891;
    tri3_xyze(1, 1) = 0.081668207847456056325;
    tri3_xyze(2, 1) = 0.75244717418508977325;
    nids.push_back(7311);
    tri3_xyze(0, 2) = 0.011334985551189688188;
    tri3_xyze(1, 2) = 0.08609776311805179827;
    tri3_xyze(2, 2) = 0.75166640582425126027;
    nids.push_back(-962);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.021882930340033308891;
    tri3_xyze(1, 0) = 0.081668207847456056325;
    tri3_xyze(2, 0) = 0.75244717418508977325;
    nids.push_back(7311);
    tri3_xyze(0, 1) = 1.0005177142313128246e-14;
    tri3_xyze(1, 1) = 0.084549150281251053696;
    tri3_xyze(2, 1) = 0.75244717418508977325;
    nids.push_back(7061);
    tri3_xyze(0, 2) = 0.011334985551189688188;
    tri3_xyze(1, 2) = 0.08609776311805179827;
    tri3_xyze(2, 2) = 0.75166640582425126027;
    nids.push_back(-962);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005549544177918912e-14;
    tri3_xyze(1, 0) = 0.090630934270712198986;
    tri3_xyze(2, 0) = 0.75088563746341274729;
    nids.push_back(7063);
    tri3_xyze(0, 1) = 1.0005930993473840787e-14;
    tri3_xyze(1, 1) = 0.096860474023532747156;
    tri3_xyze(2, 1) = 0.75009866357843346663;
    nids.push_back(7065);
    tri3_xyze(0, 2) = 0.012131586814920991557;
    tri3_xyze(1, 2) = 0.092148550443241078467;
    tri3_xyze(2, 2) = 0.75049215052092321798;
    nids.push_back(-963);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005930993473840787e-14;
    tri3_xyze(1, 0) = 0.096860474023532747156;
    tri3_xyze(2, 0) = 0.75009866357843346663;
    nids.push_back(7065);
    tri3_xyze(0, 1) = 0.025069335394958522367;
    tri3_xyze(1, 1) = 0.093560033405931525285;
    tri3_xyze(2, 1) = 0.75009866357843346663;
    nids.push_back(7315);
    tri3_xyze(0, 2) = 0.012131586814920991557;
    tri3_xyze(1, 2) = 0.092148550443241078467;
    tri3_xyze(2, 2) = 0.75049215052092321798;
    nids.push_back(-963);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025069335394958522367;
    tri3_xyze(1, 0) = 0.093560033405931525285;
    tri3_xyze(2, 0) = 0.75009866357843346663;
    nids.push_back(7315);
    tri3_xyze(0, 1) = 0.023457011864705432092;
    tri3_xyze(1, 1) = 0.087542760072787856318;
    tri3_xyze(2, 1) = 0.75088563746341285832;
    nids.push_back(7313);
    tri3_xyze(0, 2) = 0.012131586814920991557;
    tri3_xyze(1, 2) = 0.092148550443241078467;
    tri3_xyze(2, 2) = 0.75049215052092321798;
    nids.push_back(-963);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.023457011864705432092;
    tri3_xyze(1, 0) = 0.087542760072787856318;
    tri3_xyze(2, 0) = 0.75088563746341285832;
    nids.push_back(7313);
    tri3_xyze(0, 1) = 1.0005549544177918912e-14;
    tri3_xyze(1, 1) = 0.090630934270712198986;
    tri3_xyze(2, 1) = 0.75088563746341274729;
    nids.push_back(7063);
    tri3_xyze(0, 2) = 0.012131586814920991557;
    tri3_xyze(1, 2) = 0.092148550443241078467;
    tri3_xyze(2, 2) = 0.75049215052092321798;
    nids.push_back(-963);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005930993473840787e-14;
    tri3_xyze(1, 0) = 0.096860474023532747156;
    tri3_xyze(2, 0) = 0.75009866357843346663;
    nids.push_back(7065);
    tri3_xyze(0, 1) = 1.0006315474517633335e-14;
    tri3_xyze(1, 1) = 0.10313952597646400267;
    tri3_xyze(2, 1) = 0.75009866357843324458;
    nids.push_back(7067);
    tri3_xyze(0, 2) = 0.012940952255136001225;
    tri3_xyze(1, 2) = 0.098296291314451805787;
    tri3_xyze(2, 2) = 0.7500986635784333556;
    nids.push_back(-964);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0006315474517633335e-14;
    tri3_xyze(1, 0) = 0.10313952597646400267;
    tri3_xyze(2, 0) = 0.75009866357843324458;
    nids.push_back(7067);
    tri3_xyze(0, 1) = 0.026694473625565474234;
    tri3_xyze(1, 1) = 0.099625131851878920286;
    tri3_xyze(2, 1) = 0.7500986635784333556;
    nids.push_back(7317);
    tri3_xyze(0, 2) = 0.012940952255136001225;
    tri3_xyze(1, 2) = 0.098296291314451805787;
    tri3_xyze(2, 2) = 0.7500986635784333556;
    nids.push_back(-964);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026694473625565474234;
    tri3_xyze(1, 0) = 0.099625131851878920286;
    tri3_xyze(2, 0) = 0.7500986635784333556;
    nids.push_back(7317);
    tri3_xyze(0, 1) = 0.025069335394958522367;
    tri3_xyze(1, 1) = 0.093560033405931525285;
    tri3_xyze(2, 1) = 0.75009866357843346663;
    nids.push_back(7315);
    tri3_xyze(0, 2) = 0.012940952255136001225;
    tri3_xyze(1, 2) = 0.098296291314451805787;
    tri3_xyze(2, 2) = 0.7500986635784333556;
    nids.push_back(-964);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025069335394958522367;
    tri3_xyze(1, 0) = 0.093560033405931525285;
    tri3_xyze(2, 0) = 0.75009866357843346663;
    nids.push_back(7315);
    tri3_xyze(0, 1) = 1.0005930993473840787e-14;
    tri3_xyze(1, 1) = 0.096860474023532747156;
    tri3_xyze(2, 1) = 0.75009866357843346663;
    nids.push_back(7065);
    tri3_xyze(0, 2) = 0.012940952255136001225;
    tri3_xyze(1, 2) = 0.098296291314451805787;
    tri3_xyze(2, 2) = 0.7500986635784333556;
    nids.push_back(-964);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.028306797155818672063;
    tri3_xyze(1, 0) = 0.10564240518502271415;
    tri3_xyze(2, 0) = 0.75088563746341230321;
    nids.push_back(7319);
    tri3_xyze(0, 1) = 0.026694473625565474234;
    tri3_xyze(1, 1) = 0.099625131851878920286;
    tri3_xyze(2, 1) = 0.7500986635784333556;
    nids.push_back(7317);
    tri3_xyze(0, 2) = 0.013750317695351038649;
    tri3_xyze(1, 2) = 0.10444403218566257474;
    tri3_xyze(2, 2) = 0.75049215052092277389;
    nids.push_back(-965);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026694473625565474234;
    tri3_xyze(1, 0) = 0.099625131851878920286;
    tri3_xyze(2, 0) = 0.7500986635784333556;
    nids.push_back(7317);
    tri3_xyze(0, 1) = 1.0006315474517633335e-14;
    tri3_xyze(1, 1) = 0.10313952597646400267;
    tri3_xyze(2, 1) = 0.75009866357843324458;
    nids.push_back(7067);
    tri3_xyze(0, 2) = 0.013750317695351038649;
    tri3_xyze(1, 2) = 0.10444403218566257474;
    tri3_xyze(2, 2) = 0.75049215052092277389;
    nids.push_back(-965);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.8759580880385857084e-15;
    tri3_xyze(1, 0) = 0.1500000000000000222;
    tri3_xyze(2, 0) = 0.79999999999984527932;
    nids.push_back(7091);
    tri3_xyze(0, 1) = 9.7434075385316158324e-15;
    tri3_xyze(1, 1) = 0.14960573506572408253;
    tri3_xyze(2, 1) = 0.80626666167806060326;
    nids.push_back(7093);
    tri3_xyze(0, 2) = 0.019385917564247255596;
    tri3_xyze(1, 2) = 0.14725066306750683731;
    tri3_xyze(2, 2) = 0.8031333308389529968;
    nids.push_back(-977);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.7434075385316158324e-15;
    tri3_xyze(1, 0) = 0.14960573506572408253;
    tri3_xyze(2, 0) = 0.80626666167806060326;
    nids.push_back(7093);
    tri3_xyze(0, 1) = 0.038720813491581330956;
    tri3_xyze(1, 1) = 0.14450804326094299213;
    tri3_xyze(2, 1) = 0.80626666167806060326;
    nids.push_back(7343);
    tri3_xyze(0, 2) = 0.019385917564247255596;
    tri3_xyze(1, 2) = 0.14725066306750683731;
    tri3_xyze(2, 2) = 0.8031333308389529968;
    nids.push_back(-977);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.038720813491581330956;
    tri3_xyze(1, 0) = 0.14450804326094299213;
    tri3_xyze(2, 0) = 0.80626666167806060326;
    nids.push_back(7343);
    tri3_xyze(0, 1) = 0.038822856765388082112;
    tri3_xyze(1, 1) = 0.14488887394336025238;
    tri3_xyze(2, 1) = 0.79999999999984550136;
    nids.push_back(7341);
    tri3_xyze(0, 2) = 0.019385917564247255596;
    tri3_xyze(1, 2) = 0.14725066306750683731;
    tri3_xyze(2, 2) = 0.8031333308389529968;
    nids.push_back(-977);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.038822856765388082112;
    tri3_xyze(1, 0) = 0.14488887394336025238;
    tri3_xyze(2, 0) = 0.79999999999984550136;
    nids.push_back(7341);
    tri3_xyze(0, 1) = 9.8759580880385857084e-15;
    tri3_xyze(1, 1) = 0.1500000000000000222;
    tri3_xyze(2, 1) = 0.79999999999984527932;
    nids.push_back(7091);
    tri3_xyze(0, 2) = 0.019385917564247255596;
    tri3_xyze(1, 2) = 0.14725066306750683731;
    tri3_xyze(2, 2) = 0.8031333308389529968;
    nids.push_back(-977);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.7434075385316158324e-15;
    tri3_xyze(1, 0) = 0.14960573506572408253;
    tri3_xyze(2, 0) = 0.80626666167806060326;
    nids.push_back(7093);
    tri3_xyze(0, 1) = 9.7454255145857347186e-15;
    tri3_xyze(1, 1) = 0.14842915805643194793;
    tri3_xyze(2, 1) = 0.81243449435808812176;
    nids.push_back(7095);
    tri3_xyze(0, 2) = 0.019284276611286881725;
    tri3_xyze(1, 2) = 0.14647862338103720403;
    tri3_xyze(2, 2) = 0.80935057801807441802;
    nids.push_back(-978);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.7454255145857347186e-15;
    tri3_xyze(1, 0) = 0.14842915805643194793;
    tri3_xyze(2, 0) = 0.81243449435808812176;
    nids.push_back(7095);
    tri3_xyze(0, 1) = 0.038416292953546704592;
    tri3_xyze(1, 1) = 0.14337155714104973803;
    tri3_xyze(2, 1) = 0.81243449435808823278;
    nids.push_back(7345);
    tri3_xyze(0, 2) = 0.019284276611286881725;
    tri3_xyze(1, 2) = 0.14647862338103720403;
    tri3_xyze(2, 2) = 0.80935057801807441802;
    nids.push_back(-978);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.038416292953546704592;
    tri3_xyze(1, 0) = 0.14337155714104973803;
    tri3_xyze(2, 0) = 0.81243449435808823278;
    nids.push_back(7345);
    tri3_xyze(0, 1) = 0.038720813491581330956;
    tri3_xyze(1, 1) = 0.14450804326094299213;
    tri3_xyze(2, 1) = 0.80626666167806060326;
    nids.push_back(7343);
    tri3_xyze(0, 2) = 0.019284276611286881725;
    tri3_xyze(1, 2) = 0.14647862338103720403;
    tri3_xyze(2, 2) = 0.80935057801807441802;
    nids.push_back(-978);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.038720813491581330956;
    tri3_xyze(1, 0) = 0.14450804326094299213;
    tri3_xyze(2, 0) = 0.80626666167806060326;
    nids.push_back(7343);
    tri3_xyze(0, 1) = 9.7434075385316158324e-15;
    tri3_xyze(1, 1) = 0.14960573506572408253;
    tri3_xyze(2, 1) = 0.80626666167806060326;
    nids.push_back(7093);
    tri3_xyze(0, 2) = 0.019284276611286881725;
    tri3_xyze(1, 2) = 0.14647862338103720403;
    tri3_xyze(2, 2) = 0.80935057801807441802;
    nids.push_back(-978);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.7454255145857347186e-15;
    tri3_xyze(1, 0) = 0.14842915805643194793;
    tri3_xyze(2, 0) = 0.81243449435808812176;
    nids.push_back(7095);
    tri3_xyze(0, 1) = 1.0008969853489149771e-14;
    tri3_xyze(1, 1) = 0.14648882429441315023;
    tri3_xyze(2, 1) = 0.81840622763407933071;
    nids.push_back(7097);
    tri3_xyze(0, 2) = 0.019082597643911805574;
    tri3_xyze(1, 2) = 0.14494671953514751128;
    tri3_xyze(2, 2) = 0.81542036099608383726;
    nids.push_back(-979);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008969853489149771e-14;
    tri3_xyze(1, 0) = 0.14648882429441315023;
    tri3_xyze(2, 0) = 0.81840622763407933071;
    nids.push_back(7097);
    tri3_xyze(0, 1) = 0.037914097622080769612;
    tri3_xyze(1, 1) = 0.14149733864869518118;
    tri3_xyze(2, 1) = 0.81840622763407955276;
    nids.push_back(7347);
    tri3_xyze(0, 2) = 0.019082597643911805574;
    tri3_xyze(1, 2) = 0.14494671953514751128;
    tri3_xyze(2, 2) = 0.81542036099608383726;
    nids.push_back(-979);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.037914097622080769612;
    tri3_xyze(1, 0) = 0.14149733864869518118;
    tri3_xyze(2, 0) = 0.81840622763407955276;
    nids.push_back(7347);
    tri3_xyze(0, 1) = 0.038416292953546704592;
    tri3_xyze(1, 1) = 0.14337155714104973803;
    tri3_xyze(2, 1) = 0.81243449435808823278;
    nids.push_back(7345);
    tri3_xyze(0, 2) = 0.019082597643911805574;
    tri3_xyze(1, 2) = 0.14494671953514751128;
    tri3_xyze(2, 2) = 0.81542036099608383726;
    nids.push_back(-979);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.038416292953546704592;
    tri3_xyze(1, 0) = 0.14337155714104973803;
    tri3_xyze(2, 0) = 0.81243449435808823278;
    nids.push_back(7345);
    tri3_xyze(0, 1) = 9.7454255145857347186e-15;
    tri3_xyze(1, 1) = 0.14842915805643194793;
    tri3_xyze(2, 1) = 0.81243449435808812176;
    nids.push_back(7095);
    tri3_xyze(0, 2) = 0.019082597643911805574;
    tri3_xyze(1, 2) = 0.14494671953514751128;
    tri3_xyze(2, 2) = 0.81542036099608383726;
    nids.push_back(-979);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008969853489149771e-14;
    tri3_xyze(1, 0) = 0.14648882429441315023;
    tri3_xyze(2, 0) = 0.81840622763407933071;
    nids.push_back(7097);
    tri3_xyze(0, 1) = 1.0008806149422705407e-14;
    tri3_xyze(1, 1) = 0.14381533400219395635;
    tri3_xyze(2, 1) = 0.82408768370493123978;
    nids.push_back(7099);
    tri3_xyze(0, 2) = 0.018784061259914618808;
    tri3_xyze(1, 2) = 0.14267911056860244079;
    tri3_xyze(2, 2) = 0.82124695566950534076;
    nids.push_back(-980);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008806149422705407e-14;
    tri3_xyze(1, 0) = 0.14381533400219395635;
    tri3_xyze(2, 0) = 0.82408768370493123978;
    nids.push_back(7099);
    tri3_xyze(0, 1) = 0.03722214741755768691;
    tri3_xyze(1, 1) = 0.1389149453291075309;
    tri3_xyze(2, 1) = 0.82408768370493157285;
    nids.push_back(7349);
    tri3_xyze(0, 2) = 0.018784061259914618808;
    tri3_xyze(1, 2) = 0.14267911056860244079;
    tri3_xyze(2, 2) = 0.82124695566950534076;
    nids.push_back(-980);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.03722214741755768691;
    tri3_xyze(1, 0) = 0.1389149453291075309;
    tri3_xyze(2, 0) = 0.82408768370493157285;
    nids.push_back(7349);
    tri3_xyze(0, 1) = 0.037914097622080769612;
    tri3_xyze(1, 1) = 0.14149733864869518118;
    tri3_xyze(2, 1) = 0.81840622763407955276;
    nids.push_back(7347);
    tri3_xyze(0, 2) = 0.018784061259914618808;
    tri3_xyze(1, 2) = 0.14267911056860244079;
    tri3_xyze(2, 2) = 0.82124695566950534076;
    nids.push_back(-980);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.037914097622080769612;
    tri3_xyze(1, 0) = 0.14149733864869518118;
    tri3_xyze(2, 0) = 0.81840622763407955276;
    nids.push_back(7347);
    tri3_xyze(0, 1) = 1.0008969853489149771e-14;
    tri3_xyze(1, 1) = 0.14648882429441315023;
    tri3_xyze(2, 1) = 0.81840622763407933071;
    nids.push_back(7097);
    tri3_xyze(0, 2) = 0.018784061259914618808;
    tri3_xyze(1, 2) = 0.14267911056860244079;
    tri3_xyze(2, 2) = 0.82124695566950534076;
    nids.push_back(-980);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008806149422705407e-14;
    tri3_xyze(1, 0) = 0.14381533400219395635;
    tri3_xyze(2, 0) = 0.82408768370493123978;
    nids.push_back(7099);
    tri3_xyze(0, 1) = 1.000860013417727975e-14;
    tri3_xyze(1, 1) = 0.14045084971874827184;
    tri3_xyze(2, 1) = 0.82938926261446932475;
    nids.push_back(7101);
    tri3_xyze(0, 2) = 0.018393375556407874638;
    tri3_xyze(1, 2) = 0.13971155802940837809;
    tri3_xyze(2, 2) = 0.82673847315970050431;
    nids.push_back(-981);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.000860013417727975e-14;
    tri3_xyze(1, 0) = 0.14045084971874827184;
    tri3_xyze(2, 0) = 0.82938926261446932475;
    nids.push_back(7101);
    tri3_xyze(0, 1) = 0.036351354808053792933;
    tri3_xyze(1, 1) = 0.13566510306758366999;
    tri3_xyze(2, 1) = 0.8293892626144695468;
    nids.push_back(7351);
    tri3_xyze(0, 2) = 0.018393375556407874638;
    tri3_xyze(1, 2) = 0.13971155802940837809;
    tri3_xyze(2, 2) = 0.82673847315970050431;
    nids.push_back(-981);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036351354808053792933;
    tri3_xyze(1, 0) = 0.13566510306758366999;
    tri3_xyze(2, 0) = 0.8293892626144695468;
    nids.push_back(7351);
    tri3_xyze(0, 1) = 0.03722214741755768691;
    tri3_xyze(1, 1) = 0.1389149453291075309;
    tri3_xyze(2, 1) = 0.82408768370493157285;
    nids.push_back(7349);
    tri3_xyze(0, 2) = 0.018393375556407874638;
    tri3_xyze(1, 2) = 0.13971155802940837809;
    tri3_xyze(2, 2) = 0.82673847315970050431;
    nids.push_back(-981);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.03722214741755768691;
    tri3_xyze(1, 0) = 0.1389149453291075309;
    tri3_xyze(2, 0) = 0.82408768370493157285;
    nids.push_back(7349);
    tri3_xyze(0, 1) = 1.0008806149422705407e-14;
    tri3_xyze(1, 1) = 0.14381533400219395635;
    tri3_xyze(2, 1) = 0.82408768370493123978;
    nids.push_back(7099);
    tri3_xyze(0, 2) = 0.018393375556407874638;
    tri3_xyze(1, 2) = 0.13971155802940837809;
    tri3_xyze(2, 2) = 0.82673847315970050431;
    nids.push_back(-981);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.000860013417727975e-14;
    tri3_xyze(1, 0) = 0.14045084971874827184;
    tri3_xyze(2, 0) = 0.82938926261446932475;
    nids.push_back(7101);
    tri3_xyze(0, 1) = 1.0008355056736363408e-14;
    tri3_xyze(1, 1) = 0.13644843137107168229;
    tri3_xyze(2, 1) = 0.8342273552962801153;
    nids.push_back(7103);
    tri3_xyze(0, 2) = 0.017916701880320275164;
    tri3_xyze(1, 2) = 0.13609086199383835236;
    tri3_xyze(2, 2) = 0.83180830895537483105;
    nids.push_back(-982);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008355056736363408e-14;
    tri3_xyze(1, 0) = 0.13644843137107168229;
    tri3_xyze(2, 0) = 0.8342273552962801153;
    nids.push_back(7103);
    tri3_xyze(0, 1) = 0.035315452713207295954;
    tri3_xyze(1, 1) = 0.13179906381794967429;
    tri3_xyze(2, 1) = 0.83422735529628033735;
    nids.push_back(7353);
    tri3_xyze(0, 2) = 0.017916701880320275164;
    tri3_xyze(1, 2) = 0.13609086199383835236;
    tri3_xyze(2, 2) = 0.83180830895537483105;
    nids.push_back(-982);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.035315452713207295954;
    tri3_xyze(1, 0) = 0.13179906381794967429;
    tri3_xyze(2, 0) = 0.83422735529628033735;
    nids.push_back(7353);
    tri3_xyze(0, 1) = 0.036351354808053792933;
    tri3_xyze(1, 1) = 0.13566510306758366999;
    tri3_xyze(2, 1) = 0.8293892626144695468;
    nids.push_back(7351);
    tri3_xyze(0, 2) = 0.017916701880320275164;
    tri3_xyze(1, 2) = 0.13609086199383835236;
    tri3_xyze(2, 2) = 0.83180830895537483105;
    nids.push_back(-982);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036351354808053792933;
    tri3_xyze(1, 0) = 0.13566510306758366999;
    tri3_xyze(2, 0) = 0.8293892626144695468;
    nids.push_back(7351);
    tri3_xyze(0, 1) = 1.000860013417727975e-14;
    tri3_xyze(1, 1) = 0.14045084971874827184;
    tri3_xyze(2, 1) = 0.82938926261446932475;
    nids.push_back(7101);
    tri3_xyze(0, 2) = 0.017916701880320275164;
    tri3_xyze(1, 2) = 0.13609086199383835236;
    tri3_xyze(2, 2) = 0.83180830895537483105;
    nids.push_back(-982);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008355056736363408e-14;
    tri3_xyze(1, 0) = 0.13644843137107168229;
    tri3_xyze(2, 0) = 0.8342273552962801153;
    nids.push_back(7103);
    tri3_xyze(0, 1) = 1.0008074782117600136e-14;
    tri3_xyze(1, 1) = 0.13187119948743578068;
    tri3_xyze(2, 1) = 0.8385256621386353304;
    nids.push_back(7105);
    tri3_xyze(0, 2) = 0.017361557660274772041;
    tri3_xyze(1, 2) = 0.13187412300127224207;
    tri3_xyze(2, 2) = 0.83637650871745783387;
    nids.push_back(-983);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008074782117600136e-14;
    tri3_xyze(1, 0) = 0.13187119948743578068;
    tri3_xyze(2, 0) = 0.8385256621386353304;
    nids.push_back(7105);
    tri3_xyze(0, 1) = 0.034130777927871766564;
    tri3_xyze(1, 1) = 0.12737779732863185878;
    tri3_xyze(2, 1) = 0.83852566213863555245;
    nids.push_back(7355);
    tri3_xyze(0, 2) = 0.017361557660274772041;
    tri3_xyze(1, 2) = 0.13187412300127224207;
    tri3_xyze(2, 2) = 0.83637650871745783387;
    nids.push_back(-983);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.034130777927871766564;
    tri3_xyze(1, 0) = 0.12737779732863185878;
    tri3_xyze(2, 0) = 0.83852566213863555245;
    nids.push_back(7355);
    tri3_xyze(0, 1) = 0.035315452713207295954;
    tri3_xyze(1, 1) = 0.13179906381794967429;
    tri3_xyze(2, 1) = 0.83422735529628033735;
    nids.push_back(7353);
    tri3_xyze(0, 2) = 0.017361557660274772041;
    tri3_xyze(1, 2) = 0.13187412300127224207;
    tri3_xyze(2, 2) = 0.83637650871745783387;
    nids.push_back(-983);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.035315452713207295954;
    tri3_xyze(1, 0) = 0.13179906381794967429;
    tri3_xyze(2, 0) = 0.83422735529628033735;
    nids.push_back(7353);
    tri3_xyze(0, 1) = 1.0008355056736363408e-14;
    tri3_xyze(1, 1) = 0.13644843137107168229;
    tri3_xyze(2, 1) = 0.8342273552962801153;
    nids.push_back(7103);
    tri3_xyze(0, 2) = 0.017361557660274772041;
    tri3_xyze(1, 2) = 0.13187412300127224207;
    tri3_xyze(2, 2) = 0.83637650871745783387;
    nids.push_back(-983);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0008074782117600136e-14;
    tri3_xyze(1, 0) = 0.13187119948743578068;
    tri3_xyze(2, 0) = 0.8385256621386353304;
    nids.push_back(7105);
    tri3_xyze(0, 1) = 1.0007763730419158145e-14;
    tri3_xyze(1, 1) = 0.12679133974895126502;
    tri3_xyze(2, 1) = 0.84221639627494671032;
    nids.push_back(7107);
    tri3_xyze(0, 2) = 0.016736697852248565155;
    tri3_xyze(1, 2) = 0.12712784154458067754;
    tri3_xyze(2, 2) = 0.84037102920679107587;
    nids.push_back(-984);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0007763730419158145e-14;
    tri3_xyze(1, 0) = 0.12679133974895126502;
    tri3_xyze(2, 0) = 0.84221639627494671032;
    nids.push_back(7107);
    tri3_xyze(0, 1) = 0.032816013481102468408;
    tri3_xyze(1, 1) = 0.12247102961330375015;
    tri3_xyze(2, 1) = 0.84221639627494693237;
    nids.push_back(7357);
    tri3_xyze(0, 2) = 0.016736697852248565155;
    tri3_xyze(1, 2) = 0.12712784154458067754;
    tri3_xyze(2, 2) = 0.84037102920679107587;
    nids.push_back(-984);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.032816013481102468408;
    tri3_xyze(1, 0) = 0.12247102961330375015;
    tri3_xyze(2, 0) = 0.84221639627494693237;
    nids.push_back(7357);
    tri3_xyze(0, 1) = 0.034130777927871766564;
    tri3_xyze(1, 1) = 0.12737779732863185878;
    tri3_xyze(2, 1) = 0.83852566213863555245;
    nids.push_back(7355);
    tri3_xyze(0, 2) = 0.016736697852248565155;
    tri3_xyze(1, 2) = 0.12712784154458067754;
    tri3_xyze(2, 2) = 0.84037102920679107587;
    nids.push_back(-984);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.034130777927871766564;
    tri3_xyze(1, 0) = 0.12737779732863185878;
    tri3_xyze(2, 0) = 0.83852566213863555245;
    nids.push_back(7355);
    tri3_xyze(0, 1) = 1.0008074782117600136e-14;
    tri3_xyze(1, 1) = 0.13187119948743578068;
    tri3_xyze(2, 1) = 0.8385256621386353304;
    nids.push_back(7105);
    tri3_xyze(0, 2) = 0.016736697852248565155;
    tri3_xyze(1, 2) = 0.12712784154458067754;
    tri3_xyze(2, 2) = 0.84037102920679107587;
    nids.push_back(-984);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0007763730419158145e-14;
    tri3_xyze(1, 0) = 0.12679133974895126502;
    tri3_xyze(2, 0) = 0.84221639627494671032;
    nids.push_back(7107);
    tri3_xyze(0, 1) = 1.0222879288884487869e-14;
    tri3_xyze(1, 1) = 0.12128896457825505406;
    tri3_xyze(2, 1) = 0.84524135262314725736;
    nids.push_back(7109);
    tri3_xyze(0, 2) = 0.016051976868687437006;
    tri3_xyze(1, 2) = 0.12192686931762668212;
    tri3_xyze(2, 2) = 0.84372887444904698384;
    nids.push_back(-985);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0222879288884487869e-14;
    tri3_xyze(1, 0) = 0.12128896457825505406;
    tri3_xyze(2, 0) = 0.84524135262314725736;
    nids.push_back(7109);
    tri3_xyze(0, 1) = 0.031391893993627052739;
    tri3_xyze(1, 1) = 0.11715614332999664537;
    tri3_xyze(2, 1) = 0.84524135262314736838;
    nids.push_back(7359);
    tri3_xyze(0, 2) = 0.016051976868687437006;
    tri3_xyze(1, 2) = 0.12192686931762668212;
    tri3_xyze(2, 2) = 0.84372887444904698384;
    nids.push_back(-985);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.031391893993627052739;
    tri3_xyze(1, 0) = 0.11715614332999664537;
    tri3_xyze(2, 0) = 0.84524135262314736838;
    nids.push_back(7359);
    tri3_xyze(0, 1) = 0.032816013481102468408;
    tri3_xyze(1, 1) = 0.12247102961330375015;
    tri3_xyze(2, 1) = 0.84221639627494693237;
    nids.push_back(7357);
    tri3_xyze(0, 2) = 0.016051976868687437006;
    tri3_xyze(1, 2) = 0.12192686931762668212;
    tri3_xyze(2, 2) = 0.84372887444904698384;
    nids.push_back(-985);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.032816013481102468408;
    tri3_xyze(1, 0) = 0.12247102961330375015;
    tri3_xyze(2, 0) = 0.84221639627494693237;
    nids.push_back(7357);
    tri3_xyze(0, 1) = 1.0007763730419158145e-14;
    tri3_xyze(1, 1) = 0.12679133974895126502;
    tri3_xyze(2, 1) = 0.84221639627494671032;
    nids.push_back(7107);
    tri3_xyze(0, 2) = 0.016051976868687437006;
    tri3_xyze(1, 2) = 0.12192686931762668212;
    tri3_xyze(2, 2) = 0.84372887444904698384;
    nids.push_back(-985);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0222879288884487869e-14;
    tri3_xyze(1, 0) = 0.12128896457825505406;
    tri3_xyze(2, 0) = 0.84524135262314725736;
    nids.push_back(7109);
    tri3_xyze(0, 1) = 1.0109610278934578421e-14;
    tri3_xyze(1, 1) = 0.11545084971874891577;
    tri3_xyze(2, 1) = 0.84755282581460400504;
    nids.push_back(7111);
    tri3_xyze(0, 2) = 0.01531819316853454821;
    tri3_xyze(1, 2) = 0.11635322875933955256;
    tri3_xyze(2, 2) = 0.84639708921887568671;
    nids.push_back(-986);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0109610278934578421e-14;
    tri3_xyze(1, 0) = 0.11545084971874891577;
    tri3_xyze(2, 0) = 0.84755282581460400504;
    nids.push_back(7111);
    tri3_xyze(0, 1) = 0.029880878680490805671;
    tri3_xyze(1, 1) = 0.11151695741035763665;
    tri3_xyze(2, 1) = 0.84755282581460411606;
    nids.push_back(7361);
    tri3_xyze(0, 2) = 0.01531819316853454821;
    tri3_xyze(1, 2) = 0.11635322875933955256;
    tri3_xyze(2, 2) = 0.84639708921887568671;
    nids.push_back(-986);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.029880878680490805671;
    tri3_xyze(1, 0) = 0.11151695741035763665;
    tri3_xyze(2, 0) = 0.84755282581460411606;
    nids.push_back(7361);
    tri3_xyze(0, 1) = 0.031391893993627052739;
    tri3_xyze(1, 1) = 0.11715614332999664537;
    tri3_xyze(2, 1) = 0.84524135262314736838;
    nids.push_back(7359);
    tri3_xyze(0, 2) = 0.01531819316853454821;
    tri3_xyze(1, 2) = 0.11635322875933955256;
    tri3_xyze(2, 2) = 0.84639708921887568671;
    nids.push_back(-986);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.031391893993627052739;
    tri3_xyze(1, 0) = 0.11715614332999664537;
    tri3_xyze(2, 0) = 0.84524135262314736838;
    nids.push_back(7359);
    tri3_xyze(0, 1) = 1.0222879288884487869e-14;
    tri3_xyze(1, 1) = 0.12128896457825505406;
    tri3_xyze(2, 1) = 0.84524135262314725736;
    nids.push_back(7109);
    tri3_xyze(0, 2) = 0.01531819316853454821;
    tri3_xyze(1, 2) = 0.11635322875933955256;
    tri3_xyze(2, 2) = 0.84639708921887568671;
    nids.push_back(-986);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0109610278934578421e-14;
    tri3_xyze(1, 0) = 0.11545084971874891577;
    tri3_xyze(2, 0) = 0.84755282581460400504;
    nids.push_back(7111);
    tri3_xyze(0, 1) = 9.8124184358865223691e-15;
    tri3_xyze(1, 1) = 0.10936906572928781212;
    tri3_xyze(2, 1) = 0.84911436253628103099;
    nids.push_back(7113);
    tri3_xyze(0, 2) = 0.014546918959082354161;
    tri3_xyze(1, 2) = 0.11049481951085504683;
    tri3_xyze(2, 2) = 0.84833359417544251802;
    nids.push_back(-987);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.8124184358865223691e-15;
    tri3_xyze(1, 0) = 0.10936906572928781212;
    tri3_xyze(2, 0) = 0.84911436253628103099;
    nids.push_back(7113);
    tri3_xyze(0, 1) = 0.02830679715581868594;
    tri3_xyze(1, 1) = 0.10564240518502583666;
    tri3_xyze(2, 1) = 0.84911436253628114201;
    nids.push_back(7363);
    tri3_xyze(0, 2) = 0.014546918959082354161;
    tri3_xyze(1, 2) = 0.11049481951085504683;
    tri3_xyze(2, 2) = 0.84833359417544251802;
    nids.push_back(-987);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.02830679715581868594;
    tri3_xyze(1, 0) = 0.10564240518502583666;
    tri3_xyze(2, 0) = 0.84911436253628114201;
    nids.push_back(7363);
    tri3_xyze(0, 1) = 0.029880878680490805671;
    tri3_xyze(1, 1) = 0.11151695741035763665;
    tri3_xyze(2, 1) = 0.84755282581460411606;
    nids.push_back(7361);
    tri3_xyze(0, 2) = 0.014546918959082354161;
    tri3_xyze(1, 2) = 0.11049481951085504683;
    tri3_xyze(2, 2) = 0.84833359417544251802;
    nids.push_back(-987);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.029880878680490805671;
    tri3_xyze(1, 0) = 0.11151695741035763665;
    tri3_xyze(2, 0) = 0.84755282581460411606;
    nids.push_back(7361);
    tri3_xyze(0, 1) = 1.0109610278934578421e-14;
    tri3_xyze(1, 1) = 0.11545084971874891577;
    tri3_xyze(2, 1) = 0.84755282581460400504;
    nids.push_back(7111);
    tri3_xyze(0, 2) = 0.014546918959082354161;
    tri3_xyze(1, 2) = 0.11049481951085504683;
    tri3_xyze(2, 2) = 0.84833359417544251802;
    nids.push_back(-987);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.8124184358865223691e-15;
    tri3_xyze(1, 0) = 0.10936906572928781212;
    tri3_xyze(2, 0) = 0.84911436253628103099;
    nids.push_back(7113);
    tri3_xyze(0, 1) = 9.9147091733272353276e-15;
    tri3_xyze(1, 1) = 0.10313952597646713905;
    tri3_xyze(2, 1) = 0.84990133642126042268;
    nids.push_back(7115);
    tri3_xyze(0, 2) = 0.01375031769535093977;
    tri3_xyze(1, 2) = 0.10444403218566568337;
    tri3_xyze(2, 2) = 0.84950784947877078235;
    nids.push_back(-988);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.9147091733272353276e-15;
    tri3_xyze(1, 0) = 0.10313952597646713905;
    tri3_xyze(2, 0) = 0.84990133642126042268;
    nids.push_back(7115);
    tri3_xyze(0, 1) = 0.026694473625565345865;
    tri3_xyze(1, 1) = 0.099625131851881945644;
    tri3_xyze(2, 1) = 0.8499013364212605337;
    nids.push_back(7365);
    tri3_xyze(0, 2) = 0.01375031769535093977;
    tri3_xyze(1, 2) = 0.10444403218566568337;
    tri3_xyze(2, 2) = 0.84950784947877078235;
    nids.push_back(-988);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026694473625565345865;
    tri3_xyze(1, 0) = 0.099625131851881945644;
    tri3_xyze(2, 0) = 0.8499013364212605337;
    nids.push_back(7365);
    tri3_xyze(0, 1) = 0.02830679715581868594;
    tri3_xyze(1, 1) = 0.10564240518502583666;
    tri3_xyze(2, 1) = 0.84911436253628114201;
    nids.push_back(7363);
    tri3_xyze(0, 2) = 0.01375031769535093977;
    tri3_xyze(1, 2) = 0.10444403218566568337;
    tri3_xyze(2, 2) = 0.84950784947877078235;
    nids.push_back(-988);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.02830679715581868594;
    tri3_xyze(1, 0) = 0.10564240518502583666;
    tri3_xyze(2, 0) = 0.84911436253628114201;
    nids.push_back(7363);
    tri3_xyze(0, 1) = 9.8124184358865223691e-15;
    tri3_xyze(1, 1) = 0.10936906572928781212;
    tri3_xyze(2, 1) = 0.84911436253628103099;
    nids.push_back(7113);
    tri3_xyze(0, 2) = 0.01375031769535093977;
    tri3_xyze(1, 2) = 0.10444403218566568337;
    tri3_xyze(2, 2) = 0.84950784947877078235;
    nids.push_back(-988);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 9.9147091733272353276e-15;
    tri3_xyze(1, 0) = 0.10313952597646713905;
    tri3_xyze(2, 0) = 0.84990133642126042268;
    nids.push_back(7115);
    tri3_xyze(0, 1) = 1.0005930993473840787e-14;
    tri3_xyze(1, 1) = 0.096860474023535966803;
    tri3_xyze(2, 1) = 0.84990133642126064473;
    nids.push_back(7117);
    tri3_xyze(0, 2) = 0.012940952255135975205;
    tri3_xyze(1, 2) = 0.098296291314454942167;
    tri3_xyze(2, 2) = 0.84990133642126064473;
    nids.push_back(-989);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005930993473840787e-14;
    tri3_xyze(1, 0) = 0.096860474023535966803;
    tri3_xyze(2, 0) = 0.84990133642126064473;
    nids.push_back(7117);
    tri3_xyze(0, 1) = 0.02506933539495863339;
    tri3_xyze(1, 1) = 0.093560033405934731054;
    tri3_xyze(2, 1) = 0.84990133642126064473;
    nids.push_back(7367);
    tri3_xyze(0, 2) = 0.012940952255135975205;
    tri3_xyze(1, 2) = 0.098296291314454942167;
    tri3_xyze(2, 2) = 0.84990133642126064473;
    nids.push_back(-989);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.02506933539495863339;
    tri3_xyze(1, 0) = 0.093560033405934731054;
    tri3_xyze(2, 0) = 0.84990133642126064473;
    nids.push_back(7367);
    tri3_xyze(0, 1) = 0.026694473625565345865;
    tri3_xyze(1, 1) = 0.099625131851881945644;
    tri3_xyze(2, 1) = 0.8499013364212605337;
    nids.push_back(7365);
    tri3_xyze(0, 2) = 0.012940952255135975205;
    tri3_xyze(1, 2) = 0.098296291314454942167;
    tri3_xyze(2, 2) = 0.84990133642126064473;
    nids.push_back(-989);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026694473625565345865;
    tri3_xyze(1, 0) = 0.099625131851881945644;
    tri3_xyze(2, 0) = 0.8499013364212605337;
    nids.push_back(7365);
    tri3_xyze(0, 1) = 9.9147091733272353276e-15;
    tri3_xyze(1, 1) = 0.10313952597646713905;
    tri3_xyze(2, 1) = 0.84990133642126042268;
    nids.push_back(7115);
    tri3_xyze(0, 2) = 0.012940952255135975205;
    tri3_xyze(1, 2) = 0.098296291314454942167;
    tri3_xyze(2, 2) = 0.84990133642126064473;
    nids.push_back(-989);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005930993473840787e-14;
    tri3_xyze(1, 0) = 0.096860474023535966803;
    tri3_xyze(2, 0) = 0.84990133642126064473;
    nids.push_back(7117);
    tri3_xyze(0, 1) = 1.0005549544177918912e-14;
    tri3_xyze(1, 1) = 0.090630934270715279855;
    tri3_xyze(2, 1) = 0.8491143625362815861;
    nids.push_back(7119);
    tri3_xyze(0, 2) = 0.012131586814921019313;
    tri3_xyze(1, 2) = 0.092148550443244228725;
    tri3_xyze(2, 2) = 0.84950784947877111541;
    nids.push_back(-990);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005549544177918912e-14;
    tri3_xyze(1, 0) = 0.090630934270715279855;
    tri3_xyze(2, 0) = 0.8491143625362815861;
    nids.push_back(7119);
    tri3_xyze(0, 1) = 0.023457011864705432092;
    tri3_xyze(1, 1) = 0.087542760072790937187;
    tri3_xyze(2, 1) = 0.84911436253628169712;
    nids.push_back(7369);
    tri3_xyze(0, 2) = 0.012131586814921019313;
    tri3_xyze(1, 2) = 0.092148550443244228725;
    tri3_xyze(2, 2) = 0.84950784947877111541;
    nids.push_back(-990);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.023457011864705432092;
    tri3_xyze(1, 0) = 0.087542760072790937187;
    tri3_xyze(2, 0) = 0.84911436253628169712;
    nids.push_back(7369);
    tri3_xyze(0, 1) = 0.02506933539495863339;
    tri3_xyze(1, 1) = 0.093560033405934731054;
    tri3_xyze(2, 1) = 0.84990133642126064473;
    nids.push_back(7367);
    tri3_xyze(0, 2) = 0.012131586814921019313;
    tri3_xyze(1, 2) = 0.092148550443244228725;
    tri3_xyze(2, 2) = 0.84950784947877111541;
    nids.push_back(-990);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.02506933539495863339;
    tri3_xyze(1, 0) = 0.093560033405934731054;
    tri3_xyze(2, 0) = 0.84990133642126064473;
    nids.push_back(7367);
    tri3_xyze(0, 1) = 1.0005930993473840787e-14;
    tri3_xyze(1, 1) = 0.096860474023535966803;
    tri3_xyze(2, 1) = 0.84990133642126064473;
    nids.push_back(7117);
    tri3_xyze(0, 2) = 0.012131586814921019313;
    tri3_xyze(1, 2) = 0.092148550443244228725;
    tri3_xyze(2, 2) = 0.84950784947877111541;
    nids.push_back(-990);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005549544177918912e-14;
    tri3_xyze(1, 0) = 0.090630934270715279855;
    tri3_xyze(2, 0) = 0.8491143625362815861;
    nids.push_back(7119);
    tri3_xyze(0, 1) = 1.0005177142313128246e-14;
    tri3_xyze(1, 1) = 0.084549150281254134565;
    tri3_xyze(2, 1) = 0.84755282581460511526;
    nids.push_back(7121);
    tri3_xyze(0, 2) = 0.011334985551189686453;
    tri3_xyze(1, 2) = 0.086097763118054837506;
    tri3_xyze(2, 2) = 0.84833359417544329517;
    nids.push_back(-991);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005177142313128246e-14;
    tri3_xyze(1, 0) = 0.084549150281254134565;
    tri3_xyze(2, 0) = 0.84755282581460511526;
    nids.push_back(7121);
    tri3_xyze(0, 1) = 0.021882930340033298483;
    tri3_xyze(1, 1) = 0.081668207847459012294;
    tri3_xyze(2, 1) = 0.84755282581460511526;
    nids.push_back(7371);
    tri3_xyze(0, 2) = 0.011334985551189686453;
    tri3_xyze(1, 2) = 0.086097763118054837506;
    tri3_xyze(2, 2) = 0.84833359417544329517;
    nids.push_back(-991);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.021882930340033298483;
    tri3_xyze(1, 0) = 0.081668207847459012294;
    tri3_xyze(2, 0) = 0.84755282581460511526;
    nids.push_back(7371);
    tri3_xyze(0, 1) = 0.023457011864705432092;
    tri3_xyze(1, 1) = 0.087542760072790937187;
    tri3_xyze(2, 1) = 0.84911436253628169712;
    nids.push_back(7369);
    tri3_xyze(0, 2) = 0.011334985551189686453;
    tri3_xyze(1, 2) = 0.086097763118054837506;
    tri3_xyze(2, 2) = 0.84833359417544329517;
    nids.push_back(-991);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.023457011864705432092;
    tri3_xyze(1, 0) = 0.087542760072790937187;
    tri3_xyze(2, 0) = 0.84911436253628169712;
    nids.push_back(7369);
    tri3_xyze(0, 1) = 1.0005549544177918912e-14;
    tri3_xyze(1, 1) = 0.090630934270715279855;
    tri3_xyze(2, 1) = 0.8491143625362815861;
    nids.push_back(7119);
    tri3_xyze(0, 2) = 0.011334985551189686453;
    tri3_xyze(1, 2) = 0.086097763118054837506;
    tri3_xyze(2, 2) = 0.84833359417544329517;
    nids.push_back(-991);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0005177142313128246e-14;
    tri3_xyze(1, 0) = 0.084549150281254134565;
    tri3_xyze(2, 0) = 0.84755282581460511526;
    nids.push_back(7121);
    tri3_xyze(0, 1) = 1.0004819660879340093e-14;
    tri3_xyze(1, 1) = 0.078711035421747829743;
    tri3_xyze(2, 1) = 0.84524135262314858963;
    nids.push_back(7123);
    tri3_xyze(0, 2) = 0.01056371134173758955;
    tri3_xyze(1, 2) = 0.080239353869570206879;
    tri3_xyze(2, 2) = 0.84639708921887679693;
    nids.push_back(-992);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004819660879340093e-14;
    tri3_xyze(1, 0) = 0.078711035421747829743;
    tri3_xyze(2, 0) = 0.84524135262314858963;
    nids.push_back(7123);
    tri3_xyze(0, 1) = 0.020371915026897051415;
    tri3_xyze(1, 1) = 0.076029021927819850912;
    tri3_xyze(2, 1) = 0.84524135262314858963;
    nids.push_back(7373);
    tri3_xyze(0, 2) = 0.01056371134173758955;
    tri3_xyze(1, 2) = 0.080239353869570206879;
    tri3_xyze(2, 2) = 0.84639708921887679693;
    nids.push_back(-992);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.020371915026897051415;
    tri3_xyze(1, 0) = 0.076029021927819850912;
    tri3_xyze(2, 0) = 0.84524135262314858963;
    nids.push_back(7373);
    tri3_xyze(0, 1) = 0.021882930340033298483;
    tri3_xyze(1, 1) = 0.081668207847459012294;
    tri3_xyze(2, 1) = 0.84755282581460511526;
    nids.push_back(7371);
    tri3_xyze(0, 2) = 0.01056371134173758955;
    tri3_xyze(1, 2) = 0.080239353869570206879;
    tri3_xyze(2, 2) = 0.84639708921887679693;
    nids.push_back(-992);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.021882930340033298483;
    tri3_xyze(1, 0) = 0.081668207847459012294;
    tri3_xyze(2, 0) = 0.84755282581460511526;
    nids.push_back(7371);
    tri3_xyze(0, 1) = 1.0005177142313128246e-14;
    tri3_xyze(1, 1) = 0.084549150281254134565;
    tri3_xyze(2, 1) = 0.84755282581460511526;
    nids.push_back(7121);
    tri3_xyze(0, 2) = 0.01056371134173758955;
    tri3_xyze(1, 2) = 0.080239353869570206879;
    tri3_xyze(2, 2) = 0.84639708921887679693;
    nids.push_back(-992);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004819660879340093e-14;
    tri3_xyze(1, 0) = 0.078711035421747829743;
    tri3_xyze(2, 0) = 0.84524135262314858963;
    nids.push_back(7123);
    tri3_xyze(0, 1) = 1.0069505089742456541e-14;
    tri3_xyze(1, 1) = 0.073208660251051493884;
    tri3_xyze(2, 1) = 0.84221639627494848668;
    nids.push_back(7125);
    tri3_xyze(0, 2) = 0.0098299276415847250399;
    tri3_xyze(1, 2) = 0.074665713311282938536;
    tri3_xyze(2, 2) = 0.84372887444904853815;
    nids.push_back(-993);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0069505089742456541e-14;
    tri3_xyze(1, 0) = 0.073208660251051493884;
    tri3_xyze(2, 0) = 0.84221639627494848668;
    nids.push_back(7125);
    tri3_xyze(0, 1) = 0.018947795539421774524;
    tri3_xyze(1, 1) = 0.070714135644512593482;
    tri3_xyze(2, 1) = 0.8422163962749485977;
    nids.push_back(7375);
    tri3_xyze(0, 2) = 0.0098299276415847250399;
    tri3_xyze(1, 2) = 0.074665713311282938536;
    tri3_xyze(2, 2) = 0.84372887444904853815;
    nids.push_back(-993);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.018947795539421774524;
    tri3_xyze(1, 0) = 0.070714135644512593482;
    tri3_xyze(2, 0) = 0.8422163962749485977;
    nids.push_back(7375);
    tri3_xyze(0, 1) = 0.020371915026897051415;
    tri3_xyze(1, 1) = 0.076029021927819850912;
    tri3_xyze(2, 1) = 0.84524135262314858963;
    nids.push_back(7373);
    tri3_xyze(0, 2) = 0.0098299276415847250399;
    tri3_xyze(1, 2) = 0.074665713311282938536;
    tri3_xyze(2, 2) = 0.84372887444904853815;
    nids.push_back(-993);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.020371915026897051415;
    tri3_xyze(1, 0) = 0.076029021927819850912;
    tri3_xyze(2, 0) = 0.84524135262314858963;
    nids.push_back(7373);
    tri3_xyze(0, 1) = 1.0004819660879340093e-14;
    tri3_xyze(1, 1) = 0.078711035421747829743;
    tri3_xyze(2, 1) = 0.84524135262314858963;
    nids.push_back(7123);
    tri3_xyze(0, 2) = 0.0098299276415847250399;
    tri3_xyze(1, 2) = 0.074665713311282938536;
    tri3_xyze(2, 2) = 0.84372887444904853815;
    nids.push_back(-993);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0069505089742456541e-14;
    tri3_xyze(1, 0) = 0.073208660251051493884;
    tri3_xyze(2, 0) = 0.84221639627494848668;
    nids.push_back(7125);
    tri3_xyze(0, 1) = 1.0004171685873872409e-14;
    tri3_xyze(1, 1) = 0.068128800512566770053;
    tri3_xyze(2, 1) = 0.83852566213863743982;
    nids.push_back(7127);
    tri3_xyze(0, 2) = 0.009145206658023544849;
    tri3_xyze(1, 2) = 0.06946474108432876271;
    tri3_xyze(2, 2) = 0.84037102920679296325;
    nids.push_back(-994);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004171685873872409e-14;
    tri3_xyze(1, 0) = 0.068128800512566770053;
    tri3_xyze(2, 0) = 0.83852566213863743982;
    nids.push_back(7127);
    tri3_xyze(0, 1) = 0.017633031092652327182;
    tri3_xyze(1, 1) = 0.065807367929184193422;
    tri3_xyze(2, 1) = 0.83852566213863743982;
    nids.push_back(7377);
    tri3_xyze(0, 2) = 0.009145206658023544849;
    tri3_xyze(1, 2) = 0.06946474108432876271;
    tri3_xyze(2, 2) = 0.84037102920679296325;
    nids.push_back(-994);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.017633031092652327182;
    tri3_xyze(1, 0) = 0.065807367929184193422;
    tri3_xyze(2, 0) = 0.83852566213863743982;
    nids.push_back(7377);
    tri3_xyze(0, 1) = 0.018947795539421774524;
    tri3_xyze(1, 1) = 0.070714135644512593482;
    tri3_xyze(2, 1) = 0.8422163962749485977;
    nids.push_back(7375);
    tri3_xyze(0, 2) = 0.009145206658023544849;
    tri3_xyze(1, 2) = 0.06946474108432876271;
    tri3_xyze(2, 2) = 0.84037102920679296325;
    nids.push_back(-994);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.018947795539421774524;
    tri3_xyze(1, 0) = 0.070714135644512593482;
    tri3_xyze(2, 0) = 0.8422163962749485977;
    nids.push_back(7375);
    tri3_xyze(0, 1) = 1.0069505089742456541e-14;
    tri3_xyze(1, 1) = 0.073208660251051493884;
    tri3_xyze(2, 1) = 0.84221639627494848668;
    nids.push_back(7125);
    tri3_xyze(0, 2) = 0.009145206658023544849;
    tri3_xyze(1, 2) = 0.06946474108432876271;
    tri3_xyze(2, 2) = 0.84037102920679296325;
    nids.push_back(-994);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0004171685873872409e-14;
    tri3_xyze(1, 0) = 0.068128800512566770053;
    tri3_xyze(2, 0) = 0.83852566213863743982;
    nids.push_back(7127);
    tri3_xyze(0, 1) = 1.0060336543049417096e-14;
    tri3_xyze(1, 1) = 0.063551568628930479865;
    tri3_xyze(2, 1) = 0.83422735529628244677;
    nids.push_back(7129);
    tri3_xyze(0, 2) = 0.0085203468499972980638;
    tri3_xyze(1, 2) = 0.064718459627636892861;
    tri3_xyze(2, 2) = 0.8363765087174599433;
    nids.push_back(-995);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0060336543049417096e-14;
    tri3_xyze(1, 0) = 0.063551568628930479865;
    tri3_xyze(2, 0) = 0.83422735529628244677;
    nids.push_back(7129);
    tri3_xyze(0, 1) = 0.016448356307316801261;
    tri3_xyze(1, 1) = 0.061386101439866114227;
    tri3_xyze(2, 1) = 0.83422735529628255779;
    nids.push_back(7379);
    tri3_xyze(0, 2) = 0.0085203468499972980638;
    tri3_xyze(1, 2) = 0.064718459627636892861;
    tri3_xyze(2, 2) = 0.8363765087174599433;
    nids.push_back(-995);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.016448356307316801261;
    tri3_xyze(1, 0) = 0.061386101439866114227;
    tri3_xyze(2, 0) = 0.83422735529628255779;
    nids.push_back(7379);
    tri3_xyze(0, 1) = 0.017633031092652327182;
    tri3_xyze(1, 1) = 0.065807367929184193422;
    tri3_xyze(2, 1) = 0.83852566213863743982;
    nids.push_back(7377);
    tri3_xyze(0, 2) = 0.0085203468499972980638;
    tri3_xyze(1, 2) = 0.064718459627636892861;
    tri3_xyze(2, 2) = 0.8363765087174599433;
    nids.push_back(-995);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.017633031092652327182;
    tri3_xyze(1, 0) = 0.065807367929184193422;
    tri3_xyze(2, 0) = 0.83852566213863743982;
    nids.push_back(7377);
    tri3_xyze(0, 1) = 1.0004171685873872409e-14;
    tri3_xyze(1, 1) = 0.068128800512566770053;
    tri3_xyze(2, 1) = 0.83852566213863743982;
    nids.push_back(7127);
    tri3_xyze(0, 2) = 0.0085203468499972980638;
    tri3_xyze(1, 2) = 0.064718459627636892861;
    tri3_xyze(2, 2) = 0.8363765087174599433;
    nids.push_back(-995);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0060336543049417096e-14;
    tri3_xyze(1, 0) = 0.063551568628930479865;
    tri3_xyze(2, 0) = 0.83422735529628244677;
    nids.push_back(7129);
    tri3_xyze(0, 1) = 1.010942687419677122e-14;
    tri3_xyze(1, 1) = 0.059549150281253564188;
    tri3_xyze(2, 1) = 0.82938926261447187827;
    nids.push_back(7131);
    tri3_xyze(0, 2) = 0.007965202629951819227;
    tri3_xyze(1, 2) = 0.0605017206350704842;
    tri3_xyze(2, 2) = 0.83180830895537716252;
    nids.push_back(-996);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.010942687419677122e-14;
    tri3_xyze(1, 0) = 0.059549150281253564188;
    tri3_xyze(2, 0) = 0.82938926261447187827;
    nids.push_back(7131);
    tri3_xyze(0, 1) = 0.015412454212470304282;
    tri3_xyze(1, 1) = 0.057520062190231792398;
    tri3_xyze(2, 1) = 0.82938926261447198929;
    nids.push_back(7381);
    tri3_xyze(0, 2) = 0.007965202629951819227;
    tri3_xyze(1, 2) = 0.0605017206350704842;
    tri3_xyze(2, 2) = 0.83180830895537716252;
    nids.push_back(-996);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.015412454212470304282;
    tri3_xyze(1, 0) = 0.057520062190231792398;
    tri3_xyze(2, 0) = 0.82938926261447198929;
    nids.push_back(7381);
    tri3_xyze(0, 1) = 0.016448356307316801261;
    tri3_xyze(1, 1) = 0.061386101439866114227;
    tri3_xyze(2, 1) = 0.83422735529628255779;
    nids.push_back(7379);
    tri3_xyze(0, 2) = 0.007965202629951819227;
    tri3_xyze(1, 2) = 0.0605017206350704842;
    tri3_xyze(2, 2) = 0.83180830895537716252;
    nids.push_back(-996);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.016448356307316801261;
    tri3_xyze(1, 0) = 0.061386101439866114227;
    tri3_xyze(2, 0) = 0.83422735529628255779;
    nids.push_back(7379);
    tri3_xyze(0, 1) = 1.0060336543049417096e-14;
    tri3_xyze(1, 1) = 0.063551568628930479865;
    tri3_xyze(2, 1) = 0.83422735529628244677;
    nids.push_back(7129);
    tri3_xyze(0, 2) = 0.007965202629951819227;
    tri3_xyze(1, 2) = 0.0605017206350704842;
    tri3_xyze(2, 2) = 0.83180830895537716252;
    nids.push_back(-996);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.010942687419677122e-14;
    tri3_xyze(1, 0) = 0.059549150281253564188;
    tri3_xyze(2, 0) = 0.82938926261447187827;
    nids.push_back(7131);
    tri3_xyze(0, 1) = 1.0003440318568768715e-14;
    tri3_xyze(1, 1) = 0.056184665997807595184;
    tri3_xyze(2, 1) = 0.82408768370493412636;
    nids.push_back(7133);
    tri3_xyze(0, 2) = 0.0074885289538642093449;
    tri3_xyze(1, 2) = 0.056881024599500146222;
    tri3_xyze(2, 2) = 0.8267384731597029468;
    nids.push_back(-997);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003440318568768715e-14;
    tri3_xyze(1, 0) = 0.056184665997807595184;
    tri3_xyze(2, 0) = 0.82408768370493412636;
    nids.push_back(7133);
    tri3_xyze(0, 1) = 0.014541661602966418979;
    tri3_xyze(1, 1) = 0.054270219928707640056;
    tri3_xyze(2, 1) = 0.82408768370493412636;
    nids.push_back(7383);
    tri3_xyze(0, 2) = 0.0074885289538642093449;
    tri3_xyze(1, 2) = 0.056881024599500146222;
    tri3_xyze(2, 2) = 0.8267384731597029468;
    nids.push_back(-997);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.014541661602966418979;
    tri3_xyze(1, 0) = 0.054270219928707640056;
    tri3_xyze(2, 0) = 0.82408768370493412636;
    nids.push_back(7383);
    tri3_xyze(0, 1) = 0.015412454212470304282;
    tri3_xyze(1, 1) = 0.057520062190231792398;
    tri3_xyze(2, 1) = 0.82938926261447198929;
    nids.push_back(7381);
    tri3_xyze(0, 2) = 0.0074885289538642093449;
    tri3_xyze(1, 2) = 0.056881024599500146222;
    tri3_xyze(2, 2) = 0.8267384731597029468;
    nids.push_back(-997);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.015412454212470304282;
    tri3_xyze(1, 0) = 0.057520062190231792398;
    tri3_xyze(2, 0) = 0.82938926261447198929;
    nids.push_back(7381);
    tri3_xyze(0, 1) = 1.010942687419677122e-14;
    tri3_xyze(1, 1) = 0.059549150281253564188;
    tri3_xyze(2, 1) = 0.82938926261447187827;
    nids.push_back(7131);
    tri3_xyze(0, 2) = 0.0074885289538642093449;
    tri3_xyze(1, 2) = 0.056881024599500146222;
    tri3_xyze(2, 2) = 0.8267384731597029468;
    nids.push_back(-997);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003440318568768715e-14;
    tri3_xyze(1, 0) = 0.056184665997807595184;
    tri3_xyze(2, 0) = 0.82408768370493412636;
    nids.push_back(7133);
    tri3_xyze(0, 1) = 1.0003276614502322774e-14;
    tri3_xyze(1, 1) = 0.053511175705588047424;
    tri3_xyze(2, 1) = 0.81840622763408232831;
    nids.push_back(7135);
    tri3_xyze(0, 2) = 0.0070978432503574417564;
    tri3_xyze(1, 2) = 0.053913472060305729638;
    tri3_xyze(2, 2) = 0.82124695566950822734;
    nids.push_back(-998);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 1.0003276614502322774e-14;
    tri3_xyze(1, 0) = 0.053511175705588047424;
    tri3_xyze(2, 0) = 0.81840622763408232831;
    nids.push_back(7135);
    tri3_xyze(0, 1) = 0.013849711398443341481;
    tri3_xyze(1, 1) = 0.051687826609119649768;
    tri3_xyze(2, 1) = 0.81840622763408232831;
    nids.push_back(7385);
    tri3_xyze(0, 2) = 0.0070978432503574417564;
    tri3_xyze(1, 2) = 0.053913472060305729638;
    tri3_xyze(2, 2) = 0.82124695566950822734;
    nids.push_back(-998);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013849711398443341481;
    tri3_xyze(1, 0) = 0.051687826609119649768;
    tri3_xyze(2, 0) = 0.81840622763408232831;
    nids.push_back(7385);
    tri3_xyze(0, 1) = 0.014541661602966418979;
    tri3_xyze(1, 1) = 0.054270219928707640056;
    tri3_xyze(2, 1) = 0.82408768370493412636;
    nids.push_back(7383);
    tri3_xyze(0, 2) = 0.0070978432503574417564;
    tri3_xyze(1, 2) = 0.053913472060305729638;
    tri3_xyze(2, 2) = 0.82124695566950822734;
    nids.push_back(-998);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.014541661602966418979;
    tri3_xyze(1, 0) = 0.054270219928707640056;
    tri3_xyze(2, 0) = 0.82408768370493412636;
    nids.push_back(7383);
    tri3_xyze(0, 1) = 1.0003440318568768715e-14;
    tri3_xyze(1, 1) = 0.056184665997807595184;
    tri3_xyze(2, 1) = 0.82408768370493412636;
    nids.push_back(7133);
    tri3_xyze(0, 2) = 0.0070978432503574417564;
    tri3_xyze(1, 2) = 0.053913472060305729638;
    tri3_xyze(2, 2) = 0.82124695566950822734;
    nids.push_back(-998);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013347516066977408236;
    tri3_xyze(1, 0) = 0.04981360811676393413;
    tri3_xyze(2, 0) = 0.78756550564160565653;
    nids.push_back(7295);
    tri3_xyze(0, 1) = 0.013849711398443344951;
    tri3_xyze(1, 1) = 0.051687826609118511789;
    tri3_xyze(2, 1) = 0.78159377236561444757;
    nids.push_back(7297);
    tri3_xyze(0, 2) = 0.019934559072504666466;
    tri3_xyze(1, 2) = 0.048126282872743814112;
    tri3_xyze(2, 2) = 0.78457963900361016307;
    nids.push_back(-1004);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013849711398443344951;
    tri3_xyze(1, 0) = 0.051687826609118511789;
    tri3_xyze(2, 0) = 0.78159377236561444757;
    nids.push_back(7297);
    tri3_xyze(0, 1) = 0.026755587852803703469;
    tri3_xyze(1, 1) = 0.046342037547410869225;
    tri3_xyze(2, 1) = 0.78159377236561466962;
    nids.push_back(7547);
    tri3_xyze(0, 2) = 0.019934559072504666466;
    tri3_xyze(1, 2) = 0.048126282872743814112;
    tri3_xyze(2, 2) = 0.78457963900361016307;
    nids.push_back(-1004);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013849711398443344951;
    tri3_xyze(1, 0) = 0.051687826609118511789;
    tri3_xyze(2, 0) = 0.78159377236561444757;
    nids.push_back(7297);
    tri3_xyze(0, 1) = 0.014541661602966469286;
    tri3_xyze(1, 1) = 0.054270219928706113499;
    tri3_xyze(2, 1) = 0.7759123162947625385;
    nids.push_back(7299);
    tri3_xyze(0, 2) = 0.020809823463281737488;
    tri3_xyze(1, 2) = 0.050239358035619824383;
    tri3_xyze(2, 2) = 0.77875304433018854855;
    nids.push_back(-1005);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.014541661602966469286;
    tri3_xyze(1, 0) = 0.054270219928706113499;
    tri3_xyze(2, 0) = 0.7759123162947625385;
    nids.push_back(7299);
    tri3_xyze(0, 1) = 0.028092332998913432246;
    tri3_xyze(1, 1) = 0.048657348057243816897;
    tri3_xyze(2, 1) = 0.7759123162947625385;
    nids.push_back(7549);
    tri3_xyze(0, 2) = 0.020809823463281737488;
    tri3_xyze(1, 2) = 0.050239358035619824383;
    tri3_xyze(2, 2) = 0.77875304433018854855;
    nids.push_back(-1005);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.028092332998913432246;
    tri3_xyze(1, 0) = 0.048657348057243816897;
    tri3_xyze(2, 0) = 0.7759123162947625385;
    nids.push_back(7549);
    tri3_xyze(0, 1) = 0.026755587852803703469;
    tri3_xyze(1, 1) = 0.046342037547410869225;
    tri3_xyze(2, 1) = 0.78159377236561466962;
    nids.push_back(7547);
    tri3_xyze(0, 2) = 0.020809823463281737488;
    tri3_xyze(1, 2) = 0.050239358035619824383;
    tri3_xyze(2, 2) = 0.77875304433018854855;
    nids.push_back(-1005);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026755587852803703469;
    tri3_xyze(1, 0) = 0.046342037547410869225;
    tri3_xyze(2, 0) = 0.78159377236561466962;
    nids.push_back(7547);
    tri3_xyze(0, 1) = 0.013849711398443344951;
    tri3_xyze(1, 1) = 0.051687826609118511789;
    tri3_xyze(2, 1) = 0.78159377236561444757;
    nids.push_back(7297);
    tri3_xyze(0, 2) = 0.020809823463281737488;
    tri3_xyze(1, 2) = 0.050239358035619824383;
    tri3_xyze(2, 2) = 0.77875304433018854855;
    nids.push_back(-1005);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.014541661602966469286;
    tri3_xyze(1, 0) = 0.054270219928706113499;
    tri3_xyze(2, 0) = 0.7759123162947625385;
    nids.push_back(7299);
    tri3_xyze(0, 1) = 0.015412454212470300813;
    tri3_xyze(1, 1) = 0.057520062190229932775;
    tri3_xyze(2, 1) = 0.77061073738522456456;
    nids.push_back(7301);
    tri3_xyze(0, 2) = 0.021955255988746612272;
    tri3_xyze(1, 2) = 0.053004676773380249299;
    tri3_xyze(2, 2) = 0.77326152683999360704;
    nids.push_back(-1006);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.015412454212470300813;
    tri3_xyze(1, 0) = 0.057520062190229932775;
    tri3_xyze(2, 0) = 0.77061073738522456456;
    nids.push_back(7301);
    tri3_xyze(0, 1) = 0.029774575140636246745;
    tri3_xyze(1, 1) = 0.051571076917341147905;
    tri3_xyze(2, 1) = 0.77061073738522467558;
    nids.push_back(7551);
    tri3_xyze(0, 2) = 0.021955255988746612272;
    tri3_xyze(1, 2) = 0.053004676773380249299;
    tri3_xyze(2, 2) = 0.77326152683999360704;
    nids.push_back(-1006);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.029774575140636246745;
    tri3_xyze(1, 0) = 0.051571076917341147905;
    tri3_xyze(2, 0) = 0.77061073738522467558;
    nids.push_back(7551);
    tri3_xyze(0, 1) = 0.028092332998913432246;
    tri3_xyze(1, 1) = 0.048657348057243816897;
    tri3_xyze(2, 1) = 0.7759123162947625385;
    nids.push_back(7549);
    tri3_xyze(0, 2) = 0.021955255988746612272;
    tri3_xyze(1, 2) = 0.053004676773380249299;
    tri3_xyze(2, 2) = 0.77326152683999360704;
    nids.push_back(-1006);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.028092332998913432246;
    tri3_xyze(1, 0) = 0.048657348057243816897;
    tri3_xyze(2, 0) = 0.7759123162947625385;
    nids.push_back(7549);
    tri3_xyze(0, 1) = 0.014541661602966469286;
    tri3_xyze(1, 1) = 0.054270219928706113499;
    tri3_xyze(2, 1) = 0.7759123162947625385;
    nids.push_back(7299);
    tri3_xyze(0, 2) = 0.021955255988746612272;
    tri3_xyze(1, 2) = 0.053004676773380249299;
    tri3_xyze(2, 2) = 0.77326152683999360704;
    nids.push_back(-1006);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.015412454212470300813;
    tri3_xyze(1, 0) = 0.057520062190229932775;
    tri3_xyze(2, 0) = 0.77061073738522456456;
    nids.push_back(7301);
    tri3_xyze(0, 1) = 0.016448356307316835956;
    tri3_xyze(1, 1) = 0.061386101439864095008;
    tri3_xyze(2, 1) = 0.76577264470341366298;
    nids.push_back(7303);
    tri3_xyze(0, 2) = 0.023352792493724513351;
    tri3_xyze(1, 2) = 0.056378628357609289479;
    tri3_xyze(2, 2) = 0.76819169104431916928;
    nids.push_back(-1007);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.016448356307316835956;
    tri3_xyze(1, 0) = 0.061386101439864095008;
    tri3_xyze(2, 0) = 0.76577264470341366298;
    nids.push_back(7303);
    tri3_xyze(0, 1) = 0.031775784314474669889;
    tri3_xyze(1, 1) = 0.055037272883001996104;
    tri3_xyze(2, 1) = 0.76577264470341388503;
    nids.push_back(7553);
    tri3_xyze(0, 2) = 0.023352792493724513351;
    tri3_xyze(1, 2) = 0.056378628357609289479;
    tri3_xyze(2, 2) = 0.76819169104431916928;
    nids.push_back(-1007);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.031775784314474669889;
    tri3_xyze(1, 0) = 0.055037272883001996104;
    tri3_xyze(2, 0) = 0.76577264470341388503;
    nids.push_back(7553);
    tri3_xyze(0, 1) = 0.029774575140636246745;
    tri3_xyze(1, 1) = 0.051571076917341147905;
    tri3_xyze(2, 1) = 0.77061073738522467558;
    nids.push_back(7551);
    tri3_xyze(0, 2) = 0.023352792493724513351;
    tri3_xyze(1, 2) = 0.056378628357609289479;
    tri3_xyze(2, 2) = 0.76819169104431916928;
    nids.push_back(-1007);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.029774575140636246745;
    tri3_xyze(1, 0) = 0.051571076917341147905;
    tri3_xyze(2, 0) = 0.77061073738522467558;
    nids.push_back(7551);
    tri3_xyze(0, 1) = 0.015412454212470300813;
    tri3_xyze(1, 1) = 0.057520062190229932775;
    tri3_xyze(2, 1) = 0.77061073738522456456;
    nids.push_back(7301);
    tri3_xyze(0, 2) = 0.023352792493724513351;
    tri3_xyze(1, 2) = 0.056378628357609289479;
    tri3_xyze(2, 2) = 0.76819169104431916928;
    nids.push_back(-1007);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.016448356307316835956;
    tri3_xyze(1, 0) = 0.061386101439864095008;
    tri3_xyze(2, 0) = 0.76577264470341366298;
    nids.push_back(7303);
    tri3_xyze(0, 1) = 0.017633031092652306365;
    tri3_xyze(1, 1) = 0.065807367929181695421;
    tri3_xyze(2, 1) = 0.76147433786105855891;
    nids.push_back(7305);
    tri3_xyze(0, 2) = 0.024980392992684126441;
    tri3_xyze(1, 2) = 0.060308003556322667404;
    tri3_xyze(2, 2) = 0.76362349128223627748;
    nids.push_back(-1008);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.017633031092652306365;
    tri3_xyze(1, 0) = 0.065807367929181695421;
    tri3_xyze(2, 0) = 0.76147433786105855891;
    nids.push_back(7305);
    tri3_xyze(0, 1) = 0.034064400256292683145;
    tri3_xyze(1, 1) = 0.059001271973242869207;
    tri3_xyze(2, 1) = 0.76147433786105866993;
    nids.push_back(7555);
    tri3_xyze(0, 2) = 0.024980392992684126441;
    tri3_xyze(1, 2) = 0.060308003556322667404;
    tri3_xyze(2, 2) = 0.76362349128223627748;
    nids.push_back(-1008);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.034064400256292683145;
    tri3_xyze(1, 0) = 0.059001271973242869207;
    tri3_xyze(2, 0) = 0.76147433786105866993;
    nids.push_back(7555);
    tri3_xyze(0, 1) = 0.031775784314474669889;
    tri3_xyze(1, 1) = 0.055037272883001996104;
    tri3_xyze(2, 1) = 0.76577264470341388503;
    nids.push_back(7553);
    tri3_xyze(0, 2) = 0.024980392992684126441;
    tri3_xyze(1, 2) = 0.060308003556322667404;
    tri3_xyze(2, 2) = 0.76362349128223627748;
    nids.push_back(-1008);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.031775784314474669889;
    tri3_xyze(1, 0) = 0.055037272883001996104;
    tri3_xyze(2, 0) = 0.76577264470341388503;
    nids.push_back(7553);
    tri3_xyze(0, 1) = 0.016448356307316835956;
    tri3_xyze(1, 1) = 0.061386101439864095008;
    tri3_xyze(2, 1) = 0.76577264470341366298;
    nids.push_back(7303);
    tri3_xyze(0, 2) = 0.024980392992684126441;
    tri3_xyze(1, 2) = 0.060308003556322667404;
    tri3_xyze(2, 2) = 0.76362349128223627748;
    nids.push_back(-1008);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.017633031092652306365;
    tri3_xyze(1, 0) = 0.065807367929181695421;
    tri3_xyze(2, 0) = 0.76147433786105855891;
    nids.push_back(7305);
    tri3_xyze(0, 1) = 0.018947795539421625338;
    tri3_xyze(1, 1) = 0.070714135644510123235;
    tri3_xyze(2, 1) = 0.75778360372474695694;
    nids.push_back(7307);
    tri3_xyze(0, 2) = 0.026812389253475420181;
    tri3_xyze(1, 2) = 0.064730833775341750647;
    tri3_xyze(2, 2) = 0.75962897079290292446;
    nids.push_back(-1009);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.018947795539421625338;
    tri3_xyze(1, 0) = 0.070714135644510123235;
    tri3_xyze(2, 0) = 0.75778360372474695694;
    nids.push_back(7307);
    tri3_xyze(0, 1) = 0.036604330125535065876;
    tri3_xyze(1, 1) = 0.063400559554432342479;
    tri3_xyze(2, 1) = 0.75778360372474717899;
    nids.push_back(7557);
    tri3_xyze(0, 2) = 0.026812389253475420181;
    tri3_xyze(1, 2) = 0.064730833775341750647;
    tri3_xyze(2, 2) = 0.75962897079290292446;
    nids.push_back(-1009);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036604330125535065876;
    tri3_xyze(1, 0) = 0.063400559554432342479;
    tri3_xyze(2, 0) = 0.75778360372474717899;
    nids.push_back(7557);
    tri3_xyze(0, 1) = 0.034064400256292683145;
    tri3_xyze(1, 1) = 0.059001271973242869207;
    tri3_xyze(2, 1) = 0.76147433786105866993;
    nids.push_back(7555);
    tri3_xyze(0, 2) = 0.026812389253475420181;
    tri3_xyze(1, 2) = 0.064730833775341750647;
    tri3_xyze(2, 2) = 0.75962897079290292446;
    nids.push_back(-1009);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.034064400256292683145;
    tri3_xyze(1, 0) = 0.059001271973242869207;
    tri3_xyze(2, 0) = 0.76147433786105866993;
    nids.push_back(7555);
    tri3_xyze(0, 1) = 0.017633031092652306365;
    tri3_xyze(1, 1) = 0.065807367929181695421;
    tri3_xyze(2, 1) = 0.76147433786105855891;
    nids.push_back(7305);
    tri3_xyze(0, 2) = 0.026812389253475420181;
    tri3_xyze(1, 2) = 0.064730833775341750647;
    tri3_xyze(2, 2) = 0.75962897079290292446;
    nids.push_back(-1009);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.018947795539421625338;
    tri3_xyze(1, 0) = 0.070714135644510123235;
    tri3_xyze(2, 0) = 0.75778360372474695694;
    nids.push_back(7307);
    tri3_xyze(0, 1) = 0.020371915026897093048;
    tri3_xyze(1, 1) = 0.076029021927817172499;
    tri3_xyze(2, 1) = 0.75475864737654663195;
    nids.push_back(7309);
    tri3_xyze(0, 2) = 0.028819889600684225894;
    tri3_xyze(1, 2) = 0.069577368340041842454;
    tri3_xyze(2, 2) = 0.75627112555064690547;
    nids.push_back(-1010);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.020371915026897093048;
    tri3_xyze(1, 0) = 0.076029021927817172499;
    tri3_xyze(2, 0) = 0.75475864737654663195;
    nids.push_back(7309);
    tri3_xyze(0, 1) = 0.039355517710883122784;
    tri3_xyze(1, 1) = 0.068165756233407717724;
    tri3_xyze(2, 1) = 0.75475864737654696501;
    nids.push_back(7559);
    tri3_xyze(0, 2) = 0.028819889600684225894;
    tri3_xyze(1, 2) = 0.069577368340041842454;
    tri3_xyze(2, 2) = 0.75627112555064690547;
    nids.push_back(-1010);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.039355517710883122784;
    tri3_xyze(1, 0) = 0.068165756233407717724;
    tri3_xyze(2, 0) = 0.75475864737654696501;
    nids.push_back(7559);
    tri3_xyze(0, 1) = 0.036604330125535065876;
    tri3_xyze(1, 1) = 0.063400559554432342479;
    tri3_xyze(2, 1) = 0.75778360372474717899;
    nids.push_back(7557);
    tri3_xyze(0, 2) = 0.028819889600684225894;
    tri3_xyze(1, 2) = 0.069577368340041842454;
    tri3_xyze(2, 2) = 0.75627112555064690547;
    nids.push_back(-1010);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036604330125535065876;
    tri3_xyze(1, 0) = 0.063400559554432342479;
    tri3_xyze(2, 0) = 0.75778360372474717899;
    nids.push_back(7557);
    tri3_xyze(0, 1) = 0.018947795539421625338;
    tri3_xyze(1, 1) = 0.070714135644510123235;
    tri3_xyze(2, 1) = 0.75778360372474695694;
    nids.push_back(7307);
    tri3_xyze(0, 2) = 0.028819889600684225894;
    tri3_xyze(1, 2) = 0.069577368340041842454;
    tri3_xyze(2, 2) = 0.75627112555064690547;
    nids.push_back(-1010);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.020371915026897093048;
    tri3_xyze(1, 0) = 0.076029021927817172499;
    tri3_xyze(2, 0) = 0.75475864737654663195;
    nids.push_back(7309);
    tri3_xyze(0, 1) = 0.021882930340033308891;
    tri3_xyze(1, 1) = 0.081668207847456056325;
    tri3_xyze(2, 1) = 0.75244717418508977325;
    nids.push_back(7311);
    tri3_xyze(0, 2) = 0.030971234554612484674;
    tri3_xyze(1, 2) = 0.074771174505158091717;
    tri3_xyze(2, 2) = 0.75360291078081831362;
    nids.push_back(-1011);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.021882930340033308891;
    tri3_xyze(1, 0) = 0.081668207847456056325;
    tri3_xyze(2, 0) = 0.75244717418508977325;
    nids.push_back(7311);
    tri3_xyze(0, 1) = 0.042274575140636413972;
    tri3_xyze(1, 1) = 0.073221712011951434196;
    tri3_xyze(2, 1) = 0.75244717418508999529;
    nids.push_back(7561);
    tri3_xyze(0, 2) = 0.030971234554612484674;
    tri3_xyze(1, 2) = 0.074771174505158091717;
    tri3_xyze(2, 2) = 0.75360291078081831362;
    nids.push_back(-1011);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042274575140636413972;
    tri3_xyze(1, 0) = 0.073221712011951434196;
    tri3_xyze(2, 0) = 0.75244717418508999529;
    nids.push_back(7561);
    tri3_xyze(0, 1) = 0.039355517710883122784;
    tri3_xyze(1, 1) = 0.068165756233407717724;
    tri3_xyze(2, 1) = 0.75475864737654696501;
    nids.push_back(7559);
    tri3_xyze(0, 2) = 0.030971234554612484674;
    tri3_xyze(1, 2) = 0.074771174505158091717;
    tri3_xyze(2, 2) = 0.75360291078081831362;
    nids.push_back(-1011);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.039355517710883122784;
    tri3_xyze(1, 0) = 0.068165756233407717724;
    tri3_xyze(2, 0) = 0.75475864737654696501;
    nids.push_back(7559);
    tri3_xyze(0, 1) = 0.020371915026897093048;
    tri3_xyze(1, 1) = 0.076029021927817172499;
    tri3_xyze(2, 1) = 0.75475864737654663195;
    nids.push_back(7309);
    tri3_xyze(0, 2) = 0.030971234554612484674;
    tri3_xyze(1, 2) = 0.074771174505158091717;
    tri3_xyze(2, 2) = 0.75360291078081831362;
    nids.push_back(-1011);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.021882930340033308891;
    tri3_xyze(1, 0) = 0.081668207847456056325;
    tri3_xyze(2, 0) = 0.75244717418508977325;
    nids.push_back(7311);
    tri3_xyze(0, 1) = 0.023457011864705432092;
    tri3_xyze(1, 1) = 0.087542760072787856318;
    tri3_xyze(2, 1) = 0.75088563746341285832;
    nids.push_back(7313);
    tri3_xyze(0, 2) = 0.033232496120185472943;
    tri3_xyze(1, 2) = 0.080230342844837418403;
    tri3_xyze(2, 2) = 0.75166640582425137129;
    nids.push_back(-1012);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.023457011864705432092;
    tri3_xyze(1, 0) = 0.087542760072787856318;
    tri3_xyze(2, 0) = 0.75088563746341285832;
    nids.push_back(7313);
    tri3_xyze(0, 1) = 0.045315467135366757634;
    tri3_xyze(1, 1) = 0.078488691447154299019;
    tri3_xyze(2, 1) = 0.75088563746341319138;
    nids.push_back(7563);
    tri3_xyze(0, 2) = 0.033232496120185472943;
    tri3_xyze(1, 2) = 0.080230342844837418403;
    tri3_xyze(2, 2) = 0.75166640582425137129;
    nids.push_back(-1012);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.045315467135366757634;
    tri3_xyze(1, 0) = 0.078488691447154299019;
    tri3_xyze(2, 0) = 0.75088563746341319138;
    nids.push_back(7563);
    tri3_xyze(0, 1) = 0.042274575140636413972;
    tri3_xyze(1, 1) = 0.073221712011951434196;
    tri3_xyze(2, 1) = 0.75244717418508999529;
    nids.push_back(7561);
    tri3_xyze(0, 2) = 0.033232496120185472943;
    tri3_xyze(1, 2) = 0.080230342844837418403;
    tri3_xyze(2, 2) = 0.75166640582425137129;
    nids.push_back(-1012);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042274575140636413972;
    tri3_xyze(1, 0) = 0.073221712011951434196;
    tri3_xyze(2, 0) = 0.75244717418508999529;
    nids.push_back(7561);
    tri3_xyze(0, 1) = 0.021882930340033308891;
    tri3_xyze(1, 1) = 0.081668207847456056325;
    tri3_xyze(2, 1) = 0.75244717418508977325;
    nids.push_back(7311);
    tri3_xyze(0, 2) = 0.033232496120185472943;
    tri3_xyze(1, 2) = 0.080230342844837418403;
    tri3_xyze(2, 2) = 0.75166640582425137129;
    nids.push_back(-1012);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.023457011864705432092;
    tri3_xyze(1, 0) = 0.087542760072787856318;
    tri3_xyze(2, 0) = 0.75088563746341285832;
    nids.push_back(7313);
    tri3_xyze(0, 1) = 0.025069335394958522367;
    tri3_xyze(1, 1) = 0.093560033405931525285;
    tri3_xyze(2, 1) = 0.75009866357843346663;
    nids.push_back(7315);
    tri3_xyze(0, 2) = 0.035568012851701955035;
    tri3_xyze(1, 2) = 0.08586877901321390838;
    tri3_xyze(2, 2) = 0.75049215052092344003;
    nids.push_back(-1013);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025069335394958522367;
    tri3_xyze(1, 0) = 0.093560033405931525285;
    tri3_xyze(2, 0) = 0.75009866357843346663;
    nids.push_back(7315);
    tri3_xyze(0, 1) = 0.048430237011777121925;
    tri3_xyze(1, 1) = 0.083883631126981925141;
    tri3_xyze(2, 1) = 0.75009866357843379969;
    nids.push_back(7565);
    tri3_xyze(0, 2) = 0.035568012851701955035;
    tri3_xyze(1, 2) = 0.08586877901321390838;
    tri3_xyze(2, 2) = 0.75049215052092344003;
    nids.push_back(-1013);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048430237011777121925;
    tri3_xyze(1, 0) = 0.083883631126981925141;
    tri3_xyze(2, 0) = 0.75009866357843379969;
    nids.push_back(7565);
    tri3_xyze(0, 1) = 0.045315467135366757634;
    tri3_xyze(1, 1) = 0.078488691447154299019;
    tri3_xyze(2, 1) = 0.75088563746341319138;
    nids.push_back(7563);
    tri3_xyze(0, 2) = 0.035568012851701955035;
    tri3_xyze(1, 2) = 0.08586877901321390838;
    tri3_xyze(2, 2) = 0.75049215052092344003;
    nids.push_back(-1013);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.045315467135366757634;
    tri3_xyze(1, 0) = 0.078488691447154299019;
    tri3_xyze(2, 0) = 0.75088563746341319138;
    nids.push_back(7563);
    tri3_xyze(0, 1) = 0.023457011864705432092;
    tri3_xyze(1, 1) = 0.087542760072787856318;
    tri3_xyze(2, 1) = 0.75088563746341285832;
    nids.push_back(7313);
    tri3_xyze(0, 2) = 0.035568012851701955035;
    tri3_xyze(1, 2) = 0.08586877901321390838;
    tri3_xyze(2, 2) = 0.75049215052092344003;
    nids.push_back(-1013);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025069335394958522367;
    tri3_xyze(1, 0) = 0.093560033405931525285;
    tri3_xyze(2, 0) = 0.75009866357843346663;
    nids.push_back(7315);
    tri3_xyze(0, 1) = 0.026694473625565474234;
    tri3_xyze(1, 1) = 0.099625131851878920286;
    tri3_xyze(2, 1) = 0.7500986635784333556;
    nids.push_back(7317);
    tri3_xyze(0, 2) = 0.037940952255135981797;
    tri3_xyze(1, 2) = 0.091597561503673791128;
    tri3_xyze(2, 2) = 0.75009866357843357765;
    nids.push_back(-1014);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026694473625565474234;
    tri3_xyze(1, 0) = 0.099625131851878920286;
    tri3_xyze(2, 0) = 0.7500986635784333556;
    nids.push_back(7317);
    tri3_xyze(0, 1) = 0.051569762988242798252;
    tri3_xyze(1, 1) = 0.089321449629902766043;
    tri3_xyze(2, 1) = 0.75009866357843368867;
    nids.push_back(7567);
    tri3_xyze(0, 2) = 0.037940952255135981797;
    tri3_xyze(1, 2) = 0.091597561503673791128;
    tri3_xyze(2, 2) = 0.75009866357843357765;
    nids.push_back(-1014);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051569762988242798252;
    tri3_xyze(1, 0) = 0.089321449629902766043;
    tri3_xyze(2, 0) = 0.75009866357843368867;
    nids.push_back(7567);
    tri3_xyze(0, 1) = 0.048430237011777121925;
    tri3_xyze(1, 1) = 0.083883631126981925141;
    tri3_xyze(2, 1) = 0.75009866357843379969;
    nids.push_back(7565);
    tri3_xyze(0, 2) = 0.037940952255135981797;
    tri3_xyze(1, 2) = 0.091597561503673791128;
    tri3_xyze(2, 2) = 0.75009866357843357765;
    nids.push_back(-1014);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048430237011777121925;
    tri3_xyze(1, 0) = 0.083883631126981925141;
    tri3_xyze(2, 0) = 0.75009866357843379969;
    nids.push_back(7565);
    tri3_xyze(0, 1) = 0.025069335394958522367;
    tri3_xyze(1, 1) = 0.093560033405931525285;
    tri3_xyze(2, 1) = 0.75009866357843346663;
    nids.push_back(7315);
    tri3_xyze(0, 2) = 0.037940952255135981797;
    tri3_xyze(1, 2) = 0.091597561503673791128;
    tri3_xyze(2, 2) = 0.75009866357843357765;
    nids.push_back(-1014);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026694473625565474234;
    tri3_xyze(1, 0) = 0.099625131851878920286;
    tri3_xyze(2, 0) = 0.7500986635784333556;
    nids.push_back(7317);
    tri3_xyze(0, 1) = 0.028306797155818672063;
    tri3_xyze(1, 1) = 0.10564240518502271415;
    tri3_xyze(2, 1) = 0.75088563746341230321;
    nids.push_back(7319);
    tri3_xyze(0, 2) = 0.040313891658570008558;
    tri3_xyze(1, 2) = 0.097326343994133729387;
    tri3_xyze(2, 2) = 0.75049215052092299594;
    nids.push_back(-1015);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.028306797155818672063;
    tri3_xyze(1, 0) = 0.10564240518502271415;
    tri3_xyze(2, 0) = 0.75088563746341230321;
    nids.push_back(7319);
    tri3_xyze(0, 1) = 0.054684532864653093154;
    tri3_xyze(1, 1) = 0.094716389309730475432;
    tri3_xyze(2, 1) = 0.75088563746341263627;
    nids.push_back(7569);
    tri3_xyze(0, 2) = 0.040313891658570008558;
    tri3_xyze(1, 2) = 0.097326343994133729387;
    tri3_xyze(2, 2) = 0.75049215052092299594;
    nids.push_back(-1015);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.054684532864653093154;
    tri3_xyze(1, 0) = 0.094716389309730475432;
    tri3_xyze(2, 0) = 0.75088563746341263627;
    nids.push_back(7569);
    tri3_xyze(0, 1) = 0.051569762988242798252;
    tri3_xyze(1, 1) = 0.089321449629902766043;
    tri3_xyze(2, 1) = 0.75009866357843368867;
    nids.push_back(7567);
    tri3_xyze(0, 2) = 0.040313891658570008558;
    tri3_xyze(1, 2) = 0.097326343994133729387;
    tri3_xyze(2, 2) = 0.75049215052092299594;
    nids.push_back(-1015);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051569762988242798252;
    tri3_xyze(1, 0) = 0.089321449629902766043;
    tri3_xyze(2, 0) = 0.75009866357843368867;
    nids.push_back(7567);
    tri3_xyze(0, 1) = 0.026694473625565474234;
    tri3_xyze(1, 1) = 0.099625131851878920286;
    tri3_xyze(2, 1) = 0.7500986635784333556;
    nids.push_back(7317);
    tri3_xyze(0, 2) = 0.040313891658570008558;
    tri3_xyze(1, 2) = 0.097326343994133729387;
    tri3_xyze(2, 2) = 0.75049215052092299594;
    nids.push_back(-1015);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.057725424859383582532;
    tri3_xyze(1, 0) = 0.0999833687449331876;
    tri3_xyze(2, 0) = 0.75244717418508921813;
    nids.push_back(7571);
    tri3_xyze(0, 1) = 0.054684532864653093154;
    tri3_xyze(1, 1) = 0.094716389309730475432;
    tri3_xyze(2, 1) = 0.75088563746341263627;
    nids.push_back(7569);
    tri3_xyze(0, 2) = 0.042649408390086532283;
    tri3_xyze(1, 2) = 0.10296478016251023324;
    tri3_xyze(2, 2) = 0.75166640582425070516;
    nids.push_back(-1016);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.054684532864653093154;
    tri3_xyze(1, 0) = 0.094716389309730475432;
    tri3_xyze(2, 0) = 0.75088563746341263627;
    nids.push_back(7569);
    tri3_xyze(0, 1) = 0.028306797155818672063;
    tri3_xyze(1, 1) = 0.10564240518502271415;
    tri3_xyze(2, 1) = 0.75088563746341230321;
    nids.push_back(7319);
    tri3_xyze(0, 2) = 0.042649408390086532283;
    tri3_xyze(1, 2) = 0.10296478016251023324;
    tri3_xyze(2, 2) = 0.75166640582425070516;
    nids.push_back(-1016);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.038822856765388082112;
    tri3_xyze(1, 0) = 0.14488887394336025238;
    tri3_xyze(2, 0) = 0.79999999999984550136;
    nids.push_back(7341);
    tri3_xyze(0, 1) = 0.038720813491581330956;
    tri3_xyze(1, 1) = 0.14450804326094299213;
    tri3_xyze(2, 1) = 0.80626666167806060326;
    nids.push_back(7343);
    tri3_xyze(0, 2) = 0.056836634447462869302;
    tri3_xyze(1, 2) = 0.13721577372268262129;
    tri3_xyze(2, 2) = 0.80313333083895321884;
    nids.push_back(-1027);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.038720813491581330956;
    tri3_xyze(1, 0) = 0.14450804326094299213;
    tri3_xyze(2, 0) = 0.80626666167806060326;
    nids.push_back(7343);
    tri3_xyze(0, 1) = 0.07480286753287186674;
    tri3_xyze(1, 1) = 0.12956236711876151224;
    tri3_xyze(2, 1) = 0.80626666167806115837;
    nids.push_back(7593);
    tri3_xyze(0, 2) = 0.056836634447462869302;
    tri3_xyze(1, 2) = 0.13721577372268262129;
    tri3_xyze(2, 2) = 0.80313333083895321884;
    nids.push_back(-1027);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.075000000000010183521;
    tri3_xyze(1, 0) = 0.12990381056766570067;
    tri3_xyze(2, 0) = 0.79999999999984594545;
    nids.push_back(7591);
    tri3_xyze(0, 1) = 0.038822856765388082112;
    tri3_xyze(1, 1) = 0.14488887394336025238;
    tri3_xyze(2, 1) = 0.79999999999984550136;
    nids.push_back(7341);
    tri3_xyze(0, 2) = 0.056836634447462869302;
    tri3_xyze(1, 2) = 0.13721577372268262129;
    tri3_xyze(2, 2) = 0.80313333083895321884;
    nids.push_back(-1027);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.038720813491581330956;
    tri3_xyze(1, 0) = 0.14450804326094299213;
    tri3_xyze(2, 0) = 0.80626666167806060326;
    nids.push_back(7343);
    tri3_xyze(0, 1) = 0.038416292953546704592;
    tri3_xyze(1, 1) = 0.14337155714104973803;
    tri3_xyze(2, 1) = 0.81243449435808823278;
    nids.push_back(7345);
    tri3_xyze(0, 2) = 0.056538638251556375125;
    tri3_xyze(1, 2) = 0.13649634726499004778;
    tri3_xyze(2, 2) = 0.80935057801807464006;
    nids.push_back(-1028);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.038416292953546704592;
    tri3_xyze(1, 0) = 0.14337155714104973803;
    tri3_xyze(2, 0) = 0.81243449435808823278;
    nids.push_back(7345);
    tri3_xyze(0, 1) = 0.074214579028225605151;
    tri3_xyze(1, 1) = 0.1285434215392058932;
    tri3_xyze(2, 1) = 0.81243449435808878789;
    nids.push_back(7595);
    tri3_xyze(0, 2) = 0.056538638251556375125;
    tri3_xyze(1, 2) = 0.13649634726499004778;
    tri3_xyze(2, 2) = 0.80935057801807464006;
    nids.push_back(-1028);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.07480286753287186674;
    tri3_xyze(1, 0) = 0.12956236711876151224;
    tri3_xyze(2, 0) = 0.80626666167806115837;
    nids.push_back(7593);
    tri3_xyze(0, 1) = 0.038720813491581330956;
    tri3_xyze(1, 1) = 0.14450804326094299213;
    tri3_xyze(2, 1) = 0.80626666167806060326;
    nids.push_back(7343);
    tri3_xyze(0, 2) = 0.056538638251556375125;
    tri3_xyze(1, 2) = 0.13649634726499004778;
    tri3_xyze(2, 2) = 0.80935057801807464006;
    nids.push_back(-1028);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.038416292953546704592;
    tri3_xyze(1, 0) = 0.14337155714104973803;
    tri3_xyze(2, 0) = 0.81243449435808823278;
    nids.push_back(7345);
    tri3_xyze(0, 1) = 0.037914097622080769612;
    tri3_xyze(1, 1) = 0.14149733864869518118;
    tri3_xyze(2, 1) = 0.81840622763407955276;
    nids.push_back(7347);
    tri3_xyze(0, 2) = 0.055947345437767317944;
    tri3_xyze(1, 2) = 0.13506884013460696004;
    tri3_xyze(2, 2) = 0.81542036099608417032;
    nids.push_back(-1029);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.037914097622080769612;
    tri3_xyze(1, 0) = 0.14149733864869518118;
    tri3_xyze(2, 0) = 0.81840622763407955276;
    nids.push_back(7347);
    tri3_xyze(0, 1) = 0.073244412147216220177;
    tri3_xyze(1, 1) = 0.12686304320947697222;
    tri3_xyze(2, 1) = 0.81840622763407988582;
    nids.push_back(7597);
    tri3_xyze(0, 2) = 0.055947345437767317944;
    tri3_xyze(1, 2) = 0.13506884013460696004;
    tri3_xyze(2, 2) = 0.81542036099608417032;
    nids.push_back(-1029);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.074214579028225605151;
    tri3_xyze(1, 0) = 0.1285434215392058932;
    tri3_xyze(2, 0) = 0.81243449435808878789;
    nids.push_back(7595);
    tri3_xyze(0, 1) = 0.038416292953546704592;
    tri3_xyze(1, 1) = 0.14337155714104973803;
    tri3_xyze(2, 1) = 0.81243449435808823278;
    nids.push_back(7345);
    tri3_xyze(0, 2) = 0.055947345437767317944;
    tri3_xyze(1, 2) = 0.13506884013460696004;
    tri3_xyze(2, 2) = 0.81542036099608417032;
    nids.push_back(-1029);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.037914097622080769612;
    tri3_xyze(1, 0) = 0.14149733864869518118;
    tri3_xyze(2, 0) = 0.81840622763407955276;
    nids.push_back(7347);
    tri3_xyze(0, 1) = 0.03722214741755768691;
    tri3_xyze(1, 1) = 0.1389149453291075309;
    tri3_xyze(2, 1) = 0.82408768370493157285;
    nids.push_back(7349);
    tri3_xyze(0, 2) = 0.055072081046990295494;
    tri3_xyze(1, 2) = 0.13295576497173094976;
    tri3_xyze(2, 2) = 0.82124695566950578485;
    nids.push_back(-1030);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.03722214741755768691;
    tri3_xyze(1, 0) = 0.1389149453291075309;
    tri3_xyze(2, 0) = 0.82408768370493157285;
    nids.push_back(7349);
    tri3_xyze(0, 1) = 0.071907667001106526095;
    tri3_xyze(1, 1) = 0.12454773269964407312;
    tri3_xyze(2, 1) = 0.82408768370493190591;
    nids.push_back(7599);
    tri3_xyze(0, 2) = 0.055072081046990295494;
    tri3_xyze(1, 2) = 0.13295576497173094976;
    tri3_xyze(2, 2) = 0.82124695566950578485;
    nids.push_back(-1030);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.073244412147216220177;
    tri3_xyze(1, 0) = 0.12686304320947697222;
    tri3_xyze(2, 0) = 0.81840622763407988582;
    nids.push_back(7597);
    tri3_xyze(0, 1) = 0.037914097622080769612;
    tri3_xyze(1, 1) = 0.14149733864869518118;
    tri3_xyze(2, 1) = 0.81840622763407955276;
    nids.push_back(7347);
    tri3_xyze(0, 2) = 0.055072081046990295494;
    tri3_xyze(1, 2) = 0.13295576497173094976;
    tri3_xyze(2, 2) = 0.82124695566950578485;
    nids.push_back(-1030);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.03722214741755768691;
    tri3_xyze(1, 0) = 0.1389149453291075309;
    tri3_xyze(2, 0) = 0.82408768370493157285;
    nids.push_back(7349);
    tri3_xyze(0, 1) = 0.036351354808053792933;
    tri3_xyze(1, 1) = 0.13566510306758366999;
    tri3_xyze(2, 1) = 0.8293892626144695468;
    nids.push_back(7351);
    tri3_xyze(0, 2) = 0.053926648521525399893;
    tri3_xyze(1, 2) = 0.13019044623397046934;
    tri3_xyze(2, 2) = 0.82673847315970061533;
    nids.push_back(-1031);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036351354808053792933;
    tri3_xyze(1, 0) = 0.13566510306758366999;
    tri3_xyze(2, 0) = 0.8293892626144695468;
    nids.push_back(7351);
    tri3_xyze(0, 1) = 0.070225424859383614451;
    tri3_xyze(1, 1) = 0.12163400383954667272;
    tri3_xyze(2, 1) = 0.82938926261446987986;
    nids.push_back(7601);
    tri3_xyze(0, 2) = 0.053926648521525399893;
    tri3_xyze(1, 2) = 0.13019044623397046934;
    tri3_xyze(2, 2) = 0.82673847315970061533;
    nids.push_back(-1031);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.071907667001106526095;
    tri3_xyze(1, 0) = 0.12454773269964407312;
    tri3_xyze(2, 0) = 0.82408768370493190591;
    nids.push_back(7599);
    tri3_xyze(0, 1) = 0.03722214741755768691;
    tri3_xyze(1, 1) = 0.1389149453291075309;
    tri3_xyze(2, 1) = 0.82408768370493157285;
    nids.push_back(7349);
    tri3_xyze(0, 2) = 0.053926648521525399893;
    tri3_xyze(1, 2) = 0.13019044623397046934;
    tri3_xyze(2, 2) = 0.82673847315970061533;
    nids.push_back(-1031);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036351354808053792933;
    tri3_xyze(1, 0) = 0.13566510306758366999;
    tri3_xyze(2, 0) = 0.8293892626144695468;
    nids.push_back(7351);
    tri3_xyze(0, 1) = 0.035315452713207295954;
    tri3_xyze(1, 1) = 0.13179906381794967429;
    tri3_xyze(2, 1) = 0.83422735529628033735;
    nids.push_back(7353);
    tri3_xyze(0, 2) = 0.052529112016547488406;
    tri3_xyze(1, 2) = 0.12681649464974142916;
    tri3_xyze(2, 2) = 0.83180830895537516412;
    nids.push_back(-1032);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.035315452713207295954;
    tri3_xyze(1, 0) = 0.13179906381794967429;
    tri3_xyze(2, 0) = 0.83422735529628033735;
    nids.push_back(7353);
    tri3_xyze(0, 1) = 0.068224215685545236409;
    tri3_xyze(1, 1) = 0.11816780787388576901;
    tri3_xyze(2, 1) = 0.83422735529628067042;
    nids.push_back(7603);
    tri3_xyze(0, 2) = 0.052529112016547488406;
    tri3_xyze(1, 2) = 0.12681649464974142916;
    tri3_xyze(2, 2) = 0.83180830895537516412;
    nids.push_back(-1032);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.070225424859383614451;
    tri3_xyze(1, 0) = 0.12163400383954667272;
    tri3_xyze(2, 0) = 0.82938926261446987986;
    nids.push_back(7601);
    tri3_xyze(0, 1) = 0.036351354808053792933;
    tri3_xyze(1, 1) = 0.13566510306758366999;
    tri3_xyze(2, 1) = 0.8293892626144695468;
    nids.push_back(7351);
    tri3_xyze(0, 2) = 0.052529112016547488406;
    tri3_xyze(1, 2) = 0.12681649464974142916;
    tri3_xyze(2, 2) = 0.83180830895537516412;
    nids.push_back(-1032);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.035315452713207295954;
    tri3_xyze(1, 0) = 0.13179906381794967429;
    tri3_xyze(2, 0) = 0.83422735529628033735;
    nids.push_back(7353);
    tri3_xyze(0, 1) = 0.034130777927871766564;
    tri3_xyze(1, 1) = 0.12737779732863185878;
    tri3_xyze(2, 1) = 0.83852566213863555245;
    nids.push_back(7355);
    tri3_xyze(0, 2) = 0.050901511517587837152;
    tri3_xyze(1, 2) = 0.12288711945102806511;
    tri3_xyze(2, 2) = 0.83637650871745816694;
    nids.push_back(-1033);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.034130777927871766564;
    tri3_xyze(1, 0) = 0.12737779732863185878;
    tri3_xyze(2, 0) = 0.83852566213863555245;
    nids.push_back(7355);
    tri3_xyze(0, 1) = 0.065935599743727049682;
    tri3_xyze(1, 1) = 0.11420380878364494448;
    tri3_xyze(2, 1) = 0.83852566213863599653;
    nids.push_back(7605);
    tri3_xyze(0, 2) = 0.050901511517587837152;
    tri3_xyze(1, 2) = 0.12288711945102806511;
    tri3_xyze(2, 2) = 0.83637650871745816694;
    nids.push_back(-1033);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.068224215685545236409;
    tri3_xyze(1, 0) = 0.11816780787388576901;
    tri3_xyze(2, 0) = 0.83422735529628067042;
    nids.push_back(7603);
    tri3_xyze(0, 1) = 0.035315452713207295954;
    tri3_xyze(1, 1) = 0.13179906381794967429;
    tri3_xyze(2, 1) = 0.83422735529628033735;
    nids.push_back(7353);
    tri3_xyze(0, 2) = 0.050901511517587837152;
    tri3_xyze(1, 2) = 0.12288711945102806511;
    tri3_xyze(2, 2) = 0.83637650871745816694;
    nids.push_back(-1033);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.034130777927871766564;
    tri3_xyze(1, 0) = 0.12737779732863185878;
    tri3_xyze(2, 0) = 0.83852566213863555245;
    nids.push_back(7355);
    tri3_xyze(0, 1) = 0.032816013481102468408;
    tri3_xyze(1, 1) = 0.12247102961330375015;
    tri3_xyze(2, 1) = 0.84221639627494693237;
    nids.push_back(7357);
    tri3_xyze(0, 2) = 0.049069515256796536473;
    tri3_xyze(1, 2) = 0.11846428923200903738;
    tri3_xyze(2, 2) = 0.84037102920679140894;
    nids.push_back(-1034);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.032816013481102468408;
    tri3_xyze(1, 0) = 0.12247102961330375015;
    tri3_xyze(2, 0) = 0.84221639627494693237;
    nids.push_back(7357);
    tri3_xyze(0, 1) = 0.063395669874484861239;
    tri3_xyze(1, 1) = 0.10980452120245559611;
    tri3_xyze(2, 1) = 0.84221639627494737645;
    nids.push_back(7607);
    tri3_xyze(0, 2) = 0.049069515256796536473;
    tri3_xyze(1, 2) = 0.11846428923200903738;
    tri3_xyze(2, 2) = 0.84037102920679140894;
    nids.push_back(-1034);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.063395669874484861239;
    tri3_xyze(1, 0) = 0.10980452120245559611;
    tri3_xyze(2, 0) = 0.84221639627494737645;
    nids.push_back(7607);
    tri3_xyze(0, 1) = 0.065935599743727049682;
    tri3_xyze(1, 1) = 0.11420380878364494448;
    tri3_xyze(2, 1) = 0.83852566213863599653;
    nids.push_back(7605);
    tri3_xyze(0, 2) = 0.049069515256796536473;
    tri3_xyze(1, 2) = 0.11846428923200903738;
    tri3_xyze(2, 2) = 0.84037102920679140894;
    nids.push_back(-1034);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.065935599743727049682;
    tri3_xyze(1, 0) = 0.11420380878364494448;
    tri3_xyze(2, 0) = 0.83852566213863599653;
    nids.push_back(7605);
    tri3_xyze(0, 1) = 0.034130777927871766564;
    tri3_xyze(1, 1) = 0.12737779732863185878;
    tri3_xyze(2, 1) = 0.83852566213863555245;
    nids.push_back(7355);
    tri3_xyze(0, 2) = 0.049069515256796536473;
    tri3_xyze(1, 2) = 0.11846428923200903738;
    tri3_xyze(2, 2) = 0.84037102920679140894;
    nids.push_back(-1034);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.032816013481102468408;
    tri3_xyze(1, 0) = 0.12247102961330375015;
    tri3_xyze(2, 0) = 0.84221639627494693237;
    nids.push_back(7357);
    tri3_xyze(0, 1) = 0.031391893993627052739;
    tri3_xyze(1, 1) = 0.11715614332999664537;
    tri3_xyze(2, 1) = 0.84524135262314736838;
    nids.push_back(7359);
    tri3_xyze(0, 2) = 0.047062014909587779332;
    tri3_xyze(1, 2) = 0.1136177546673089872;
    tri3_xyze(2, 2) = 0.84372887444904742793;
    nids.push_back(-1035);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.031391893993627052739;
    tri3_xyze(1, 0) = 0.11715614332999664537;
    tri3_xyze(2, 0) = 0.84524135262314736838;
    nids.push_back(7359);
    tri3_xyze(0, 1) = 0.060644482289136734943;
    tri3_xyze(1, 1) = 0.10503932452347999882;
    tri3_xyze(2, 1) = 0.84524135262314770145;
    nids.push_back(7609);
    tri3_xyze(0, 2) = 0.047062014909587779332;
    tri3_xyze(1, 2) = 0.1136177546673089872;
    tri3_xyze(2, 2) = 0.84372887444904742793;
    nids.push_back(-1035);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.060644482289136734943;
    tri3_xyze(1, 0) = 0.10503932452347999882;
    tri3_xyze(2, 0) = 0.84524135262314770145;
    nids.push_back(7609);
    tri3_xyze(0, 1) = 0.063395669874484861239;
    tri3_xyze(1, 1) = 0.10980452120245559611;
    tri3_xyze(2, 1) = 0.84221639627494737645;
    nids.push_back(7607);
    tri3_xyze(0, 2) = 0.047062014909587779332;
    tri3_xyze(1, 2) = 0.1136177546673089872;
    tri3_xyze(2, 2) = 0.84372887444904742793;
    nids.push_back(-1035);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.063395669874484861239;
    tri3_xyze(1, 0) = 0.10980452120245559611;
    tri3_xyze(2, 0) = 0.84221639627494737645;
    nids.push_back(7607);
    tri3_xyze(0, 1) = 0.032816013481102468408;
    tri3_xyze(1, 1) = 0.12247102961330375015;
    tri3_xyze(2, 1) = 0.84221639627494693237;
    nids.push_back(7357);
    tri3_xyze(0, 2) = 0.047062014909587779332;
    tri3_xyze(1, 2) = 0.1136177546673089872;
    tri3_xyze(2, 2) = 0.84372887444904742793;
    nids.push_back(-1035);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.031391893993627052739;
    tri3_xyze(1, 0) = 0.11715614332999664537;
    tri3_xyze(2, 0) = 0.84524135262314736838;
    nids.push_back(7359);
    tri3_xyze(0, 1) = 0.029880878680490805671;
    tri3_xyze(1, 1) = 0.11151695741035763665;
    tri3_xyze(2, 1) = 0.84755282581460411606;
    nids.push_back(7361);
    tri3_xyze(0, 2) = 0.044910669955659558716;
    tri3_xyze(1, 2) = 0.1084239485021926408;
    tri3_xyze(2, 2) = 0.84639708921887590876;
    nids.push_back(-1036);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.029880878680490805671;
    tri3_xyze(1, 0) = 0.11151695741035763665;
    tri3_xyze(2, 0) = 0.84755282581460411606;
    nids.push_back(7361);
    tri3_xyze(0, 1) = 0.057725424859383624165;
    tri3_xyze(1, 1) = 0.099983368744936254591;
    tri3_xyze(2, 1) = 0.84755282581460467117;
    nids.push_back(7611);
    tri3_xyze(0, 2) = 0.044910669955659558716;
    tri3_xyze(1, 2) = 0.1084239485021926408;
    tri3_xyze(2, 2) = 0.84639708921887590876;
    nids.push_back(-1036);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.057725424859383624165;
    tri3_xyze(1, 0) = 0.099983368744936254591;
    tri3_xyze(2, 0) = 0.84755282581460467117;
    nids.push_back(7611);
    tri3_xyze(0, 1) = 0.060644482289136734943;
    tri3_xyze(1, 1) = 0.10503932452347999882;
    tri3_xyze(2, 1) = 0.84524135262314770145;
    nids.push_back(7609);
    tri3_xyze(0, 2) = 0.044910669955659558716;
    tri3_xyze(1, 2) = 0.1084239485021926408;
    tri3_xyze(2, 2) = 0.84639708921887590876;
    nids.push_back(-1036);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.060644482289136734943;
    tri3_xyze(1, 0) = 0.10503932452347999882;
    tri3_xyze(2, 0) = 0.84524135262314770145;
    nids.push_back(7609);
    tri3_xyze(0, 1) = 0.031391893993627052739;
    tri3_xyze(1, 1) = 0.11715614332999664537;
    tri3_xyze(2, 1) = 0.84524135262314736838;
    nids.push_back(7359);
    tri3_xyze(0, 2) = 0.044910669955659558716;
    tri3_xyze(1, 2) = 0.1084239485021926408;
    tri3_xyze(2, 2) = 0.84639708921887590876;
    nids.push_back(-1036);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.029880878680490805671;
    tri3_xyze(1, 0) = 0.11151695741035763665;
    tri3_xyze(2, 0) = 0.84755282581460411606;
    nids.push_back(7361);
    tri3_xyze(0, 1) = 0.02830679715581868594;
    tri3_xyze(1, 1) = 0.10564240518502583666;
    tri3_xyze(2, 1) = 0.84911436253628114201;
    nids.push_back(7363);
    tri3_xyze(0, 2) = 0.042649408390086546161;
    tri3_xyze(1, 2) = 0.10296478016251331411;
    tri3_xyze(2, 2) = 0.84833359417544285108;
    nids.push_back(-1037);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.02830679715581868594;
    tri3_xyze(1, 0) = 0.10564240518502583666;
    tri3_xyze(2, 0) = 0.84911436253628114201;
    nids.push_back(7363);
    tri3_xyze(0, 1) = 0.054684532864653072337;
    tri3_xyze(1, 1) = 0.094716389309733514668;
    tri3_xyze(2, 1) = 0.84911436253628147508;
    nids.push_back(7613);
    tri3_xyze(0, 2) = 0.042649408390086546161;
    tri3_xyze(1, 2) = 0.10296478016251331411;
    tri3_xyze(2, 2) = 0.84833359417544285108;
    nids.push_back(-1037);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.054684532864653072337;
    tri3_xyze(1, 0) = 0.094716389309733514668;
    tri3_xyze(2, 0) = 0.84911436253628147508;
    nids.push_back(7613);
    tri3_xyze(0, 1) = 0.057725424859383624165;
    tri3_xyze(1, 1) = 0.099983368744936254591;
    tri3_xyze(2, 1) = 0.84755282581460467117;
    nids.push_back(7611);
    tri3_xyze(0, 2) = 0.042649408390086546161;
    tri3_xyze(1, 2) = 0.10296478016251331411;
    tri3_xyze(2, 2) = 0.84833359417544285108;
    nids.push_back(-1037);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.057725424859383624165;
    tri3_xyze(1, 0) = 0.099983368744936254591;
    tri3_xyze(2, 0) = 0.84755282581460467117;
    nids.push_back(7611);
    tri3_xyze(0, 1) = 0.029880878680490805671;
    tri3_xyze(1, 1) = 0.11151695741035763665;
    tri3_xyze(2, 1) = 0.84755282581460411606;
    nids.push_back(7361);
    tri3_xyze(0, 2) = 0.042649408390086546161;
    tri3_xyze(1, 2) = 0.10296478016251331411;
    tri3_xyze(2, 2) = 0.84833359417544285108;
    nids.push_back(-1037);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.02830679715581868594;
    tri3_xyze(1, 0) = 0.10564240518502583666;
    tri3_xyze(2, 0) = 0.84911436253628114201;
    nids.push_back(7363);
    tri3_xyze(0, 1) = 0.026694473625565345865;
    tri3_xyze(1, 1) = 0.099625131851881945644;
    tri3_xyze(2, 1) = 0.8499013364212605337;
    nids.push_back(7365);
    tri3_xyze(0, 2) = 0.040313891658569966925;
    tri3_xyze(1, 2) = 0.097326343994136796378;
    tri3_xyze(2, 2) = 0.84950784947877100439;
    nids.push_back(-1038);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026694473625565345865;
    tri3_xyze(1, 0) = 0.099625131851881945644;
    tri3_xyze(2, 0) = 0.8499013364212605337;
    nids.push_back(7365);
    tri3_xyze(0, 1) = 0.051569762988242770496;
    tri3_xyze(1, 1) = 0.089321449629905846912;
    tri3_xyze(2, 1) = 0.84990133642126086677;
    nids.push_back(7615);
    tri3_xyze(0, 2) = 0.040313891658569966925;
    tri3_xyze(1, 2) = 0.097326343994136796378;
    tri3_xyze(2, 2) = 0.84950784947877100439;
    nids.push_back(-1038);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051569762988242770496;
    tri3_xyze(1, 0) = 0.089321449629905846912;
    tri3_xyze(2, 0) = 0.84990133642126086677;
    nids.push_back(7615);
    tri3_xyze(0, 1) = 0.054684532864653072337;
    tri3_xyze(1, 1) = 0.094716389309733514668;
    tri3_xyze(2, 1) = 0.84911436253628147508;
    nids.push_back(7613);
    tri3_xyze(0, 2) = 0.040313891658569966925;
    tri3_xyze(1, 2) = 0.097326343994136796378;
    tri3_xyze(2, 2) = 0.84950784947877100439;
    nids.push_back(-1038);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.054684532864653072337;
    tri3_xyze(1, 0) = 0.094716389309733514668;
    tri3_xyze(2, 0) = 0.84911436253628147508;
    nids.push_back(7613);
    tri3_xyze(0, 1) = 0.02830679715581868594;
    tri3_xyze(1, 1) = 0.10564240518502583666;
    tri3_xyze(2, 1) = 0.84911436253628114201;
    nids.push_back(7363);
    tri3_xyze(0, 2) = 0.040313891658569966925;
    tri3_xyze(1, 2) = 0.097326343994136796378;
    tri3_xyze(2, 2) = 0.84950784947877100439;
    nids.push_back(-1038);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026694473625565345865;
    tri3_xyze(1, 0) = 0.099625131851881945644;
    tri3_xyze(2, 0) = 0.8499013364212605337;
    nids.push_back(7365);
    tri3_xyze(0, 1) = 0.02506933539495863339;
    tri3_xyze(1, 1) = 0.093560033405934731054;
    tri3_xyze(2, 1) = 0.84990133642126064473;
    nids.push_back(7367);
    tri3_xyze(0, 2) = 0.037940952255135967919;
    tri3_xyze(1, 2) = 0.091597561503676885875;
    tri3_xyze(2, 2) = 0.84990133642126086677;
    nids.push_back(-1039);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.02506933539495863339;
    tri3_xyze(1, 0) = 0.093560033405934731054;
    tri3_xyze(2, 0) = 0.84990133642126064473;
    nids.push_back(7367);
    tri3_xyze(0, 1) = 0.048430237011777101108;
    tri3_xyze(1, 1) = 0.083883631126985019888;
    tri3_xyze(2, 1) = 0.84990133642126097779;
    nids.push_back(7617);
    tri3_xyze(0, 2) = 0.037940952255135967919;
    tri3_xyze(1, 2) = 0.091597561503676885875;
    tri3_xyze(2, 2) = 0.84990133642126086677;
    nids.push_back(-1039);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048430237011777101108;
    tri3_xyze(1, 0) = 0.083883631126985019888;
    tri3_xyze(2, 0) = 0.84990133642126097779;
    nids.push_back(7617);
    tri3_xyze(0, 1) = 0.051569762988242770496;
    tri3_xyze(1, 1) = 0.089321449629905846912;
    tri3_xyze(2, 1) = 0.84990133642126086677;
    nids.push_back(7615);
    tri3_xyze(0, 2) = 0.037940952255135967919;
    tri3_xyze(1, 2) = 0.091597561503676885875;
    tri3_xyze(2, 2) = 0.84990133642126086677;
    nids.push_back(-1039);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051569762988242770496;
    tri3_xyze(1, 0) = 0.089321449629905846912;
    tri3_xyze(2, 0) = 0.84990133642126086677;
    nids.push_back(7615);
    tri3_xyze(0, 1) = 0.026694473625565345865;
    tri3_xyze(1, 1) = 0.099625131851881945644;
    tri3_xyze(2, 1) = 0.8499013364212605337;
    nids.push_back(7365);
    tri3_xyze(0, 2) = 0.037940952255135967919;
    tri3_xyze(1, 2) = 0.091597561503676885875;
    tri3_xyze(2, 2) = 0.84990133642126086677;
    nids.push_back(-1039);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.02506933539495863339;
    tri3_xyze(1, 0) = 0.093560033405934731054;
    tri3_xyze(2, 0) = 0.84990133642126064473;
    nids.push_back(7367);
    tri3_xyze(0, 1) = 0.023457011864705432092;
    tri3_xyze(1, 1) = 0.087542760072790937187;
    tri3_xyze(2, 1) = 0.84911436253628169712;
    nids.push_back(7369);
    tri3_xyze(0, 2) = 0.035568012851701996668;
    tri3_xyze(1, 2) = 0.085868779013217017004;
    tri3_xyze(2, 2) = 0.84950784947877133746;
    nids.push_back(-1040);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.023457011864705432092;
    tri3_xyze(1, 0) = 0.087542760072790937187;
    tri3_xyze(2, 0) = 0.84911436253628169712;
    nids.push_back(7369);
    tri3_xyze(0, 1) = 0.045315467135366820084;
    tri3_xyze(1, 1) = 0.078488691447157338255;
    tri3_xyze(2, 1) = 0.84911436253628203019;
    nids.push_back(7619);
    tri3_xyze(0, 2) = 0.035568012851701996668;
    tri3_xyze(1, 2) = 0.085868779013217017004;
    tri3_xyze(2, 2) = 0.84950784947877133746;
    nids.push_back(-1040);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.045315467135366820084;
    tri3_xyze(1, 0) = 0.078488691447157338255;
    tri3_xyze(2, 0) = 0.84911436253628203019;
    nids.push_back(7619);
    tri3_xyze(0, 1) = 0.048430237011777101108;
    tri3_xyze(1, 1) = 0.083883631126985019888;
    tri3_xyze(2, 1) = 0.84990133642126097779;
    nids.push_back(7617);
    tri3_xyze(0, 2) = 0.035568012851701996668;
    tri3_xyze(1, 2) = 0.085868779013217017004;
    tri3_xyze(2, 2) = 0.84950784947877133746;
    nids.push_back(-1040);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048430237011777101108;
    tri3_xyze(1, 0) = 0.083883631126985019888;
    tri3_xyze(2, 0) = 0.84990133642126097779;
    nids.push_back(7617);
    tri3_xyze(0, 1) = 0.02506933539495863339;
    tri3_xyze(1, 1) = 0.093560033405934731054;
    tri3_xyze(2, 1) = 0.84990133642126064473;
    nids.push_back(7367);
    tri3_xyze(0, 2) = 0.035568012851701996668;
    tri3_xyze(1, 2) = 0.085868779013217017004;
    tri3_xyze(2, 2) = 0.84950784947877133746;
    nids.push_back(-1040);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.023457011864705432092;
    tri3_xyze(1, 0) = 0.087542760072790937187;
    tri3_xyze(2, 0) = 0.84911436253628169712;
    nids.push_back(7369);
    tri3_xyze(0, 1) = 0.021882930340033298483;
    tri3_xyze(1, 1) = 0.081668207847459012294;
    tri3_xyze(2, 1) = 0.84755282581460511526;
    nids.push_back(7371);
    tri3_xyze(0, 2) = 0.033232496120185459065;
    tri3_xyze(1, 2) = 0.080230342844840457639;
    tri3_xyze(2, 2) = 0.84833359417544351722;
    nids.push_back(-1041);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.021882930340033298483;
    tri3_xyze(1, 0) = 0.081668207847459012294;
    tri3_xyze(2, 0) = 0.84755282581460511526;
    nids.push_back(7371);
    tri3_xyze(0, 1) = 0.042274575140636268256;
    tri3_xyze(1, 1) = 0.073221712011954515065;
    tri3_xyze(2, 1) = 0.84755282581460544833;
    nids.push_back(7621);
    tri3_xyze(0, 2) = 0.033232496120185459065;
    tri3_xyze(1, 2) = 0.080230342844840457639;
    tri3_xyze(2, 2) = 0.84833359417544351722;
    nids.push_back(-1041);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042274575140636268256;
    tri3_xyze(1, 0) = 0.073221712011954515065;
    tri3_xyze(2, 0) = 0.84755282581460544833;
    nids.push_back(7621);
    tri3_xyze(0, 1) = 0.045315467135366820084;
    tri3_xyze(1, 1) = 0.078488691447157338255;
    tri3_xyze(2, 1) = 0.84911436253628203019;
    nids.push_back(7619);
    tri3_xyze(0, 2) = 0.033232496120185459065;
    tri3_xyze(1, 2) = 0.080230342844840457639;
    tri3_xyze(2, 2) = 0.84833359417544351722;
    nids.push_back(-1041);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.045315467135366820084;
    tri3_xyze(1, 0) = 0.078488691447157338255;
    tri3_xyze(2, 0) = 0.84911436253628203019;
    nids.push_back(7619);
    tri3_xyze(0, 1) = 0.023457011864705432092;
    tri3_xyze(1, 1) = 0.087542760072790937187;
    tri3_xyze(2, 1) = 0.84911436253628169712;
    nids.push_back(7369);
    tri3_xyze(0, 2) = 0.033232496120185459065;
    tri3_xyze(1, 2) = 0.080230342844840457639;
    tri3_xyze(2, 2) = 0.84833359417544351722;
    nids.push_back(-1041);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.021882930340033298483;
    tri3_xyze(1, 0) = 0.081668207847459012294;
    tri3_xyze(2, 0) = 0.84755282581460511526;
    nids.push_back(7371);
    tri3_xyze(0, 1) = 0.020371915026897051415;
    tri3_xyze(1, 1) = 0.076029021927819850912;
    tri3_xyze(2, 1) = 0.84524135262314858963;
    nids.push_back(7373);
    tri3_xyze(0, 2) = 0.030971234554612432632;
    tri3_xyze(1, 2) = 0.074771174505160992174;
    tri3_xyze(2, 2) = 0.84639708921887701898;
    nids.push_back(-1042);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.020371915026897051415;
    tri3_xyze(1, 0) = 0.076029021927819850912;
    tri3_xyze(2, 0) = 0.84524135262314858963;
    nids.push_back(7373);
    tri3_xyze(0, 1) = 0.039355517710883122784;
    tri3_xyze(1, 1) = 0.068165756233410548792;
    tri3_xyze(2, 1) = 0.8452413526231489227;
    nids.push_back(7623);
    tri3_xyze(0, 2) = 0.030971234554612432632;
    tri3_xyze(1, 2) = 0.074771174505160992174;
    tri3_xyze(2, 2) = 0.84639708921887701898;
    nids.push_back(-1042);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.039355517710883122784;
    tri3_xyze(1, 0) = 0.068165756233410548792;
    tri3_xyze(2, 0) = 0.8452413526231489227;
    nids.push_back(7623);
    tri3_xyze(0, 1) = 0.042274575140636268256;
    tri3_xyze(1, 1) = 0.073221712011954515065;
    tri3_xyze(2, 1) = 0.84755282581460544833;
    nids.push_back(7621);
    tri3_xyze(0, 2) = 0.030971234554612432632;
    tri3_xyze(1, 2) = 0.074771174505160992174;
    tri3_xyze(2, 2) = 0.84639708921887701898;
    nids.push_back(-1042);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042274575140636268256;
    tri3_xyze(1, 0) = 0.073221712011954515065;
    tri3_xyze(2, 0) = 0.84755282581460544833;
    nids.push_back(7621);
    tri3_xyze(0, 1) = 0.021882930340033298483;
    tri3_xyze(1, 1) = 0.081668207847459012294;
    tri3_xyze(2, 1) = 0.84755282581460511526;
    nids.push_back(7371);
    tri3_xyze(0, 2) = 0.030971234554612432632;
    tri3_xyze(1, 2) = 0.074771174505160992174;
    tri3_xyze(2, 2) = 0.84639708921887701898;
    nids.push_back(-1042);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.020371915026897051415;
    tri3_xyze(1, 0) = 0.076029021927819850912;
    tri3_xyze(2, 0) = 0.84524135262314858963;
    nids.push_back(7373);
    tri3_xyze(0, 1) = 0.018947795539421774524;
    tri3_xyze(1, 1) = 0.070714135644512593482;
    tri3_xyze(2, 1) = 0.8422163962749485977;
    nids.push_back(7375);
    tri3_xyze(0, 2) = 0.028819889600684246711;
    tri3_xyze(1, 2) = 0.069577368340044451478;
    tri3_xyze(2, 2) = 0.8437288744490487602;
    nids.push_back(-1043);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.018947795539421774524;
    tri3_xyze(1, 0) = 0.070714135644512593482;
    tri3_xyze(2, 0) = 0.8422163962749485977;
    nids.push_back(7375);
    tri3_xyze(0, 1) = 0.036604330125535038121;
    tri3_xyze(1, 1) = 0.063400559554434840481;
    tri3_xyze(2, 1) = 0.84221639627494881974;
    nids.push_back(7625);
    tri3_xyze(0, 2) = 0.028819889600684246711;
    tri3_xyze(1, 2) = 0.069577368340044451478;
    tri3_xyze(2, 2) = 0.8437288744490487602;
    nids.push_back(-1043);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036604330125535038121;
    tri3_xyze(1, 0) = 0.063400559554434840481;
    tri3_xyze(2, 0) = 0.84221639627494881974;
    nids.push_back(7625);
    tri3_xyze(0, 1) = 0.039355517710883122784;
    tri3_xyze(1, 1) = 0.068165756233410548792;
    tri3_xyze(2, 1) = 0.8452413526231489227;
    nids.push_back(7623);
    tri3_xyze(0, 2) = 0.028819889600684246711;
    tri3_xyze(1, 2) = 0.069577368340044451478;
    tri3_xyze(2, 2) = 0.8437288744490487602;
    nids.push_back(-1043);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.039355517710883122784;
    tri3_xyze(1, 0) = 0.068165756233410548792;
    tri3_xyze(2, 0) = 0.8452413526231489227;
    nids.push_back(7623);
    tri3_xyze(0, 1) = 0.020371915026897051415;
    tri3_xyze(1, 1) = 0.076029021927819850912;
    tri3_xyze(2, 1) = 0.84524135262314858963;
    nids.push_back(7373);
    tri3_xyze(0, 2) = 0.028819889600684246711;
    tri3_xyze(1, 2) = 0.069577368340044451478;
    tri3_xyze(2, 2) = 0.8437288744490487602;
    nids.push_back(-1043);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.018947795539421774524;
    tri3_xyze(1, 0) = 0.070714135644512593482;
    tri3_xyze(2, 0) = 0.8422163962749485977;
    nids.push_back(7375);
    tri3_xyze(0, 1) = 0.017633031092652327182;
    tri3_xyze(1, 1) = 0.065807367929184193422;
    tri3_xyze(2, 1) = 0.83852566213863743982;
    nids.push_back(7377);
    tri3_xyze(0, 2) = 0.026812389253475461814;
    tri3_xyze(1, 2) = 0.064730833775344248648;
    tri3_xyze(2, 2) = 0.84037102920679307427;
    nids.push_back(-1044);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.017633031092652327182;
    tri3_xyze(1, 0) = 0.065807367929184193422;
    tri3_xyze(2, 0) = 0.83852566213863743982;
    nids.push_back(7377);
    tri3_xyze(0, 1) = 0.034064400256292717839;
    tri3_xyze(1, 1) = 0.059001271973245346392;
    tri3_xyze(2, 1) = 0.83852566213863766187;
    nids.push_back(7627);
    tri3_xyze(0, 2) = 0.026812389253475461814;
    tri3_xyze(1, 2) = 0.064730833775344248648;
    tri3_xyze(2, 2) = 0.84037102920679307427;
    nids.push_back(-1044);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.034064400256292717839;
    tri3_xyze(1, 0) = 0.059001271973245346392;
    tri3_xyze(2, 0) = 0.83852566213863766187;
    nids.push_back(7627);
    tri3_xyze(0, 1) = 0.036604330125535038121;
    tri3_xyze(1, 1) = 0.063400559554434840481;
    tri3_xyze(2, 1) = 0.84221639627494881974;
    nids.push_back(7625);
    tri3_xyze(0, 2) = 0.026812389253475461814;
    tri3_xyze(1, 2) = 0.064730833775344248648;
    tri3_xyze(2, 2) = 0.84037102920679307427;
    nids.push_back(-1044);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036604330125535038121;
    tri3_xyze(1, 0) = 0.063400559554434840481;
    tri3_xyze(2, 0) = 0.84221639627494881974;
    nids.push_back(7625);
    tri3_xyze(0, 1) = 0.018947795539421774524;
    tri3_xyze(1, 1) = 0.070714135644512593482;
    tri3_xyze(2, 1) = 0.8422163962749485977;
    nids.push_back(7375);
    tri3_xyze(0, 2) = 0.026812389253475461814;
    tri3_xyze(1, 2) = 0.064730833775344248648;
    tri3_xyze(2, 2) = 0.84037102920679307427;
    nids.push_back(-1044);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.017633031092652327182;
    tri3_xyze(1, 0) = 0.065807367929184193422;
    tri3_xyze(2, 0) = 0.83852566213863743982;
    nids.push_back(7377);
    tri3_xyze(0, 1) = 0.016448356307316801261;
    tri3_xyze(1, 1) = 0.061386101439866114227;
    tri3_xyze(2, 1) = 0.83422735529628255779;
    nids.push_back(7379);
    tri3_xyze(0, 2) = 0.02498039299268412991;
    tri3_xyze(1, 2) = 0.060308003556324943362;
    tri3_xyze(2, 2) = 0.83637650871746016534;
    nids.push_back(-1045);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.016448356307316801261;
    tri3_xyze(1, 0) = 0.061386101439866114227;
    tri3_xyze(2, 0) = 0.83422735529628255779;
    nids.push_back(7379);
    tri3_xyze(0, 1) = 0.031775784314474669889;
    tri3_xyze(1, 1) = 0.055037272883004147161;
    tri3_xyze(2, 1) = 0.83422735529628266882;
    nids.push_back(7629);
    tri3_xyze(0, 2) = 0.02498039299268412991;
    tri3_xyze(1, 2) = 0.060308003556324943362;
    tri3_xyze(2, 2) = 0.83637650871746016534;
    nids.push_back(-1045);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.031775784314474669889;
    tri3_xyze(1, 0) = 0.055037272883004147161;
    tri3_xyze(2, 0) = 0.83422735529628266882;
    nids.push_back(7629);
    tri3_xyze(0, 1) = 0.034064400256292717839;
    tri3_xyze(1, 1) = 0.059001271973245346392;
    tri3_xyze(2, 1) = 0.83852566213863766187;
    nids.push_back(7627);
    tri3_xyze(0, 2) = 0.02498039299268412991;
    tri3_xyze(1, 2) = 0.060308003556324943362;
    tri3_xyze(2, 2) = 0.83637650871746016534;
    nids.push_back(-1045);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.034064400256292717839;
    tri3_xyze(1, 0) = 0.059001271973245346392;
    tri3_xyze(2, 0) = 0.83852566213863766187;
    nids.push_back(7627);
    tri3_xyze(0, 1) = 0.017633031092652327182;
    tri3_xyze(1, 1) = 0.065807367929184193422;
    tri3_xyze(2, 1) = 0.83852566213863743982;
    nids.push_back(7377);
    tri3_xyze(0, 2) = 0.02498039299268412991;
    tri3_xyze(1, 2) = 0.060308003556324943362;
    tri3_xyze(2, 2) = 0.83637650871746016534;
    nids.push_back(-1045);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.016448356307316801261;
    tri3_xyze(1, 0) = 0.061386101439866114227;
    tri3_xyze(2, 0) = 0.83422735529628255779;
    nids.push_back(7379);
    tri3_xyze(0, 1) = 0.015412454212470304282;
    tri3_xyze(1, 1) = 0.057520062190231792398;
    tri3_xyze(2, 1) = 0.82938926261447198929;
    nids.push_back(7381);
    tri3_xyze(0, 2) = 0.02335279249372451682;
    tri3_xyze(1, 2) = 0.056378628357611260125;
    tri3_xyze(2, 2) = 0.83180830895537738456;
    nids.push_back(-1046);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.015412454212470304282;
    tri3_xyze(1, 0) = 0.057520062190231792398;
    tri3_xyze(2, 0) = 0.82938926261447198929;
    nids.push_back(7381);
    tri3_xyze(0, 1) = 0.029774575140636295317;
    tri3_xyze(1, 1) = 0.051571076917342986712;
    tri3_xyze(2, 1) = 0.82938926261447210031;
    nids.push_back(7631);
    tri3_xyze(0, 2) = 0.02335279249372451682;
    tri3_xyze(1, 2) = 0.056378628357611260125;
    tri3_xyze(2, 2) = 0.83180830895537738456;
    nids.push_back(-1046);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.029774575140636295317;
    tri3_xyze(1, 0) = 0.051571076917342986712;
    tri3_xyze(2, 0) = 0.82938926261447210031;
    nids.push_back(7631);
    tri3_xyze(0, 1) = 0.031775784314474669889;
    tri3_xyze(1, 1) = 0.055037272883004147161;
    tri3_xyze(2, 1) = 0.83422735529628266882;
    nids.push_back(7629);
    tri3_xyze(0, 2) = 0.02335279249372451682;
    tri3_xyze(1, 2) = 0.056378628357611260125;
    tri3_xyze(2, 2) = 0.83180830895537738456;
    nids.push_back(-1046);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.031775784314474669889;
    tri3_xyze(1, 0) = 0.055037272883004147161;
    tri3_xyze(2, 0) = 0.83422735529628266882;
    nids.push_back(7629);
    tri3_xyze(0, 1) = 0.016448356307316801261;
    tri3_xyze(1, 1) = 0.061386101439866114227;
    tri3_xyze(2, 1) = 0.83422735529628255779;
    nids.push_back(7379);
    tri3_xyze(0, 2) = 0.02335279249372451682;
    tri3_xyze(1, 2) = 0.056378628357611260125;
    tri3_xyze(2, 2) = 0.83180830895537738456;
    nids.push_back(-1046);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.015412454212470304282;
    tri3_xyze(1, 0) = 0.057520062190231792398;
    tri3_xyze(2, 0) = 0.82938926261447198929;
    nids.push_back(7381);
    tri3_xyze(0, 1) = 0.014541661602966418979;
    tri3_xyze(1, 1) = 0.054270219928707640056;
    tri3_xyze(2, 1) = 0.82408768370493412636;
    nids.push_back(7383);
    tri3_xyze(0, 2) = 0.021955255988746612272;
    tri3_xyze(1, 2) = 0.053004676773381914634;
    tri3_xyze(2, 2) = 0.82673847315970316885;
    nids.push_back(-1047);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.014541661602966418979;
    tri3_xyze(1, 0) = 0.054270219928707640056;
    tri3_xyze(2, 0) = 0.82408768370493412636;
    nids.push_back(7383);
    tri3_xyze(0, 1) = 0.028092332998913428777;
    tri3_xyze(1, 1) = 0.048657348057245225492;
    tri3_xyze(2, 1) = 0.82408768370493423738;
    nids.push_back(7633);
    tri3_xyze(0, 2) = 0.021955255988746612272;
    tri3_xyze(1, 2) = 0.053004676773381914634;
    tri3_xyze(2, 2) = 0.82673847315970316885;
    nids.push_back(-1047);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.028092332998913428777;
    tri3_xyze(1, 0) = 0.048657348057245225492;
    tri3_xyze(2, 0) = 0.82408768370493423738;
    nids.push_back(7633);
    tri3_xyze(0, 1) = 0.029774575140636295317;
    tri3_xyze(1, 1) = 0.051571076917342986712;
    tri3_xyze(2, 1) = 0.82938926261447210031;
    nids.push_back(7631);
    tri3_xyze(0, 2) = 0.021955255988746612272;
    tri3_xyze(1, 2) = 0.053004676773381914634;
    tri3_xyze(2, 2) = 0.82673847315970316885;
    nids.push_back(-1047);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.029774575140636295317;
    tri3_xyze(1, 0) = 0.051571076917342986712;
    tri3_xyze(2, 0) = 0.82938926261447210031;
    nids.push_back(7631);
    tri3_xyze(0, 1) = 0.015412454212470304282;
    tri3_xyze(1, 1) = 0.057520062190231792398;
    tri3_xyze(2, 1) = 0.82938926261447198929;
    nids.push_back(7381);
    tri3_xyze(0, 2) = 0.021955255988746612272;
    tri3_xyze(1, 2) = 0.053004676773381914634;
    tri3_xyze(2, 2) = 0.82673847315970316885;
    nids.push_back(-1047);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.014541661602966418979;
    tri3_xyze(1, 0) = 0.054270219928707640056;
    tri3_xyze(2, 0) = 0.82408768370493412636;
    nids.push_back(7383);
    tri3_xyze(0, 1) = 0.013849711398443341481;
    tri3_xyze(1, 1) = 0.051687826609119649768;
    tri3_xyze(2, 1) = 0.81840622763408232831;
    nids.push_back(7385);
    tri3_xyze(0, 2) = 0.020809823463281716671;
    tri3_xyze(1, 2) = 0.050239358035621121956;
    tri3_xyze(2, 2) = 0.82124695566950833836;
    nids.push_back(-1048);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013849711398443341481;
    tri3_xyze(1, 0) = 0.051687826609119649768;
    tri3_xyze(2, 0) = 0.81840622763408232831;
    nids.push_back(7385);
    tri3_xyze(0, 1) = 0.026755587852803689591;
    tri3_xyze(1, 1) = 0.046342037547411993326;
    tri3_xyze(2, 1) = 0.81840622763408255036;
    nids.push_back(7635);
    tri3_xyze(0, 2) = 0.020809823463281716671;
    tri3_xyze(1, 2) = 0.050239358035621121956;
    tri3_xyze(2, 2) = 0.82124695566950833836;
    nids.push_back(-1048);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026755587852803689591;
    tri3_xyze(1, 0) = 0.046342037547411993326;
    tri3_xyze(2, 0) = 0.81840622763408255036;
    nids.push_back(7635);
    tri3_xyze(0, 1) = 0.028092332998913428777;
    tri3_xyze(1, 1) = 0.048657348057245225492;
    tri3_xyze(2, 1) = 0.82408768370493423738;
    nids.push_back(7633);
    tri3_xyze(0, 2) = 0.020809823463281716671;
    tri3_xyze(1, 2) = 0.050239358035621121956;
    tri3_xyze(2, 2) = 0.82124695566950833836;
    nids.push_back(-1048);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.028092332998913428777;
    tri3_xyze(1, 0) = 0.048657348057245225492;
    tri3_xyze(2, 0) = 0.82408768370493423738;
    nids.push_back(7633);
    tri3_xyze(0, 1) = 0.014541661602966418979;
    tri3_xyze(1, 1) = 0.054270219928707640056;
    tri3_xyze(2, 1) = 0.82408768370493412636;
    nids.push_back(7383);
    tri3_xyze(0, 2) = 0.020809823463281716671;
    tri3_xyze(1, 2) = 0.050239358035621121956;
    tri3_xyze(2, 2) = 0.82124695566950833836;
    nids.push_back(-1048);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013849711398443341481;
    tri3_xyze(1, 0) = 0.051687826609119649768;
    tri3_xyze(2, 0) = 0.81840622763408232831;
    nids.push_back(7385);
    tri3_xyze(0, 1) = 0.013347516066977403032;
    tri3_xyze(1, 1) = 0.04981360811676469047;
    tri3_xyze(2, 1) = 0.81243449435809123038;
    nids.push_back(7387);
    tri3_xyze(0, 2) = 0.019934559072504659527;
    tri3_xyze(1, 2) = 0.048126282872744771679;
    tri3_xyze(2, 2) = 0.81542036099608694588;
    nids.push_back(-1049);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013347516066977403032;
    tri3_xyze(1, 0) = 0.04981360811676469047;
    tri3_xyze(2, 0) = 0.81243449435809123038;
    nids.push_back(7387);
    tri3_xyze(0, 1) = 0.025785420971794200534;
    tri3_xyze(1, 1) = 0.044661659217682739276;
    tri3_xyze(2, 1) = 0.8124344943580913414;
    nids.push_back(7637);
    tri3_xyze(0, 2) = 0.019934559072504659527;
    tri3_xyze(1, 2) = 0.048126282872744771679;
    tri3_xyze(2, 2) = 0.81542036099608694588;
    nids.push_back(-1049);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025785420971794200534;
    tri3_xyze(1, 0) = 0.044661659217682739276;
    tri3_xyze(2, 0) = 0.8124344943580913414;
    nids.push_back(7637);
    tri3_xyze(0, 1) = 0.026755587852803689591;
    tri3_xyze(1, 1) = 0.046342037547411993326;
    tri3_xyze(2, 1) = 0.81840622763408255036;
    nids.push_back(7635);
    tri3_xyze(0, 2) = 0.019934559072504659527;
    tri3_xyze(1, 2) = 0.048126282872744771679;
    tri3_xyze(2, 2) = 0.81542036099608694588;
    nids.push_back(-1049);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026755587852803689591;
    tri3_xyze(1, 0) = 0.046342037547411993326;
    tri3_xyze(2, 0) = 0.81840622763408255036;
    nids.push_back(7635);
    tri3_xyze(0, 1) = 0.013849711398443341481;
    tri3_xyze(1, 1) = 0.051687826609119649768;
    tri3_xyze(2, 1) = 0.81840622763408232831;
    nids.push_back(7385);
    tri3_xyze(0, 2) = 0.019934559072504659527;
    tri3_xyze(1, 2) = 0.048126282872744771679;
    tri3_xyze(2, 2) = 0.81542036099608694588;
    nids.push_back(-1049);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013347516066977403032;
    tri3_xyze(1, 0) = 0.04981360811676469047;
    tri3_xyze(2, 0) = 0.81243449435809123038;
    nids.push_back(7387);
    tri3_xyze(0, 1) = 0.013042995528942823505;
    tri3_xyze(1, 1) = 0.048677121996871054732;
    tri3_xyze(2, 1) = 0.80626666167806371188;
    nids.push_back(7389);
    tri3_xyze(0, 2) = 0.019343266258715605815;
    tri3_xyze(1, 2) = 0.046698775742361281482;
    tri3_xyze(2, 2) = 0.80935057801807752664;
    nids.push_back(-1050);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.013042995528942823505;
    tri3_xyze(1, 0) = 0.048677121996871054732;
    tri3_xyze(2, 0) = 0.80626666167806371188;
    nids.push_back(7389);
    tri3_xyze(0, 1) = 0.025197132467147994456;
    tri3_xyze(1, 1) = 0.043642713638126676146;
    tri3_xyze(2, 1) = 0.80626666167806393393;
    nids.push_back(7639);
    tri3_xyze(0, 2) = 0.019343266258715605815;
    tri3_xyze(1, 2) = 0.046698775742361281482;
    tri3_xyze(2, 2) = 0.80935057801807752664;
    nids.push_back(-1050);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025197132467147994456;
    tri3_xyze(1, 0) = 0.043642713638126676146;
    tri3_xyze(2, 0) = 0.80626666167806393393;
    nids.push_back(7639);
    tri3_xyze(0, 1) = 0.025785420971794200534;
    tri3_xyze(1, 1) = 0.044661659217682739276;
    tri3_xyze(2, 1) = 0.8124344943580913414;
    nids.push_back(7637);
    tri3_xyze(0, 2) = 0.019343266258715605815;
    tri3_xyze(1, 2) = 0.046698775742361281482;
    tri3_xyze(2, 2) = 0.80935057801807752664;
    nids.push_back(-1050);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025785420971794200534;
    tri3_xyze(1, 0) = 0.044661659217682739276;
    tri3_xyze(2, 0) = 0.8124344943580913414;
    nids.push_back(7637);
    tri3_xyze(0, 1) = 0.013347516066977403032;
    tri3_xyze(1, 1) = 0.04981360811676469047;
    tri3_xyze(2, 1) = 0.81243449435809123038;
    nids.push_back(7387);
    tri3_xyze(0, 2) = 0.019343266258715605815;
    tri3_xyze(1, 2) = 0.046698775742361281482;
    tri3_xyze(2, 2) = 0.80935057801807752664;
    nids.push_back(-1050);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025197132467147994456;
    tri3_xyze(1, 0) = 0.043642713638126676146;
    tri3_xyze(2, 0) = 0.80626666167806393393;
    nids.push_back(7639);
    tri3_xyze(0, 1) = 0.025000000000009965639;
    tri3_xyze(1, 1) = 0.043301270189221953422;
    tri3_xyze(2, 1) = 0.79999999999984872101;
    nids.push_back(7541);
    tri3_xyze(0, 2) = 0.030296649498610835061;
    tri3_xyze(1, 2) = 0.039483362338653592294;
    tri3_xyze(2, 2) = 0.80313333083895643849;
    nids.push_back(-1052);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025000000000009965639;
    tri3_xyze(1, 0) = 0.043301270189221953422;
    tri3_xyze(2, 0) = 0.79999999999984872101;
    nids.push_back(7541);
    tri3_xyze(0, 1) = 0.035355339059337340113;
    tri3_xyze(1, 1) = 0.035355339059327417495;
    tri3_xyze(2, 1) = 0.79999999999984894306;
    nids.push_back(7791);
    tri3_xyze(0, 2) = 0.030296649498610835061;
    tri3_xyze(1, 2) = 0.039483362338653592294;
    tri3_xyze(2, 2) = 0.80313333083895643849;
    nids.push_back(-1052);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.035355339059337340113;
    tri3_xyze(1, 0) = 0.035355339059327417495;
    tri3_xyze(2, 0) = 0.79999999999984894306;
    nids.push_back(7791);
    tri3_xyze(0, 1) = 0.035634126467948050443;
    tri3_xyze(1, 1) = 0.035634126467938308236;
    tri3_xyze(2, 1) = 0.80626666167806415597;
    nids.push_back(7889);
    tri3_xyze(0, 2) = 0.030296649498610835061;
    tri3_xyze(1, 2) = 0.039483362338653592294;
    tri3_xyze(2, 2) = 0.80313333083895643849;
    nids.push_back(-1052);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.035634126467948050443;
    tri3_xyze(1, 0) = 0.035634126467938308236;
    tri3_xyze(2, 0) = 0.80626666167806415597;
    nids.push_back(7889);
    tri3_xyze(0, 1) = 0.025197132467147994456;
    tri3_xyze(1, 1) = 0.043642713638126676146;
    tri3_xyze(2, 1) = 0.80626666167806393393;
    nids.push_back(7639);
    tri3_xyze(0, 2) = 0.030296649498610835061;
    tri3_xyze(1, 2) = 0.039483362338653592294;
    tri3_xyze(2, 2) = 0.80313333083895643849;
    nids.push_back(-1052);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.028092332998913432246;
    tri3_xyze(1, 0) = 0.048657348057243816897;
    tri3_xyze(2, 0) = 0.7759123162947625385;
    nids.push_back(7549);
    tri3_xyze(0, 1) = 0.029774575140636246745;
    tri3_xyze(1, 1) = 0.051571076917341147905;
    tri3_xyze(2, 1) = 0.77061073738522467558;
    nids.push_back(7551);
    tri3_xyze(0, 2) = 0.03492576861077266448;
    tri3_xyze(1, 2) = 0.04551614781952607186;
    tri3_xyze(2, 2) = 0.77326152683999382909;
    nids.push_back(-1056);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.029774575140636246745;
    tri3_xyze(1, 0) = 0.051571076917341147905;
    tri3_xyze(2, 0) = 0.77061073738522467558;
    nids.push_back(7551);
    tri3_xyze(0, 1) = 0.042107607977780542319;
    tri3_xyze(1, 1) = 0.0421076079777696205;
    tri3_xyze(2, 1) = 0.77061073738522500864;
    nids.push_back(7801);
    tri3_xyze(0, 2) = 0.03492576861077266448;
    tri3_xyze(1, 2) = 0.04551614781952607186;
    tri3_xyze(2, 2) = 0.77326152683999382909;
    nids.push_back(-1056);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.029774575140636246745;
    tri3_xyze(1, 0) = 0.051571076917341147905;
    tri3_xyze(2, 0) = 0.77061073738522467558;
    nids.push_back(7551);
    tri3_xyze(0, 1) = 0.031775784314474669889;
    tri3_xyze(1, 1) = 0.055037272883001996104;
    tri3_xyze(2, 1) = 0.76577264470341388503;
    nids.push_back(7553);
    tri3_xyze(0, 2) = 0.037148928141364938316;
    tri3_xyze(1, 2) = 0.048413425727667491749;
    tri3_xyze(2, 2) = 0.76819169104431939132;
    nids.push_back(-1057);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.031775784314474669889;
    tri3_xyze(1, 0) = 0.055037272883001996104;
    tri3_xyze(2, 0) = 0.76577264470341388503;
    nids.push_back(7553);
    tri3_xyze(0, 1) = 0.044937745132568283901;
    tri3_xyze(1, 1) = 0.044937745132557209427;
    tri3_xyze(2, 1) = 0.76577264470341410707;
    nids.push_back(7803);
    tri3_xyze(0, 2) = 0.037148928141364938316;
    tri3_xyze(1, 2) = 0.048413425727667491749;
    tri3_xyze(2, 2) = 0.76819169104431939132;
    nids.push_back(-1057);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042107607977780542319;
    tri3_xyze(1, 0) = 0.0421076079777696205;
    tri3_xyze(2, 0) = 0.77061073738522500864;
    nids.push_back(7801);
    tri3_xyze(0, 1) = 0.029774575140636246745;
    tri3_xyze(1, 1) = 0.051571076917341147905;
    tri3_xyze(2, 1) = 0.77061073738522467558;
    nids.push_back(7551);
    tri3_xyze(0, 2) = 0.037148928141364938316;
    tri3_xyze(1, 2) = 0.048413425727667491749;
    tri3_xyze(2, 2) = 0.76819169104431939132;
    nids.push_back(-1057);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.031775784314474669889;
    tri3_xyze(1, 0) = 0.055037272883001996104;
    tri3_xyze(2, 0) = 0.76577264470341388503;
    nids.push_back(7553);
    tri3_xyze(0, 1) = 0.034064400256292683145;
    tri3_xyze(1, 1) = 0.059001271973242869207;
    tri3_xyze(2, 1) = 0.76147433786105866993;
    nids.push_back(7555);
    tri3_xyze(0, 2) = 0.039738066634971550006;
    tri3_xyze(1, 2) = 0.051787656706335366552;
    tri3_xyze(2, 2) = 0.7636234912822363885;
    nids.push_back(-1058);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.034064400256292683145;
    tri3_xyze(1, 0) = 0.059001271973242869207;
    tri3_xyze(2, 0) = 0.76147433786105866993;
    nids.push_back(7555);
    tri3_xyze(0, 1) = 0.048174336836550576968;
    tri3_xyze(1, 1) = 0.048174336836539370654;
    tri3_xyze(2, 1) = 0.76147433786105911402;
    nids.push_back(7805);
    tri3_xyze(0, 2) = 0.039738066634971550006;
    tri3_xyze(1, 2) = 0.051787656706335366552;
    tri3_xyze(2, 2) = 0.7636234912822363885;
    nids.push_back(-1058);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048174336836550576968;
    tri3_xyze(1, 0) = 0.048174336836539370654;
    tri3_xyze(2, 0) = 0.76147433786105911402;
    nids.push_back(7805);
    tri3_xyze(0, 1) = 0.044937745132568283901;
    tri3_xyze(1, 1) = 0.044937745132557209427;
    tri3_xyze(2, 1) = 0.76577264470341410707;
    nids.push_back(7803);
    tri3_xyze(0, 2) = 0.039738066634971550006;
    tri3_xyze(1, 2) = 0.051787656706335366552;
    tri3_xyze(2, 2) = 0.7636234912822363885;
    nids.push_back(-1058);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.044937745132568283901;
    tri3_xyze(1, 0) = 0.044937745132557209427;
    tri3_xyze(2, 0) = 0.76577264470341410707;
    nids.push_back(7803);
    tri3_xyze(0, 1) = 0.031775784314474669889;
    tri3_xyze(1, 1) = 0.055037272883001996104;
    tri3_xyze(2, 1) = 0.76577264470341388503;
    nids.push_back(7553);
    tri3_xyze(0, 2) = 0.039738066634971550006;
    tri3_xyze(1, 2) = 0.051787656706335366552;
    tri3_xyze(2, 2) = 0.7636234912822363885;
    nids.push_back(-1058);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.034064400256292683145;
    tri3_xyze(1, 0) = 0.059001271973242869207;
    tri3_xyze(2, 0) = 0.76147433786105866993;
    nids.push_back(7555);
    tri3_xyze(0, 1) = 0.036604330125535065876;
    tri3_xyze(1, 1) = 0.063400559554432342479;
    tri3_xyze(2, 1) = 0.75778360372474717899;
    nids.push_back(7557);
    tri3_xyze(0, 2) = 0.042652351830871987337;
    tri3_xyze(1, 2) = 0.055585627117328249847;
    tri3_xyze(2, 2) = 0.7596289707929031465;
    nids.push_back(-1059);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036604330125535065876;
    tri3_xyze(1, 0) = 0.063400559554432342479;
    tri3_xyze(2, 0) = 0.75778360372474717899;
    nids.push_back(7557);
    tri3_xyze(0, 1) = 0.051766340105109630298;
    tri3_xyze(1, 1) = 0.051766340105098403168;
    tri3_xyze(2, 1) = 0.75778360372474762308;
    nids.push_back(7807);
    tri3_xyze(0, 2) = 0.042652351830871987337;
    tri3_xyze(1, 2) = 0.055585627117328249847;
    tri3_xyze(2, 2) = 0.7596289707929031465;
    nids.push_back(-1059);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051766340105109630298;
    tri3_xyze(1, 0) = 0.051766340105098403168;
    tri3_xyze(2, 0) = 0.75778360372474762308;
    nids.push_back(7807);
    tri3_xyze(0, 1) = 0.048174336836550576968;
    tri3_xyze(1, 1) = 0.048174336836539370654;
    tri3_xyze(2, 1) = 0.76147433786105911402;
    nids.push_back(7805);
    tri3_xyze(0, 2) = 0.042652351830871987337;
    tri3_xyze(1, 2) = 0.055585627117328249847;
    tri3_xyze(2, 2) = 0.7596289707929031465;
    nids.push_back(-1059);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048174336836550576968;
    tri3_xyze(1, 0) = 0.048174336836539370654;
    tri3_xyze(2, 0) = 0.76147433786105911402;
    nids.push_back(7805);
    tri3_xyze(0, 1) = 0.034064400256292683145;
    tri3_xyze(1, 1) = 0.059001271973242869207;
    tri3_xyze(2, 1) = 0.76147433786105866993;
    nids.push_back(7555);
    tri3_xyze(0, 2) = 0.042652351830871987337;
    tri3_xyze(1, 2) = 0.055585627117328249847;
    tri3_xyze(2, 2) = 0.7596289707929031465;
    nids.push_back(-1059);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036604330125535065876;
    tri3_xyze(1, 0) = 0.063400559554432342479;
    tri3_xyze(2, 0) = 0.75778360372474717899;
    nids.push_back(7557);
    tri3_xyze(0, 1) = 0.039355517710883122784;
    tri3_xyze(1, 1) = 0.068165756233407717724;
    tri3_xyze(2, 1) = 0.75475864737654696501;
    nids.push_back(7559);
    tri3_xyze(0, 2) = 0.045845823710617301938;
    tri3_xyze(1, 2) = 0.059747440698467130238;
    tri3_xyze(2, 2) = 0.75627112555064734956;
    nids.push_back(-1060);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.039355517710883122784;
    tri3_xyze(1, 0) = 0.068165756233407717724;
    tri3_xyze(2, 0) = 0.75475864737654696501;
    nids.push_back(7559);
    tri3_xyze(0, 1) = 0.055657106900941374916;
    tri3_xyze(1, 1) = 0.05565710690093005758;
    tri3_xyze(2, 1) = 0.75475864737654729808;
    nids.push_back(7809);
    tri3_xyze(0, 2) = 0.045845823710617301938;
    tri3_xyze(1, 2) = 0.059747440698467130238;
    tri3_xyze(2, 2) = 0.75627112555064734956;
    nids.push_back(-1060);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.055657106900941374916;
    tri3_xyze(1, 0) = 0.05565710690093005758;
    tri3_xyze(2, 0) = 0.75475864737654729808;
    nids.push_back(7809);
    tri3_xyze(0, 1) = 0.051766340105109630298;
    tri3_xyze(1, 1) = 0.051766340105098403168;
    tri3_xyze(2, 1) = 0.75778360372474762308;
    nids.push_back(7807);
    tri3_xyze(0, 2) = 0.045845823710617301938;
    tri3_xyze(1, 2) = 0.059747440698467130238;
    tri3_xyze(2, 2) = 0.75627112555064734956;
    nids.push_back(-1060);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051766340105109630298;
    tri3_xyze(1, 0) = 0.051766340105098403168;
    tri3_xyze(2, 0) = 0.75778360372474762308;
    nids.push_back(7807);
    tri3_xyze(0, 1) = 0.036604330125535065876;
    tri3_xyze(1, 1) = 0.063400559554432342479;
    tri3_xyze(2, 1) = 0.75778360372474717899;
    nids.push_back(7557);
    tri3_xyze(0, 2) = 0.045845823710617301938;
    tri3_xyze(1, 2) = 0.059747440698467130238;
    tri3_xyze(2, 2) = 0.75627112555064734956;
    nids.push_back(-1060);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.039355517710883122784;
    tri3_xyze(1, 0) = 0.068165756233407717724;
    tri3_xyze(2, 0) = 0.75475864737654696501;
    nids.push_back(7559);
    tri3_xyze(0, 1) = 0.042274575140636413972;
    tri3_xyze(1, 1) = 0.073221712011951434196;
    tri3_xyze(2, 1) = 0.75244717418508999529;
    nids.push_back(7561);
    tri3_xyze(0, 2) = 0.049268119314976266399;
    tri3_xyze(1, 2) = 0.064207463163430497644;
    tri3_xyze(2, 2) = 0.75360291078081864669;
    nids.push_back(-1061);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042274575140636413972;
    tri3_xyze(1, 0) = 0.073221712011951434196;
    tri3_xyze(2, 0) = 0.75244717418508999529;
    nids.push_back(7561);
    tri3_xyze(0, 1) = 0.059785277507444167799;
    tri3_xyze(1, 1) = 0.059785277507432788013;
    tri3_xyze(2, 1) = 0.75244717418509043938;
    nids.push_back(7811);
    tri3_xyze(0, 2) = 0.049268119314976266399;
    tri3_xyze(1, 2) = 0.064207463163430497644;
    tri3_xyze(2, 2) = 0.75360291078081864669;
    nids.push_back(-1061);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.059785277507444167799;
    tri3_xyze(1, 0) = 0.059785277507432788013;
    tri3_xyze(2, 0) = 0.75244717418509043938;
    nids.push_back(7811);
    tri3_xyze(0, 1) = 0.055657106900941374916;
    tri3_xyze(1, 1) = 0.05565710690093005758;
    tri3_xyze(2, 1) = 0.75475864737654729808;
    nids.push_back(7809);
    tri3_xyze(0, 2) = 0.049268119314976266399;
    tri3_xyze(1, 2) = 0.064207463163430497644;
    tri3_xyze(2, 2) = 0.75360291078081864669;
    nids.push_back(-1061);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.055657106900941374916;
    tri3_xyze(1, 0) = 0.05565710690093005758;
    tri3_xyze(2, 0) = 0.75475864737654729808;
    nids.push_back(7809);
    tri3_xyze(0, 1) = 0.039355517710883122784;
    tri3_xyze(1, 1) = 0.068165756233407717724;
    tri3_xyze(2, 1) = 0.75475864737654696501;
    nids.push_back(7559);
    tri3_xyze(0, 2) = 0.049268119314976266399;
    tri3_xyze(1, 2) = 0.064207463163430497644;
    tri3_xyze(2, 2) = 0.75360291078081864669;
    nids.push_back(-1061);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042274575140636413972;
    tri3_xyze(1, 0) = 0.073221712011951434196;
    tri3_xyze(2, 0) = 0.75244717418508999529;
    nids.push_back(7561);
    tri3_xyze(0, 1) = 0.045315467135366757634;
    tri3_xyze(1, 1) = 0.078488691447154299019;
    tri3_xyze(2, 1) = 0.75088563746341319138;
    nids.push_back(7563);
    tri3_xyze(0, 2) = 0.052865266997887808142;
    tri3_xyze(1, 2) = 0.068895357293657755182;
    tri3_xyze(2, 2) = 0.75166640582425181538;
    nids.push_back(-1062);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.045315467135366757634;
    tri3_xyze(1, 0) = 0.078488691447154299019;
    tri3_xyze(2, 0) = 0.75088563746341319138;
    nids.push_back(7563);
    tri3_xyze(0, 1) = 0.064085748208103893164;
    tri3_xyze(1, 1) = 0.064085748208092471745;
    tri3_xyze(2, 1) = 0.75088563746341363547;
    nids.push_back(7813);
    tri3_xyze(0, 2) = 0.052865266997887808142;
    tri3_xyze(1, 2) = 0.068895357293657755182;
    tri3_xyze(2, 2) = 0.75166640582425181538;
    nids.push_back(-1062);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.059785277507444167799;
    tri3_xyze(1, 0) = 0.059785277507432788013;
    tri3_xyze(2, 0) = 0.75244717418509043938;
    nids.push_back(7811);
    tri3_xyze(0, 1) = 0.042274575140636413972;
    tri3_xyze(1, 1) = 0.073221712011951434196;
    tri3_xyze(2, 1) = 0.75244717418508999529;
    nids.push_back(7561);
    tri3_xyze(0, 2) = 0.052865266997887808142;
    tri3_xyze(1, 2) = 0.068895357293657755182;
    tri3_xyze(2, 2) = 0.75166640582425181538;
    nids.push_back(-1062);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.045315467135366757634;
    tri3_xyze(1, 0) = 0.078488691447154299019;
    tri3_xyze(2, 0) = 0.75088563746341319138;
    nids.push_back(7563);
    tri3_xyze(0, 1) = 0.048430237011777121925;
    tri3_xyze(1, 1) = 0.083883631126981925141;
    tri3_xyze(2, 1) = 0.75009866357843379969;
    nids.push_back(7565);
    tri3_xyze(0, 2) = 0.056580537591560599309;
    tri3_xyze(1, 2) = 0.073737192198302936585;
    tri3_xyze(2, 2) = 0.75049215052092366207;
    nids.push_back(-1063);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048430237011777121925;
    tri3_xyze(1, 0) = 0.083883631126981925141;
    tri3_xyze(2, 0) = 0.75009866357843379969;
    nids.push_back(7565);
    tri3_xyze(0, 1) = 0.068490698010994610634;
    tri3_xyze(1, 1) = 0.068490698010983036559;
    tri3_xyze(2, 1) = 0.7500986635784343548;
    nids.push_back(7815);
    tri3_xyze(0, 2) = 0.056580537591560599309;
    tri3_xyze(1, 2) = 0.073737192198302936585;
    tri3_xyze(2, 2) = 0.75049215052092366207;
    nids.push_back(-1063);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.064085748208103893164;
    tri3_xyze(1, 0) = 0.064085748208092471745;
    tri3_xyze(2, 0) = 0.75088563746341363547;
    nids.push_back(7813);
    tri3_xyze(0, 1) = 0.045315467135366757634;
    tri3_xyze(1, 1) = 0.078488691447154299019;
    tri3_xyze(2, 1) = 0.75088563746341319138;
    nids.push_back(7563);
    tri3_xyze(0, 2) = 0.056580537591560599309;
    tri3_xyze(1, 2) = 0.073737192198302936585;
    tri3_xyze(2, 2) = 0.75049215052092366207;
    nids.push_back(-1063);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048430237011777121925;
    tri3_xyze(1, 0) = 0.083883631126981925141;
    tri3_xyze(2, 0) = 0.75009866357843379969;
    nids.push_back(7565);
    tri3_xyze(0, 1) = 0.051569762988242798252;
    tri3_xyze(1, 1) = 0.089321449629902766043;
    tri3_xyze(2, 1) = 0.75009866357843368867;
    nids.push_back(7567);
    tri3_xyze(0, 2) = 0.060355339059337327623;
    tri3_xyze(1, 2) = 0.078656609248547754154;
    tri3_xyze(2, 2) = 0.75009866357843391071;
    nids.push_back(-1064);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.068490698010994610634;
    tri3_xyze(1, 0) = 0.068490698010983036559;
    tri3_xyze(2, 0) = 0.7500986635784343548;
    nids.push_back(7815);
    tri3_xyze(0, 1) = 0.048430237011777121925;
    tri3_xyze(1, 1) = 0.083883631126981925141;
    tri3_xyze(2, 1) = 0.75009866357843379969;
    nids.push_back(7565);
    tri3_xyze(0, 2) = 0.060355339059337327623;
    tri3_xyze(1, 2) = 0.078656609248547754154;
    tri3_xyze(2, 2) = 0.75009866357843391071;
    nids.push_back(-1064);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.070225424859383614451;
    tri3_xyze(1, 0) = 0.12163400383954667272;
    tri3_xyze(2, 0) = 0.82938926261446987986;
    nids.push_back(7601);
    tri3_xyze(0, 1) = 0.099313748259548931263;
    tri3_xyze(1, 1) = 0.09931374825953985519;
    tri3_xyze(2, 1) = 0.82938926261447065702;
    nids.push_back(7851);
    tri3_xyze(0, 2) = 0.085784909507902018522;
    tri3_xyze(1, 2) = 0.11179707067757263528;
    tri3_xyze(2, 2) = 0.82673847315970128147;
    nids.push_back(-1081);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.099313748259548931263;
    tri3_xyze(1, 0) = 0.09931374825953985519;
    tri3_xyze(2, 0) = 0.82938926261447065702;
    nids.push_back(7851);
    tri3_xyze(0, 1) = 0.10169279791156897452;
    tri3_xyze(1, 1) = 0.1016927979115599262;
    tri3_xyze(2, 1) = 0.82408768370493257205;
    nids.push_back(7849);
    tri3_xyze(0, 2) = 0.085784909507902018522;
    tri3_xyze(1, 2) = 0.11179707067757263528;
    tri3_xyze(2, 2) = 0.82673847315970128147;
    nids.push_back(-1081);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.068224215685545236409;
    tri3_xyze(1, 0) = 0.11816780787388576901;
    tri3_xyze(2, 0) = 0.83422735529628067042;
    nids.push_back(7603);
    tri3_xyze(0, 1) = 0.096483611104761266009;
    tri3_xyze(1, 1) = 0.096483611104752328713;
    tri3_xyze(2, 1) = 0.83422735529628144757;
    nids.push_back(7853);
    tri3_xyze(0, 2) = 0.083561749977309765502;
    tri3_xyze(1, 2) = 0.10889979276943116682;
    tri3_xyze(2, 2) = 0.83180830895537560821;
    nids.push_back(-1082);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.096483611104761266009;
    tri3_xyze(1, 0) = 0.096483611104752328713;
    tri3_xyze(2, 0) = 0.83422735529628144757;
    nids.push_back(7853);
    tri3_xyze(0, 1) = 0.099313748259548931263;
    tri3_xyze(1, 1) = 0.09931374825953985519;
    tri3_xyze(2, 1) = 0.82938926261447065702;
    nids.push_back(7851);
    tri3_xyze(0, 2) = 0.083561749977309765502;
    tri3_xyze(1, 2) = 0.10889979276943116682;
    tri3_xyze(2, 2) = 0.83180830895537560821;
    nids.push_back(-1082);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.099313748259548931263;
    tri3_xyze(1, 0) = 0.09931374825953985519;
    tri3_xyze(2, 0) = 0.82938926261447065702;
    nids.push_back(7851);
    tri3_xyze(0, 1) = 0.070225424859383614451;
    tri3_xyze(1, 1) = 0.12163400383954667272;
    tri3_xyze(2, 1) = 0.82938926261446987986;
    nids.push_back(7601);
    tri3_xyze(0, 2) = 0.083561749977309765502;
    tri3_xyze(1, 2) = 0.10889979276943116682;
    tri3_xyze(2, 2) = 0.83180830895537560821;
    nids.push_back(-1082);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.065935599743727049682;
    tri3_xyze(1, 0) = 0.11420380878364494448;
    tri3_xyze(2, 0) = 0.83852566213863599653;
    nids.push_back(7605);
    tri3_xyze(0, 1) = 0.093247019400778799469;
    tri3_xyze(1, 1) = 0.093247019400770181363;
    tri3_xyze(2, 1) = 0.83852566213863666267;
    nids.push_back(7855);
    tri3_xyze(0, 2) = 0.080972611483703091362;
    tri3_xyze(1, 2) = 0.10552556179076330589;
    tri3_xyze(2, 2) = 0.83637650871745872205;
    nids.push_back(-1083);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.093247019400778799469;
    tri3_xyze(1, 0) = 0.093247019400770181363;
    tri3_xyze(2, 0) = 0.83852566213863666267;
    nids.push_back(7855);
    tri3_xyze(0, 1) = 0.096483611104761266009;
    tri3_xyze(1, 1) = 0.096483611104752328713;
    tri3_xyze(2, 1) = 0.83422735529628144757;
    nids.push_back(7853);
    tri3_xyze(0, 2) = 0.080972611483703091362;
    tri3_xyze(1, 2) = 0.10552556179076330589;
    tri3_xyze(2, 2) = 0.83637650871745872205;
    nids.push_back(-1083);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.096483611104761266009;
    tri3_xyze(1, 0) = 0.096483611104752328713;
    tri3_xyze(2, 0) = 0.83422735529628144757;
    nids.push_back(7853);
    tri3_xyze(0, 1) = 0.068224215685545236409;
    tri3_xyze(1, 1) = 0.11816780787388576901;
    tri3_xyze(2, 1) = 0.83422735529628067042;
    nids.push_back(7603);
    tri3_xyze(0, 2) = 0.080972611483703091362;
    tri3_xyze(1, 2) = 0.10552556179076330589;
    tri3_xyze(2, 2) = 0.83637650871745872205;
    nids.push_back(-1083);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.063395669874484861239;
    tri3_xyze(1, 0) = 0.10980452120245559611;
    tri3_xyze(2, 0) = 0.84221639627494737645;
    nids.push_back(7607);
    tri3_xyze(0, 1) = 0.0896550161322196143;
    tri3_xyze(1, 1) = 0.08965501613221139865;
    tri3_xyze(2, 1) = 0.84221639627494804259;
    nids.push_back(7857);
    tri3_xyze(0, 2) = 0.078058326287802570764;
    tri3_xyze(1, 2) = 0.10172759137977052668;
    tri3_xyze(2, 2) = 0.84037102920679196405;
    nids.push_back(-1084);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0896550161322196143;
    tri3_xyze(1, 0) = 0.08965501613221139865;
    tri3_xyze(2, 0) = 0.84221639627494804259;
    nids.push_back(7857);
    tri3_xyze(0, 1) = 0.093247019400778799469;
    tri3_xyze(1, 1) = 0.093247019400770181363;
    tri3_xyze(2, 1) = 0.83852566213863666267;
    nids.push_back(7855);
    tri3_xyze(0, 2) = 0.078058326287802570764;
    tri3_xyze(1, 2) = 0.10172759137977052668;
    tri3_xyze(2, 2) = 0.84037102920679196405;
    nids.push_back(-1084);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.093247019400778799469;
    tri3_xyze(1, 0) = 0.093247019400770181363;
    tri3_xyze(2, 0) = 0.83852566213863666267;
    nids.push_back(7855);
    tri3_xyze(0, 1) = 0.065935599743727049682;
    tri3_xyze(1, 1) = 0.11420380878364494448;
    tri3_xyze(2, 1) = 0.83852566213863599653;
    nids.push_back(7605);
    tri3_xyze(0, 2) = 0.078058326287802570764;
    tri3_xyze(1, 2) = 0.10172759137977052668;
    tri3_xyze(2, 2) = 0.84037102920679196405;
    nids.push_back(-1084);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.085764249336388140299;
    tri3_xyze(1, 0) = 0.085764249336379549948;
    tri3_xyze(2, 0) = 0.84524135262314836758;
    nids.push_back(7859);
    tri3_xyze(0, 1) = 0.0896550161322196143;
    tri3_xyze(1, 1) = 0.08965501613221139865;
    tri3_xyze(2, 1) = 0.84221639627494804259;
    nids.push_back(7857);
    tri3_xyze(0, 2) = 0.074864854408057332491;
    tri3_xyze(1, 2) = 0.097565777798631625473;
    tri3_xyze(2, 2) = 0.84372887444904787202;
    nids.push_back(-1085);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0896550161322196143;
    tri3_xyze(1, 0) = 0.08965501613221139865;
    tri3_xyze(2, 0) = 0.84221639627494804259;
    nids.push_back(7857);
    tri3_xyze(0, 1) = 0.063395669874484861239;
    tri3_xyze(1, 1) = 0.10980452120245559611;
    tri3_xyze(2, 1) = 0.84221639627494737645;
    nids.push_back(7607);
    tri3_xyze(0, 2) = 0.074864854408057332491;
    tri3_xyze(1, 2) = 0.097565777798631625473;
    tri3_xyze(2, 2) = 0.84372887444904787202;
    nids.push_back(-1085);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.060644482289136734943;
    tri3_xyze(1, 0) = 0.10503932452347999882;
    tri3_xyze(2, 0) = 0.84524135262314770145;
    nids.push_back(7609);
    tri3_xyze(0, 1) = 0.057725424859383624165;
    tri3_xyze(1, 1) = 0.099983368744936254591;
    tri3_xyze(2, 1) = 0.84755282581460467117;
    nids.push_back(7611);
    tri3_xyze(0, 2) = 0.071442558803698444359;
    tri3_xyze(1, 2) = 0.093105755333668147045;
    tri3_xyze(2, 2) = 0.84639708921887635285;
    nids.push_back(-1086);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.057725424859383624165;
    tri3_xyze(1, 0) = 0.099983368744936254591;
    tri3_xyze(2, 0) = 0.84755282581460467117;
    nids.push_back(7611);
    tri3_xyze(0, 1) = 0.081636078729885278027;
    tri3_xyze(1, 1) = 0.081636078729876770943;
    tri3_xyze(2, 1) = 0.84755282581460522628;
    nids.push_back(7861);
    tri3_xyze(0, 2) = 0.071442558803698444359;
    tri3_xyze(1, 2) = 0.093105755333668147045;
    tri3_xyze(2, 2) = 0.84639708921887635285;
    nids.push_back(-1086);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.081636078729885278027;
    tri3_xyze(1, 0) = 0.081636078729876770943;
    tri3_xyze(2, 0) = 0.84755282581460522628;
    nids.push_back(7861);
    tri3_xyze(0, 1) = 0.085764249336388140299;
    tri3_xyze(1, 1) = 0.085764249336379549948;
    tri3_xyze(2, 1) = 0.84524135262314836758;
    nids.push_back(7859);
    tri3_xyze(0, 2) = 0.071442558803698444359;
    tri3_xyze(1, 2) = 0.093105755333668147045;
    tri3_xyze(2, 2) = 0.84639708921887635285;
    nids.push_back(-1086);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.085764249336388140299;
    tri3_xyze(1, 0) = 0.085764249336379549948;
    tri3_xyze(2, 0) = 0.84524135262314836758;
    nids.push_back(7859);
    tri3_xyze(0, 1) = 0.060644482289136734943;
    tri3_xyze(1, 1) = 0.10503932452347999882;
    tri3_xyze(2, 1) = 0.84524135262314770145;
    nids.push_back(7609);
    tri3_xyze(0, 2) = 0.071442558803698444359;
    tri3_xyze(1, 2) = 0.093105755333668147045;
    tri3_xyze(2, 2) = 0.84639708921887635285;
    nids.push_back(-1086);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.057725424859383624165;
    tri3_xyze(1, 0) = 0.099983368744936254591;
    tri3_xyze(2, 0) = 0.84755282581460467117;
    nids.push_back(7611);
    tri3_xyze(0, 1) = 0.054684532864653072337;
    tri3_xyze(1, 1) = 0.094716389309733514668;
    tri3_xyze(2, 1) = 0.84911436253628147508;
    nids.push_back(7613);
    tri3_xyze(0, 2) = 0.067845411120786860981;
    tri3_xyze(1, 2) = 0.088417861203440931139;
    tri3_xyze(2, 2) = 0.84833359417544329517;
    nids.push_back(-1087);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.054684532864653072337;
    tri3_xyze(1, 0) = 0.094716389309733514668;
    tri3_xyze(2, 0) = 0.84911436253628147508;
    nids.push_back(7613);
    tri3_xyze(0, 1) = 0.077335608029225483273;
    tri3_xyze(1, 1) = 0.077335608029217170478;
    tri3_xyze(2, 1) = 0.84911436253628203019;
    nids.push_back(7863);
    tri3_xyze(0, 2) = 0.067845411120786860981;
    tri3_xyze(1, 2) = 0.088417861203440931139;
    tri3_xyze(2, 2) = 0.84833359417544329517;
    nids.push_back(-1087);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.077335608029225483273;
    tri3_xyze(1, 0) = 0.077335608029217170478;
    tri3_xyze(2, 0) = 0.84911436253628203019;
    nids.push_back(7863);
    tri3_xyze(0, 1) = 0.081636078729885278027;
    tri3_xyze(1, 1) = 0.081636078729876770943;
    tri3_xyze(2, 1) = 0.84755282581460522628;
    nids.push_back(7861);
    tri3_xyze(0, 2) = 0.067845411120786860981;
    tri3_xyze(1, 2) = 0.088417861203440931139;
    tri3_xyze(2, 2) = 0.84833359417544329517;
    nids.push_back(-1087);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.081636078729885278027;
    tri3_xyze(1, 0) = 0.081636078729876770943;
    tri3_xyze(2, 0) = 0.84755282581460522628;
    nids.push_back(7861);
    tri3_xyze(0, 1) = 0.057725424859383624165;
    tri3_xyze(1, 1) = 0.099983368744936254591;
    tri3_xyze(2, 1) = 0.84755282581460467117;
    nids.push_back(7611);
    tri3_xyze(0, 2) = 0.067845411120786860981;
    tri3_xyze(1, 2) = 0.088417861203440931139;
    tri3_xyze(2, 2) = 0.84833359417544329517;
    nids.push_back(-1087);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.054684532864653072337;
    tri3_xyze(1, 0) = 0.094716389309733514668;
    tri3_xyze(2, 0) = 0.84911436253628147508;
    nids.push_back(7613);
    tri3_xyze(0, 1) = 0.051569762988242770496;
    tri3_xyze(1, 1) = 0.089321449629905846912;
    tri3_xyze(2, 1) = 0.84990133642126086677;
    nids.push_back(7615);
    tri3_xyze(0, 2) = 0.064130140527114048998;
    tri3_xyze(1, 2) = 0.083576026298795749736;
    tri3_xyze(2, 2) = 0.84950784947877133746;
    nids.push_back(-1088);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051569762988242770496;
    tri3_xyze(1, 0) = 0.089321449629905846912;
    tri3_xyze(2, 0) = 0.84990133642126086677;
    nids.push_back(7615);
    tri3_xyze(0, 1) = 0.072930658226334904581;
    tri3_xyze(1, 1) = 0.072930658226326466886;
    tri3_xyze(2, 1) = 0.84990133642126131086;
    nids.push_back(7865);
    tri3_xyze(0, 2) = 0.064130140527114048998;
    tri3_xyze(1, 2) = 0.083576026298795749736;
    tri3_xyze(2, 2) = 0.84950784947877133746;
    nids.push_back(-1088);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.072930658226334904581;
    tri3_xyze(1, 0) = 0.072930658226326466886;
    tri3_xyze(2, 0) = 0.84990133642126131086;
    nids.push_back(7865);
    tri3_xyze(0, 1) = 0.077335608029225483273;
    tri3_xyze(1, 1) = 0.077335608029217170478;
    tri3_xyze(2, 1) = 0.84911436253628203019;
    nids.push_back(7863);
    tri3_xyze(0, 2) = 0.064130140527114048998;
    tri3_xyze(1, 2) = 0.083576026298795749736;
    tri3_xyze(2, 2) = 0.84950784947877133746;
    nids.push_back(-1088);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.077335608029225483273;
    tri3_xyze(1, 0) = 0.077335608029217170478;
    tri3_xyze(2, 0) = 0.84911436253628203019;
    nids.push_back(7863);
    tri3_xyze(0, 1) = 0.054684532864653072337;
    tri3_xyze(1, 1) = 0.094716389309733514668;
    tri3_xyze(2, 1) = 0.84911436253628147508;
    nids.push_back(7613);
    tri3_xyze(0, 2) = 0.064130140527114048998;
    tri3_xyze(1, 2) = 0.083576026298795749736;
    tri3_xyze(2, 2) = 0.84950784947877133746;
    nids.push_back(-1088);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051569762988242770496;
    tri3_xyze(1, 0) = 0.089321449629905846912;
    tri3_xyze(2, 0) = 0.84990133642126086677;
    nids.push_back(7615);
    tri3_xyze(0, 1) = 0.048430237011777101108;
    tri3_xyze(1, 1) = 0.083883631126985019888;
    tri3_xyze(2, 1) = 0.84990133642126097779;
    nids.push_back(7617);
    tri3_xyze(0, 2) = 0.060355339059337327623;
    tri3_xyze(1, 2) = 0.078656609248550862779;
    tri3_xyze(2, 2) = 0.84990133642126108882;
    nids.push_back(-1089);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048430237011777101108;
    tri3_xyze(1, 0) = 0.083883631126985019888;
    tri3_xyze(2, 0) = 0.84990133642126097779;
    nids.push_back(7617);
    tri3_xyze(0, 1) = 0.068490698010994513489;
    tri3_xyze(1, 1) = 0.068490698010986075794;
    tri3_xyze(2, 1) = 0.84990133642126153291;
    nids.push_back(7867);
    tri3_xyze(0, 2) = 0.060355339059337327623;
    tri3_xyze(1, 2) = 0.078656609248550862779;
    tri3_xyze(2, 2) = 0.84990133642126108882;
    nids.push_back(-1089);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.068490698010994513489;
    tri3_xyze(1, 0) = 0.068490698010986075794;
    tri3_xyze(2, 0) = 0.84990133642126153291;
    nids.push_back(7867);
    tri3_xyze(0, 1) = 0.072930658226334904581;
    tri3_xyze(1, 1) = 0.072930658226326466886;
    tri3_xyze(2, 1) = 0.84990133642126131086;
    nids.push_back(7865);
    tri3_xyze(0, 2) = 0.060355339059337327623;
    tri3_xyze(1, 2) = 0.078656609248550862779;
    tri3_xyze(2, 2) = 0.84990133642126108882;
    nids.push_back(-1089);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.072930658226334904581;
    tri3_xyze(1, 0) = 0.072930658226326466886;
    tri3_xyze(2, 0) = 0.84990133642126131086;
    nids.push_back(7865);
    tri3_xyze(0, 1) = 0.051569762988242770496;
    tri3_xyze(1, 1) = 0.089321449629905846912;
    tri3_xyze(2, 1) = 0.84990133642126086677;
    nids.push_back(7615);
    tri3_xyze(0, 2) = 0.060355339059337327623;
    tri3_xyze(1, 2) = 0.078656609248550862779;
    tri3_xyze(2, 2) = 0.84990133642126108882;
    nids.push_back(-1089);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048430237011777101108;
    tri3_xyze(1, 0) = 0.083883631126985019888;
    tri3_xyze(2, 0) = 0.84990133642126097779;
    nids.push_back(7617);
    tri3_xyze(0, 1) = 0.045315467135366820084;
    tri3_xyze(1, 1) = 0.078488691447157338255;
    tri3_xyze(2, 1) = 0.84911436253628203019;
    nids.push_back(7619);
    tri3_xyze(0, 2) = 0.056580537591560578492;
    tri3_xyze(1, 2) = 0.073737192198305989699;
    tri3_xyze(2, 2) = 0.84950784947877178155;
    nids.push_back(-1090);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.045315467135366820084;
    tri3_xyze(1, 0) = 0.078488691447157338255;
    tri3_xyze(2, 0) = 0.84911436253628203019;
    nids.push_back(7619);
    tri3_xyze(0, 1) = 0.064085748208103893164;
    tri3_xyze(1, 1) = 0.064085748208095552614;
    tri3_xyze(2, 1) = 0.84911436253628247428;
    nids.push_back(7869);
    tri3_xyze(0, 2) = 0.056580537591560578492;
    tri3_xyze(1, 2) = 0.073737192198305989699;
    tri3_xyze(2, 2) = 0.84950784947877178155;
    nids.push_back(-1090);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.064085748208103893164;
    tri3_xyze(1, 0) = 0.064085748208095552614;
    tri3_xyze(2, 0) = 0.84911436253628247428;
    nids.push_back(7869);
    tri3_xyze(0, 1) = 0.068490698010994513489;
    tri3_xyze(1, 1) = 0.068490698010986075794;
    tri3_xyze(2, 1) = 0.84990133642126153291;
    nids.push_back(7867);
    tri3_xyze(0, 2) = 0.056580537591560578492;
    tri3_xyze(1, 2) = 0.073737192198305989699;
    tri3_xyze(2, 2) = 0.84950784947877178155;
    nids.push_back(-1090);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.068490698010994513489;
    tri3_xyze(1, 0) = 0.068490698010986075794;
    tri3_xyze(2, 0) = 0.84990133642126153291;
    nids.push_back(7867);
    tri3_xyze(0, 1) = 0.048430237011777101108;
    tri3_xyze(1, 1) = 0.083883631126985019888;
    tri3_xyze(2, 1) = 0.84990133642126097779;
    nids.push_back(7617);
    tri3_xyze(0, 2) = 0.056580537591560578492;
    tri3_xyze(1, 2) = 0.073737192198305989699;
    tri3_xyze(2, 2) = 0.84950784947877178155;
    nids.push_back(-1090);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.045315467135366820084;
    tri3_xyze(1, 0) = 0.078488691447157338255;
    tri3_xyze(2, 0) = 0.84911436253628203019;
    nids.push_back(7619);
    tri3_xyze(0, 1) = 0.042274575140636268256;
    tri3_xyze(1, 1) = 0.073221712011954515065;
    tri3_xyze(2, 1) = 0.84755282581460544833;
    nids.push_back(7621);
    tri3_xyze(0, 2) = 0.052865266997887794265;
    tri3_xyze(1, 2) = 0.068895357293660794418;
    tri3_xyze(2, 2) = 0.8483335941754439613;
    nids.push_back(-1091);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042274575140636268256;
    tri3_xyze(1, 0) = 0.073221712011954515065;
    tri3_xyze(2, 0) = 0.84755282581460544833;
    nids.push_back(7621);
    tri3_xyze(0, 1) = 0.059785277507444195555;
    tri3_xyze(1, 1) = 0.059785277507435799493;
    tri3_xyze(2, 1) = 0.84755282581460589242;
    nids.push_back(7871);
    tri3_xyze(0, 2) = 0.052865266997887794265;
    tri3_xyze(1, 2) = 0.068895357293660794418;
    tri3_xyze(2, 2) = 0.8483335941754439613;
    nids.push_back(-1091);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.059785277507444195555;
    tri3_xyze(1, 0) = 0.059785277507435799493;
    tri3_xyze(2, 0) = 0.84755282581460589242;
    nids.push_back(7871);
    tri3_xyze(0, 1) = 0.064085748208103893164;
    tri3_xyze(1, 1) = 0.064085748208095552614;
    tri3_xyze(2, 1) = 0.84911436253628247428;
    nids.push_back(7869);
    tri3_xyze(0, 2) = 0.052865266997887794265;
    tri3_xyze(1, 2) = 0.068895357293660794418;
    tri3_xyze(2, 2) = 0.8483335941754439613;
    nids.push_back(-1091);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.064085748208103893164;
    tri3_xyze(1, 0) = 0.064085748208095552614;
    tri3_xyze(2, 0) = 0.84911436253628247428;
    nids.push_back(7869);
    tri3_xyze(0, 1) = 0.045315467135366820084;
    tri3_xyze(1, 1) = 0.078488691447157338255;
    tri3_xyze(2, 1) = 0.84911436253628203019;
    nids.push_back(7619);
    tri3_xyze(0, 2) = 0.052865266997887794265;
    tri3_xyze(1, 2) = 0.068895357293660794418;
    tri3_xyze(2, 2) = 0.8483335941754439613;
    nids.push_back(-1091);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042274575140636268256;
    tri3_xyze(1, 0) = 0.073221712011954515065;
    tri3_xyze(2, 0) = 0.84755282581460544833;
    nids.push_back(7621);
    tri3_xyze(0, 1) = 0.039355517710883122784;
    tri3_xyze(1, 1) = 0.068165756233410548792;
    tri3_xyze(2, 1) = 0.8452413526231489227;
    nids.push_back(7623);
    tri3_xyze(0, 2) = 0.049268119314976245582;
    tri3_xyze(1, 2) = 0.064207463163433439735;
    tri3_xyze(2, 2) = 0.84639708921887735205;
    nids.push_back(-1092);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.039355517710883122784;
    tri3_xyze(1, 0) = 0.068165756233410548792;
    tri3_xyze(2, 0) = 0.8452413526231489227;
    nids.push_back(7623);
    tri3_xyze(0, 1) = 0.055657106900941374916;
    tri3_xyze(1, 1) = 0.055657106900932902527;
    tri3_xyze(2, 1) = 0.84524135262314925576;
    nids.push_back(7873);
    tri3_xyze(0, 2) = 0.049268119314976245582;
    tri3_xyze(1, 2) = 0.064207463163433439735;
    tri3_xyze(2, 2) = 0.84639708921887735205;
    nids.push_back(-1092);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.055657106900941374916;
    tri3_xyze(1, 0) = 0.055657106900932902527;
    tri3_xyze(2, 0) = 0.84524135262314925576;
    nids.push_back(7873);
    tri3_xyze(0, 1) = 0.059785277507444195555;
    tri3_xyze(1, 1) = 0.059785277507435799493;
    tri3_xyze(2, 1) = 0.84755282581460589242;
    nids.push_back(7871);
    tri3_xyze(0, 2) = 0.049268119314976245582;
    tri3_xyze(1, 2) = 0.064207463163433439735;
    tri3_xyze(2, 2) = 0.84639708921887735205;
    nids.push_back(-1092);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.059785277507444195555;
    tri3_xyze(1, 0) = 0.059785277507435799493;
    tri3_xyze(2, 0) = 0.84755282581460589242;
    nids.push_back(7871);
    tri3_xyze(0, 1) = 0.042274575140636268256;
    tri3_xyze(1, 1) = 0.073221712011954515065;
    tri3_xyze(2, 1) = 0.84755282581460544833;
    nids.push_back(7621);
    tri3_xyze(0, 2) = 0.049268119314976245582;
    tri3_xyze(1, 2) = 0.064207463163433439735;
    tri3_xyze(2, 2) = 0.84639708921887735205;
    nids.push_back(-1092);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.039355517710883122784;
    tri3_xyze(1, 0) = 0.068165756233410548792;
    tri3_xyze(2, 0) = 0.8452413526231489227;
    nids.push_back(7623);
    tri3_xyze(0, 1) = 0.036604330125535038121;
    tri3_xyze(1, 1) = 0.063400559554434840481;
    tri3_xyze(2, 1) = 0.84221639627494881974;
    nids.push_back(7625);
    tri3_xyze(0, 2) = 0.045845823710617281121;
    tri3_xyze(1, 2) = 0.059747440698469822529;
    tri3_xyze(2, 2) = 0.84372887444904898224;
    nids.push_back(-1093);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036604330125535038121;
    tri3_xyze(1, 0) = 0.063400559554434840481;
    tri3_xyze(2, 0) = 0.84221639627494881974;
    nids.push_back(7625);
    tri3_xyze(0, 1) = 0.051766340105109574787;
    tri3_xyze(1, 1) = 0.051766340105100998314;
    tri3_xyze(2, 1) = 0.84221639627494915281;
    nids.push_back(7875);
    tri3_xyze(0, 2) = 0.045845823710617281121;
    tri3_xyze(1, 2) = 0.059747440698469822529;
    tri3_xyze(2, 2) = 0.84372887444904898224;
    nids.push_back(-1093);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051766340105109574787;
    tri3_xyze(1, 0) = 0.051766340105100998314;
    tri3_xyze(2, 0) = 0.84221639627494915281;
    nids.push_back(7875);
    tri3_xyze(0, 1) = 0.055657106900941374916;
    tri3_xyze(1, 1) = 0.055657106900932902527;
    tri3_xyze(2, 1) = 0.84524135262314925576;
    nids.push_back(7873);
    tri3_xyze(0, 2) = 0.045845823710617281121;
    tri3_xyze(1, 2) = 0.059747440698469822529;
    tri3_xyze(2, 2) = 0.84372887444904898224;
    nids.push_back(-1093);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.055657106900941374916;
    tri3_xyze(1, 0) = 0.055657106900932902527;
    tri3_xyze(2, 0) = 0.84524135262314925576;
    nids.push_back(7873);
    tri3_xyze(0, 1) = 0.039355517710883122784;
    tri3_xyze(1, 1) = 0.068165756233410548792;
    tri3_xyze(2, 1) = 0.8452413526231489227;
    nids.push_back(7623);
    tri3_xyze(0, 2) = 0.045845823710617281121;
    tri3_xyze(1, 2) = 0.059747440698469822529;
    tri3_xyze(2, 2) = 0.84372887444904898224;
    nids.push_back(-1093);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036604330125535038121;
    tri3_xyze(1, 0) = 0.063400559554434840481;
    tri3_xyze(2, 0) = 0.84221639627494881974;
    nids.push_back(7625);
    tri3_xyze(0, 1) = 0.034064400256292717839;
    tri3_xyze(1, 1) = 0.059001271973245346392;
    tri3_xyze(2, 1) = 0.83852566213863766187;
    nids.push_back(7627);
    tri3_xyze(0, 2) = 0.042652351830871987337;
    tri3_xyze(1, 2) = 0.055585627117330782543;
    tri3_xyze(2, 2) = 0.84037102920679340734;
    nids.push_back(-1094);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.034064400256292717839;
    tri3_xyze(1, 0) = 0.059001271973245346392;
    tri3_xyze(2, 0) = 0.83852566213863766187;
    nids.push_back(7627);
    tri3_xyze(0, 1) = 0.04817433683655062554;
    tri3_xyze(1, 1) = 0.048174336836541924167;
    tri3_xyze(2, 1) = 0.83852566213863799494;
    nids.push_back(7877);
    tri3_xyze(0, 2) = 0.042652351830871987337;
    tri3_xyze(1, 2) = 0.055585627117330782543;
    tri3_xyze(2, 2) = 0.84037102920679340734;
    nids.push_back(-1094);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.04817433683655062554;
    tri3_xyze(1, 0) = 0.048174336836541924167;
    tri3_xyze(2, 0) = 0.83852566213863799494;
    nids.push_back(7877);
    tri3_xyze(0, 1) = 0.051766340105109574787;
    tri3_xyze(1, 1) = 0.051766340105100998314;
    tri3_xyze(2, 1) = 0.84221639627494915281;
    nids.push_back(7875);
    tri3_xyze(0, 2) = 0.042652351830871987337;
    tri3_xyze(1, 2) = 0.055585627117330782543;
    tri3_xyze(2, 2) = 0.84037102920679340734;
    nids.push_back(-1094);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.051766340105109574787;
    tri3_xyze(1, 0) = 0.051766340105100998314;
    tri3_xyze(2, 0) = 0.84221639627494915281;
    nids.push_back(7875);
    tri3_xyze(0, 1) = 0.036604330125535038121;
    tri3_xyze(1, 1) = 0.063400559554434840481;
    tri3_xyze(2, 1) = 0.84221639627494881974;
    nids.push_back(7625);
    tri3_xyze(0, 2) = 0.042652351830871987337;
    tri3_xyze(1, 2) = 0.055585627117330782543;
    tri3_xyze(2, 2) = 0.84037102920679340734;
    nids.push_back(-1094);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.034064400256292717839;
    tri3_xyze(1, 0) = 0.059001271973245346392;
    tri3_xyze(2, 0) = 0.83852566213863766187;
    nids.push_back(7627);
    tri3_xyze(0, 1) = 0.031775784314474669889;
    tri3_xyze(1, 1) = 0.055037272883004147161;
    tri3_xyze(2, 1) = 0.83422735529628266882;
    nids.push_back(7629);
    tri3_xyze(0, 2) = 0.039738066634971550006;
    tri3_xyze(1, 2) = 0.05178765670633769802;
    tri3_xyze(2, 2) = 0.83637650871746027637;
    nids.push_back(-1095);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.031775784314474669889;
    tri3_xyze(1, 0) = 0.055037272883004147161;
    tri3_xyze(2, 0) = 0.83422735529628266882;
    nids.push_back(7629);
    tri3_xyze(0, 1) = 0.044937745132568200634;
    tri3_xyze(1, 1) = 0.044937745132559360484;
    tri3_xyze(2, 1) = 0.83422735529628300188;
    nids.push_back(7879);
    tri3_xyze(0, 2) = 0.039738066634971550006;
    tri3_xyze(1, 2) = 0.05178765670633769802;
    tri3_xyze(2, 2) = 0.83637650871746027637;
    nids.push_back(-1095);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.044937745132568200634;
    tri3_xyze(1, 0) = 0.044937745132559360484;
    tri3_xyze(2, 0) = 0.83422735529628300188;
    nids.push_back(7879);
    tri3_xyze(0, 1) = 0.04817433683655062554;
    tri3_xyze(1, 1) = 0.048174336836541924167;
    tri3_xyze(2, 1) = 0.83852566213863799494;
    nids.push_back(7877);
    tri3_xyze(0, 2) = 0.039738066634971550006;
    tri3_xyze(1, 2) = 0.05178765670633769802;
    tri3_xyze(2, 2) = 0.83637650871746027637;
    nids.push_back(-1095);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.04817433683655062554;
    tri3_xyze(1, 0) = 0.048174336836541924167;
    tri3_xyze(2, 0) = 0.83852566213863799494;
    nids.push_back(7877);
    tri3_xyze(0, 1) = 0.034064400256292717839;
    tri3_xyze(1, 1) = 0.059001271973245346392;
    tri3_xyze(2, 1) = 0.83852566213863766187;
    nids.push_back(7627);
    tri3_xyze(0, 2) = 0.039738066634971550006;
    tri3_xyze(1, 2) = 0.05178765670633769802;
    tri3_xyze(2, 2) = 0.83637650871746027637;
    nids.push_back(-1095);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.031775784314474669889;
    tri3_xyze(1, 0) = 0.055037272883004147161;
    tri3_xyze(2, 0) = 0.83422735529628266882;
    nids.push_back(7629);
    tri3_xyze(0, 1) = 0.029774575140636295317;
    tri3_xyze(1, 1) = 0.051571076917342986712;
    tri3_xyze(2, 1) = 0.82938926261447210031;
    nids.push_back(7631);
    tri3_xyze(0, 2) = 0.037148928141364924438;
    tri3_xyze(1, 2) = 0.048413425727669504028;
    tri3_xyze(2, 2) = 0.83180830895537749559;
    nids.push_back(-1096);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.029774575140636295317;
    tri3_xyze(1, 0) = 0.051571076917342986712;
    tri3_xyze(2, 0) = 0.82938926261447210031;
    nids.push_back(7631);
    tri3_xyze(0, 1) = 0.042107607977780528441;
    tri3_xyze(1, 1) = 0.042107607977771528696;
    tri3_xyze(2, 1) = 0.82938926261447243338;
    nids.push_back(7881);
    tri3_xyze(0, 2) = 0.037148928141364924438;
    tri3_xyze(1, 2) = 0.048413425727669504028;
    tri3_xyze(2, 2) = 0.83180830895537749559;
    nids.push_back(-1096);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042107607977780528441;
    tri3_xyze(1, 0) = 0.042107607977771528696;
    tri3_xyze(2, 0) = 0.82938926261447243338;
    nids.push_back(7881);
    tri3_xyze(0, 1) = 0.044937745132568200634;
    tri3_xyze(1, 1) = 0.044937745132559360484;
    tri3_xyze(2, 1) = 0.83422735529628300188;
    nids.push_back(7879);
    tri3_xyze(0, 2) = 0.037148928141364924438;
    tri3_xyze(1, 2) = 0.048413425727669504028;
    tri3_xyze(2, 2) = 0.83180830895537749559;
    nids.push_back(-1096);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.044937745132568200634;
    tri3_xyze(1, 0) = 0.044937745132559360484;
    tri3_xyze(2, 0) = 0.83422735529628300188;
    nids.push_back(7879);
    tri3_xyze(0, 1) = 0.031775784314474669889;
    tri3_xyze(1, 1) = 0.055037272883004147161;
    tri3_xyze(2, 1) = 0.83422735529628266882;
    nids.push_back(7629);
    tri3_xyze(0, 2) = 0.037148928141364924438;
    tri3_xyze(1, 2) = 0.048413425727669504028;
    tri3_xyze(2, 2) = 0.83180830895537749559;
    nids.push_back(-1096);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.029774575140636295317;
    tri3_xyze(1, 0) = 0.051571076917342986712;
    tri3_xyze(2, 0) = 0.82938926261447210031;
    nids.push_back(7631);
    tri3_xyze(0, 1) = 0.028092332998913428777;
    tri3_xyze(1, 1) = 0.048657348057245225492;
    tri3_xyze(2, 1) = 0.82408768370493423738;
    nids.push_back(7633);
    tri3_xyze(0, 2) = 0.03492576861077266448;
    tri3_xyze(1, 2) = 0.045516147819527744134;
    tri3_xyze(2, 2) = 0.82673847315970327987;
    nids.push_back(-1097);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.028092332998913428777;
    tri3_xyze(1, 0) = 0.048657348057245225492;
    tri3_xyze(2, 0) = 0.82408768370493423738;
    nids.push_back(7633);
    tri3_xyze(0, 1) = 0.039728558325760408854;
    tri3_xyze(1, 1) = 0.039728558325751242575;
    tri3_xyze(2, 1) = 0.82408768370493457045;
    nids.push_back(7883);
    tri3_xyze(0, 2) = 0.03492576861077266448;
    tri3_xyze(1, 2) = 0.045516147819527744134;
    tri3_xyze(2, 2) = 0.82673847315970327987;
    nids.push_back(-1097);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.039728558325760408854;
    tri3_xyze(1, 0) = 0.039728558325751242575;
    tri3_xyze(2, 0) = 0.82408768370493457045;
    nids.push_back(7883);
    tri3_xyze(0, 1) = 0.042107607977780528441;
    tri3_xyze(1, 1) = 0.042107607977771528696;
    tri3_xyze(2, 1) = 0.82938926261447243338;
    nids.push_back(7881);
    tri3_xyze(0, 2) = 0.03492576861077266448;
    tri3_xyze(1, 2) = 0.045516147819527744134;
    tri3_xyze(2, 2) = 0.82673847315970327987;
    nids.push_back(-1097);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.042107607977780528441;
    tri3_xyze(1, 0) = 0.042107607977771528696;
    tri3_xyze(2, 0) = 0.82938926261447243338;
    nids.push_back(7881);
    tri3_xyze(0, 1) = 0.029774575140636295317;
    tri3_xyze(1, 1) = 0.051571076917342986712;
    tri3_xyze(2, 1) = 0.82938926261447210031;
    nids.push_back(7631);
    tri3_xyze(0, 2) = 0.03492576861077266448;
    tri3_xyze(1, 2) = 0.045516147819527744134;
    tri3_xyze(2, 2) = 0.82673847315970327987;
    nids.push_back(-1097);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.028092332998913428777;
    tri3_xyze(1, 0) = 0.048657348057245225492;
    tri3_xyze(2, 0) = 0.82408768370493423738;
    nids.push_back(7633);
    tri3_xyze(0, 1) = 0.026755587852803689591;
    tri3_xyze(1, 1) = 0.046342037547411993326;
    tri3_xyze(2, 1) = 0.81840622763408255036;
    nids.push_back(7635);
    tri3_xyze(0, 2) = 0.033103648597043303092;
    tri3_xyze(1, 2) = 0.043141514785273699095;
    tri3_xyze(2, 2) = 0.82124695566950844938;
    nids.push_back(-1098);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026755587852803689591;
    tri3_xyze(1, 0) = 0.046342037547411993326;
    tri3_xyze(2, 0) = 0.81840622763408255036;
    nids.push_back(7635);
    tri3_xyze(0, 1) = 0.037838115210695681678;
    tri3_xyze(1, 1) = 0.037838115210686328049;
    tri3_xyze(2, 1) = 0.8184062276340827724;
    nids.push_back(7885);
    tri3_xyze(0, 2) = 0.033103648597043303092;
    tri3_xyze(1, 2) = 0.043141514785273699095;
    tri3_xyze(2, 2) = 0.82124695566950844938;
    nids.push_back(-1098);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.037838115210695681678;
    tri3_xyze(1, 0) = 0.037838115210686328049;
    tri3_xyze(2, 0) = 0.8184062276340827724;
    nids.push_back(7885);
    tri3_xyze(0, 1) = 0.039728558325760408854;
    tri3_xyze(1, 1) = 0.039728558325751242575;
    tri3_xyze(2, 1) = 0.82408768370493457045;
    nids.push_back(7883);
    tri3_xyze(0, 2) = 0.033103648597043303092;
    tri3_xyze(1, 2) = 0.043141514785273699095;
    tri3_xyze(2, 2) = 0.82124695566950844938;
    nids.push_back(-1098);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.039728558325760408854;
    tri3_xyze(1, 0) = 0.039728558325751242575;
    tri3_xyze(2, 0) = 0.82408768370493457045;
    nids.push_back(7883);
    tri3_xyze(0, 1) = 0.028092332998913428777;
    tri3_xyze(1, 1) = 0.048657348057245225492;
    tri3_xyze(2, 1) = 0.82408768370493423738;
    nids.push_back(7633);
    tri3_xyze(0, 2) = 0.033103648597043303092;
    tri3_xyze(1, 2) = 0.043141514785273699095;
    tri3_xyze(2, 2) = 0.82124695566950844938;
    nids.push_back(-1098);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.026755587852803689591;
    tri3_xyze(1, 0) = 0.046342037547411993326;
    tri3_xyze(2, 0) = 0.81840622763408255036;
    nids.push_back(7635);
    tri3_xyze(0, 1) = 0.025785420971794200534;
    tri3_xyze(1, 1) = 0.044661659217682739276;
    tri3_xyze(2, 1) = 0.8124344943580913414;
    nids.push_back(7637);
    tri3_xyze(0, 2) = 0.031711304021275116272;
    tri3_xyze(1, 2) = 0.041326976006394587626;
    tri3_xyze(2, 2) = 0.8154203609960870569;
    nids.push_back(-1099);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025785420971794200534;
    tri3_xyze(1, 0) = 0.044661659217682739276;
    tri3_xyze(2, 0) = 0.8124344943580913414;
    nids.push_back(7637);
    tri3_xyze(0, 1) = 0.036466092049806886344;
    tri3_xyze(1, 1) = 0.036466092049797275976;
    tri3_xyze(2, 1) = 0.81243449435809167447;
    nids.push_back(7887);
    tri3_xyze(0, 2) = 0.031711304021275116272;
    tri3_xyze(1, 2) = 0.041326976006394587626;
    tri3_xyze(2, 2) = 0.8154203609960870569;
    nids.push_back(-1099);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036466092049806886344;
    tri3_xyze(1, 0) = 0.036466092049797275976;
    tri3_xyze(2, 0) = 0.81243449435809167447;
    nids.push_back(7887);
    tri3_xyze(0, 1) = 0.037838115210695681678;
    tri3_xyze(1, 1) = 0.037838115210686328049;
    tri3_xyze(2, 1) = 0.8184062276340827724;
    nids.push_back(7885);
    tri3_xyze(0, 2) = 0.031711304021275116272;
    tri3_xyze(1, 2) = 0.041326976006394587626;
    tri3_xyze(2, 2) = 0.8154203609960870569;
    nids.push_back(-1099);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.037838115210695681678;
    tri3_xyze(1, 0) = 0.037838115210686328049;
    tri3_xyze(2, 0) = 0.8184062276340827724;
    nids.push_back(7885);
    tri3_xyze(0, 1) = 0.026755587852803689591;
    tri3_xyze(1, 1) = 0.046342037547411993326;
    tri3_xyze(2, 1) = 0.81840622763408255036;
    nids.push_back(7635);
    tri3_xyze(0, 2) = 0.031711304021275116272;
    tri3_xyze(1, 2) = 0.041326976006394587626;
    tri3_xyze(2, 2) = 0.8154203609960870569;
    nids.push_back(-1099);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025785420971794200534;
    tri3_xyze(1, 0) = 0.044661659217682739276;
    tri3_xyze(2, 0) = 0.8124344943580913414;
    nids.push_back(7637);
    tri3_xyze(0, 1) = 0.025197132467147994456;
    tri3_xyze(1, 1) = 0.043642713638126676146;
    tri3_xyze(2, 1) = 0.80626666167806393393;
    nids.push_back(7639);
    tri3_xyze(0, 2) = 0.030770692989174282944;
    tri3_xyze(1, 2) = 0.040101147843386253378;
    tri3_xyze(2, 2) = 0.80935057801807785971;
    nids.push_back(-1100);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.025197132467147994456;
    tri3_xyze(1, 0) = 0.043642713638126676146;
    tri3_xyze(2, 0) = 0.80626666167806393393;
    nids.push_back(7639);
    tri3_xyze(0, 1) = 0.035634126467948050443;
    tri3_xyze(1, 1) = 0.035634126467938308236;
    tri3_xyze(2, 1) = 0.80626666167806415597;
    nids.push_back(7889);
    tri3_xyze(0, 2) = 0.030770692989174282944;
    tri3_xyze(1, 2) = 0.040101147843386253378;
    tri3_xyze(2, 2) = 0.80935057801807785971;
    nids.push_back(-1100);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.035634126467948050443;
    tri3_xyze(1, 0) = 0.035634126467938308236;
    tri3_xyze(2, 0) = 0.80626666167806415597;
    nids.push_back(7889);
    tri3_xyze(0, 1) = 0.036466092049806886344;
    tri3_xyze(1, 1) = 0.036466092049797275976;
    tri3_xyze(2, 1) = 0.81243449435809167447;
    nids.push_back(7887);
    tri3_xyze(0, 2) = 0.030770692989174282944;
    tri3_xyze(1, 2) = 0.040101147843386253378;
    tri3_xyze(2, 2) = 0.80935057801807785971;
    nids.push_back(-1100);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.036466092049806886344;
    tri3_xyze(1, 0) = 0.036466092049797275976;
    tri3_xyze(2, 0) = 0.81243449435809167447;
    nids.push_back(7887);
    tri3_xyze(0, 1) = 0.025785420971794200534;
    tri3_xyze(1, 1) = 0.044661659217682739276;
    tri3_xyze(2, 1) = 0.8124344943580913414;
    nids.push_back(7637);
    tri3_xyze(0, 2) = 0.030770692989174282944;
    tri3_xyze(1, 2) = 0.040101147843386253378;
    tri3_xyze(2, 2) = 0.80935057801807785971;
    nids.push_back(-1100);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.035634126467948050443;
    tri3_xyze(1, 0) = 0.035634126467938308236;
    tri3_xyze(2, 0) = 0.80626666167806415597;
    nids.push_back(7889);
    tri3_xyze(0, 1) = 0.035355339059337340113;
    tri3_xyze(1, 1) = 0.035355339059327417495;
    tri3_xyze(2, 1) = 0.79999999999984894306;
    nids.push_back(7791);
    tri3_xyze(0, 2) = 0.039483362338663431645;
    tri3_xyze(1, 2) = 0.030296649498601013056;
    tri3_xyze(2, 2) = 0.80313333083895677156;
    nids.push_back(-1102);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.035355339059337340113;
    tri3_xyze(1, 0) = 0.035355339059327417495;
    tri3_xyze(2, 0) = 0.79999999999984894306;
    nids.push_back(7791);
    tri3_xyze(0, 1) = 0.043301270189231896857;
    tri3_xyze(1, 1) = 0.025000000000000056899;
    tri3_xyze(2, 1) = 0.79999999999984927612;
    nids.push_back(8041);
    tri3_xyze(0, 2) = 0.039483362338663431645;
    tri3_xyze(1, 2) = 0.030296649498601013056;
    tri3_xyze(2, 2) = 0.80313333083895677156;
    nids.push_back(-1102);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.043301270189231896857;
    tri3_xyze(1, 0) = 0.025000000000000056899;
    tri3_xyze(2, 0) = 0.79999999999984927612;
    nids.push_back(8041);
    tri3_xyze(0, 1) = 0.043642713638136439169;
    tri3_xyze(1, 1) = 0.025197132467138269596;
    tri3_xyze(2, 1) = 0.80626666167806460006;
    nids.push_back(8139);
    tri3_xyze(0, 2) = 0.039483362338663431645;
    tri3_xyze(1, 2) = 0.030296649498601013056;
    tri3_xyze(2, 2) = 0.80313333083895677156;
    nids.push_back(-1102);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.043642713638136439169;
    tri3_xyze(1, 0) = 0.025197132467138269596;
    tri3_xyze(2, 0) = 0.80626666167806460006;
    nids.push_back(8139);
    tri3_xyze(0, 1) = 0.035634126467948050443;
    tri3_xyze(1, 1) = 0.035634126467938308236;
    tri3_xyze(2, 1) = 0.80626666167806415597;
    nids.push_back(7889);
    tri3_xyze(0, 2) = 0.039483362338663431645;
    tri3_xyze(1, 2) = 0.030296649498601013056;
    tri3_xyze(2, 2) = 0.80313333083895677156;
    nids.push_back(-1102);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.048174336836550576968;
    tri3_xyze(1, 0) = 0.048174336836539370654;
    tri3_xyze(2, 0) = 0.76147433786105911402;
    nids.push_back(7805);
    tri3_xyze(0, 1) = 0.051766340105109630298;
    tri3_xyze(1, 1) = 0.051766340105098403168;
    tri3_xyze(2, 1) = 0.75778360372474762308;
    nids.push_back(7807);
    tri3_xyze(0, 2) = 0.05558562711733945616;
    tri3_xyze(1, 2) = 0.042652351830860815718;
    tri3_xyze(2, 2) = 0.75962897079290359059;
    nids.push_back(-1109);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.050000000000000016653;
    hex8_xyze(1, 0) = 0.050000000000000016653;
    hex8_xyze(2, 0) = 0.75;
    nids.push_back(1865);
    hex8_xyze(0, 1) = 0.050000000000000016653;
    hex8_xyze(1, 1) = 0.10000000000000000555;
    hex8_xyze(2, 1) = 0.75;
    nids.push_back(1866);
    hex8_xyze(0, 2) = 0;
    hex8_xyze(1, 2) = 0.10000000000000000555;
    hex8_xyze(2, 2) = 0.75;
    nids.push_back(1877);
    hex8_xyze(0, 3) = 0;
    hex8_xyze(1, 3) = 0.050000000000000016653;
    hex8_xyze(2, 3) = 0.75;
    nids.push_back(1876);
    hex8_xyze(0, 4) = 0.050000000000000016653;
    hex8_xyze(1, 4) = 0.050000000000000016653;
    hex8_xyze(2, 4) = 0.80000000000000015543;
    nids.push_back(1986);
    hex8_xyze(0, 5) = 0.050000000000000016653;
    hex8_xyze(1, 5) = 0.099999999999999977796;
    hex8_xyze(2, 5) = 0.80000000000000015543;
    nids.push_back(1987);
    hex8_xyze(0, 6) = 0;
    hex8_xyze(1, 6) = 0.10000000000000000555;
    hex8_xyze(2, 6) = 0.80000000000000015543;
    nids.push_back(1998);
    hex8_xyze(0, 7) = 0;
    hex8_xyze(1, 7) = 0.050000000000000016653;
    hex8_xyze(2, 7) = 0.80000000000000015543;
    nids.push_back(1997);

    intersection.add_element(6922, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.10000000000000000555;
    hex8_xyze(1, 0) = 0.050000000000000016653;
    hex8_xyze(2, 0) = 0.80000000000000015543;
    nids.push_back(1975);
    hex8_xyze(0, 1) = 0.10000000000000000555;
    hex8_xyze(1, 1) = 0.099999999999999977796;
    hex8_xyze(2, 1) = 0.80000000000000015543;
    nids.push_back(1976);
    hex8_xyze(0, 2) = 0.050000000000000016653;
    hex8_xyze(1, 2) = 0.099999999999999977796;
    hex8_xyze(2, 2) = 0.80000000000000015543;
    nids.push_back(1987);
    hex8_xyze(0, 3) = 0.050000000000000016653;
    hex8_xyze(1, 3) = 0.050000000000000016653;
    hex8_xyze(2, 3) = 0.80000000000000015543;
    nids.push_back(1986);
    hex8_xyze(0, 4) = 0.10000000000000003331;
    hex8_xyze(1, 4) = 0.050000000000000016653;
    hex8_xyze(2, 4) = 0.84999999999999986677;
    nids.push_back(2096);
    hex8_xyze(0, 5) = 0.10000000000000000555;
    hex8_xyze(1, 5) = 0.10000000000000000555;
    hex8_xyze(2, 5) = 0.85000000000000008882;
    nids.push_back(2097);
    hex8_xyze(0, 6) = 0.050000000000000002776;
    hex8_xyze(1, 6) = 0.10000000000000000555;
    hex8_xyze(2, 6) = 0.84999999999999986677;
    nids.push_back(2108);
    hex8_xyze(0, 7) = 0.050000000000000016653;
    hex8_xyze(1, 7) = 0.050000000000000016653;
    hex8_xyze(2, 7) = 0.84999999999999986677;
    nids.push_back(2107);

    intersection.add_element(7012, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.050000000000000016653;
    hex8_xyze(1, 0) = -2.0816681711721688214e-18;
    hex8_xyze(2, 0) = 0.80000000000000015543;
    nids.push_back(1985);
    hex8_xyze(0, 1) = 0.050000000000000016653;
    hex8_xyze(1, 1) = 0.050000000000000016653;
    hex8_xyze(2, 1) = 0.80000000000000015543;
    nids.push_back(1986);
    hex8_xyze(0, 2) = 0;
    hex8_xyze(1, 2) = 0.050000000000000016653;
    hex8_xyze(2, 2) = 0.80000000000000015543;
    nids.push_back(1997);
    hex8_xyze(0, 3) = 0;
    hex8_xyze(1, 3) = 0;
    hex8_xyze(2, 3) = 0.80000000000000015543;
    nids.push_back(1996);
    hex8_xyze(0, 4) = 0.050000000000000016653;
    hex8_xyze(1, 4) = -2.4286128663675290077e-18;
    hex8_xyze(2, 4) = 0.8499999999999999778;
    nids.push_back(2106);
    hex8_xyze(0, 5) = 0.050000000000000016653;
    hex8_xyze(1, 5) = 0.050000000000000016653;
    hex8_xyze(2, 5) = 0.84999999999999986677;
    nids.push_back(2107);
    hex8_xyze(0, 6) = 0;
    hex8_xyze(1, 6) = 0.050000000000000016653;
    hex8_xyze(2, 6) = 0.8499999999999999778;
    nids.push_back(2118);
    hex8_xyze(0, 7) = 0;
    hex8_xyze(1, 7) = 0;
    hex8_xyze(2, 7) = 0.8499999999999999778;
    nids.push_back(2117);

    intersection.add_element(7021, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.050000000000000016653;
    hex8_xyze(1, 0) = 0.050000000000000016653;
    hex8_xyze(2, 0) = 0.80000000000000015543;
    nids.push_back(1986);
    hex8_xyze(0, 1) = 0.050000000000000016653;
    hex8_xyze(1, 1) = 0.099999999999999977796;
    hex8_xyze(2, 1) = 0.80000000000000015543;
    nids.push_back(1987);
    hex8_xyze(0, 2) = 0;
    hex8_xyze(1, 2) = 0.10000000000000000555;
    hex8_xyze(2, 2) = 0.80000000000000015543;
    nids.push_back(1998);
    hex8_xyze(0, 3) = 0;
    hex8_xyze(1, 3) = 0.050000000000000016653;
    hex8_xyze(2, 3) = 0.80000000000000015543;
    nids.push_back(1997);
    hex8_xyze(0, 4) = 0.050000000000000016653;
    hex8_xyze(1, 4) = 0.050000000000000016653;
    hex8_xyze(2, 4) = 0.84999999999999986677;
    nids.push_back(2107);
    hex8_xyze(0, 5) = 0.050000000000000002776;
    hex8_xyze(1, 5) = 0.10000000000000000555;
    hex8_xyze(2, 5) = 0.84999999999999986677;
    nids.push_back(2108);
    hex8_xyze(0, 6) = 0;
    hex8_xyze(1, 6) = 0.10000000000000000555;
    hex8_xyze(2, 6) = 0.8499999999999999778;
    nids.push_back(2119);
    hex8_xyze(0, 7) = 0;
    hex8_xyze(1, 7) = 0.050000000000000016653;
    hex8_xyze(2, 7) = 0.8499999999999999778;
    nids.push_back(2118);

    intersection.add_element(7022, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.050000000000000016653;
    hex8_xyze(1, 0) = 0.099999999999999977796;
    hex8_xyze(2, 0) = 0.80000000000000015543;
    nids.push_back(1987);
    hex8_xyze(0, 1) = 0.050000000000000016653;
    hex8_xyze(1, 1) = 0.1500000000000000222;
    hex8_xyze(2, 1) = 0.80000000000000015543;
    nids.push_back(1988);
    hex8_xyze(0, 2) = 0;
    hex8_xyze(1, 2) = 0.15000000000000004996;
    hex8_xyze(2, 2) = 0.80000000000000015543;
    nids.push_back(1999);
    hex8_xyze(0, 3) = 0;
    hex8_xyze(1, 3) = 0.10000000000000000555;
    hex8_xyze(2, 3) = 0.80000000000000015543;
    nids.push_back(1998);
    hex8_xyze(0, 4) = 0.050000000000000002776;
    hex8_xyze(1, 4) = 0.10000000000000000555;
    hex8_xyze(2, 4) = 0.84999999999999986677;
    nids.push_back(2108);
    hex8_xyze(0, 5) = 0.050000000000000002776;
    hex8_xyze(1, 5) = 0.1500000000000000222;
    hex8_xyze(2, 5) = 0.8499999999999999778;
    nids.push_back(2109);
    hex8_xyze(0, 6) = 0;
    hex8_xyze(1, 6) = 0.15000000000000004996;
    hex8_xyze(2, 6) = 0.8499999999999999778;
    nids.push_back(2120);
    hex8_xyze(0, 7) = 0;
    hex8_xyze(1, 7) = 0.10000000000000000555;
    hex8_xyze(2, 7) = 0.8499999999999999778;
    nids.push_back(2119);

    intersection.add_element(7023, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0;
    hex8_xyze(1, 0) = 0.050000000000000016653;
    hex8_xyze(2, 0) = 0.80000000000000015543;
    nids.push_back(1997);
    hex8_xyze(0, 1) = 0;
    hex8_xyze(1, 1) = 0.10000000000000000555;
    hex8_xyze(2, 1) = 0.80000000000000015543;
    nids.push_back(1998);
    hex8_xyze(0, 2) = -0.050000000000000016653;
    hex8_xyze(1, 2) = 0.10000000000000000555;
    hex8_xyze(2, 2) = 0.80000000000000015543;
    nids.push_back(2009);
    hex8_xyze(0, 3) = -0.050000000000000030531;
    hex8_xyze(1, 3) = 0.050000000000000016653;
    hex8_xyze(2, 3) = 0.80000000000000015543;
    nids.push_back(2008);
    hex8_xyze(0, 4) = 0;
    hex8_xyze(1, 4) = 0.050000000000000016653;
    hex8_xyze(2, 4) = 0.8499999999999999778;
    nids.push_back(2118);
    hex8_xyze(0, 5) = 0;
    hex8_xyze(1, 5) = 0.10000000000000000555;
    hex8_xyze(2, 5) = 0.8499999999999999778;
    nids.push_back(2119);
    hex8_xyze(0, 6) = -0.050000000000000016653;
    hex8_xyze(1, 6) = 0.10000000000000000555;
    hex8_xyze(2, 6) = 0.84999999999999986677;
    nids.push_back(2130);
    hex8_xyze(0, 7) = -0.050000000000000016653;
    hex8_xyze(1, 7) = 0.050000000000000002776;
    hex8_xyze(2, 7) = 0.85000000000000008882;
    nids.push_back(2129);

    intersection.add_element(7032, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  intersection.cut_test_cut(
      true, Cut::VCellGaussPts_DirectDivergence, Cut::BCellGaussPts_Tessellation);
  intersection.cut_finalize(
      true, Cut::VCellGaussPts_DirectDivergence, Cut::BCellGaussPts_Tessellation, false, true);
}
