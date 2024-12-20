// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fsi_nox_fixpoint.hpp"

#include <NOX_Abstract_Group.H>
#include <NOX_GlobalData.H>

FOUR_C_NAMESPACE_OPEN

NOX::FSI::FixPoint::FixPoint(
    const Teuchos::RCP<::NOX::Utils>& utils, Teuchos::ParameterList& params)
    : utils_(utils)
{
}



bool NOX::FSI::FixPoint::reset(
    const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params)
{
  utils_ = gd->getUtils();
  return true;
}


bool NOX::FSI::FixPoint::compute(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& group,
    const ::NOX::Solver::Generic& solver)
{
  ::NOX::Abstract::Group::ReturnType status;

  // Compute F at current solution
  status = group.computeF();
  if (status != ::NOX::Abstract::Group::Ok) throw_error("compute", "Unable to compute F");

  // The residual is the direction.
  dir.update(1.0, group.getF());

  return true;
}


bool NOX::FSI::FixPoint::compute(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& group,
    const ::NOX::Solver::LineSearchBased& solver)
{
  return ::NOX::Direction::Generic::compute(dir, group, solver);
}


void NOX::FSI::FixPoint::throw_error(const std::string& functionName, const std::string& errorMsg)
{
  if (utils_->isPrintType(::NOX::Utils::Error))
    utils_->err() << "FixPoint::" << functionName << " - " << errorMsg << std::endl;
  throw "NOX Error";
}

FOUR_C_NAMESPACE_CLOSE
