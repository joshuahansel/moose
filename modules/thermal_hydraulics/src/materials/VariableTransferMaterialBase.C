//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableTransferMaterialBase.h"
#include "SubProblem.h"
#include "SystemBase.h"
#include "NearestNodeLocator.h"
#include "PenetrationLocator.h"
#include "Assembly.h"
#include "ADUtils.h"
#include "MooseVariableBase.h"

#include "libmesh/node.h"
#include "libmesh/elem.h"

registerMooseObject("ThermalHydraulicsApp", VariableTransferMaterialBase);

InputParameters
VariableTransferMaterialBase::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<BoundaryName>("primary_boundary",
                                        "The boundary of the 2D structure  to get the value from.");
  params.addRequiredParam<BoundaryName>(
      "secondary_boundary",
      "The boundary coincident with the 1D domain block that we are transferring data to (e.g. the "
      "block this material is executing on).");

  // Have to use std::string to circumvent block restrictable testing
  params.addRequiredParam<std::string>("paired_variable", "The variable to get the value of.");

  return params;
}

VariableTransferMaterialBase::VariableTransferMaterialBase(const InputParameters & parameters)
  : Material(parameters),
    _penetration_locator(getPenetrationLocator(getParam<BoundaryName>("primary_boundary"),
                                               getParam<BoundaryName>("secondary_boundary"),
                                               Order(FIRST))),
    _nearest_node(_penetration_locator._nearest_node),
    _nl_sys(_subproblem.systemBaseNonlinear()),
    _serialized_solution(_nl_sys.currentSolution()),
    _paired_variable(
        _subproblem
            .getVariable(_tid, getParam<std::string>("paired_variable"), Moose::VAR_NONLINEAR)
            .number())
{
  _penetration_locator.setCheckWhetherReasonable(false);
}

std::vector<ADReal>
VariableTransferMaterialBase::getNodalVariableValues() const
{
  std::vector<ADReal> nodal_values;
  for (const auto i : _current_elem->node_index_range())
  {
    const Node & nd = _current_elem->node_ref(i);

    // Assumes the variable you are coupling to is from the nonlinear system for now.
    const Node * const nearest = _nearest_node.nearestNode(nd.id());
    mooseAssert(nearest, "I do not have the nearest node for you");
    const auto dof_number = nearest->dof_number(_nl_sys.number(), _paired_variable, 0);
    nodal_values.push_back((*_serialized_solution)(dof_number));
    Moose::derivInsert(nodal_values.back().derivatives(), dof_number, 1.0);
  }

  return nodal_values;
}
