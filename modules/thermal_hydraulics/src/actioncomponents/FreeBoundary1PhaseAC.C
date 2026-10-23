//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FreeBoundary1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"
#include "MooseUtils.h"

// Registered under the classic Component's name ('FreeBoundary1Phase'); see FlowChannel1PhaseAC.C
// for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           FreeBoundary1PhaseAC,
                           "FreeBoundary1Phase",
                           "add_user_object");
registerMooseActionAliased("ThermalHydraulicsApp",
                           FreeBoundary1PhaseAC,
                           "FreeBoundary1Phase",
                           "add_bc");
registerActionComponentAliased("ThermalHydraulicsApp", FreeBoundary1PhaseAC, "FreeBoundary1Phase");

InputParameters
FreeBoundary1PhaseAC::validParams()
{
  InputParameters params = FlowBoundary1PhaseAC::validParams();
  params.addClassDescription("Component to create a free flow boundary for single-phase flow.");
  return params;
}

FreeBoundary1PhaseAC::FreeBoundary1PhaseAC(const InputParameters & params)
  : ActionComponent(params), FlowBoundary1PhaseAC(params)
{
}

void
FreeBoundary1PhaseAC::addUserObjects()
{
  ExecFlagEnum userobject_execute_on(MooseUtils::getDefaultExecFlagEnum());
  userobject_execute_on = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};

  const std::string class_name = "ADBoundaryFlux3EqnFreeOutflow";
  InputParameters params = _factory.getValidParams(class_name);
  params.set<UserObjectName>("fluid_properties") = connectedPhysics().fluidPropertiesName();
  params.set<ExecFlagEnum>("execute_on") = userobject_execute_on;
  getProblem().addUserObject(class_name, _boundary_uo_name, params);
}
