//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HSBoundarySpecifiedTemperatureAC.h"

// Registered under the classic Component's name; see FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundarySpecifiedTemperatureAC,
                           "HSBoundarySpecifiedTemperature",
                           "add_bc");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundarySpecifiedTemperatureAC,
                           "HSBoundarySpecifiedTemperature",
                           "check_integrity");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HSBoundarySpecifiedTemperatureAC,
                               "HSBoundarySpecifiedTemperature");

InputParameters
HSBoundarySpecifiedTemperatureAC::validParams()
{
  InputParameters params = HSBoundaryAC::validParams();
  params.addRequiredParam<FunctionName>("T", "Prescribed temperature [K]");
  params.addClassDescription("Applies Dirichlet boundary conditions on a heat structure");
  return params;
}

HSBoundarySpecifiedTemperatureAC::HSBoundarySpecifiedTemperatureAC(const InputParameters & params)
  : ActionComponent(params), HSBoundaryAC(params)
{
}

void
HSBoundarySpecifiedTemperatureAC::addBCs()
{
  const std::string class_name = "ADFunctionDirichletBC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") = temperatureName();
  params.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params.set<FunctionName>("function") = getParam<FunctionName>("T");
  getProblem().addBoundaryCondition(class_name, genName(name(), "bc"), params);
}
