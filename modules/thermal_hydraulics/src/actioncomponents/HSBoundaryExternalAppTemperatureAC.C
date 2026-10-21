//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HSBoundaryExternalAppTemperatureAC.h"

// Registered under the classic Component's name; see FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryExternalAppTemperatureAC,
                           "HSBoundaryExternalAppTemperature",
                           "add_aux_variable");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryExternalAppTemperatureAC,
                           "HSBoundaryExternalAppTemperature",
                           "add_bc");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryExternalAppTemperatureAC,
                           "HSBoundaryExternalAppTemperature",
                           "check_integrity");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HSBoundaryExternalAppTemperatureAC,
                               "HSBoundaryExternalAppTemperature");

InputParameters
HSBoundaryExternalAppTemperatureAC::validParams()
{
  InputParameters params = HSBoundaryAC::validParams();
  params.addParam<VariableName>(
      "T_ext",
      "T_ext",
      "Name of the variable that will store the values computed by the external application");
  params.addClassDescription("Heat structure boundary condition to set temperature values computed "
                             "by an external application");
  return params;
}

HSBoundaryExternalAppTemperatureAC::HSBoundaryExternalAppTemperatureAC(
    const InputParameters & params)
  : ActionComponent(params), HSBoundaryAC(params)
{
  addRequiredTask("add_aux_variable");
}

void
HSBoundaryExternalAppTemperatureAC::addAuxVariables()
{
  requestAuxVariable(getParam<VariableName>("T_ext"), "LAGRANGE", "FIRST");
}

void
HSBoundaryExternalAppTemperatureAC::addBCs()
{
  const std::string class_name = "ADMatchedValueBC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") = temperatureName();
  params.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params.set<std::vector<VariableName>>("v") = {getParam<VariableName>("T_ext")};
  getProblem().addBoundaryCondition(class_name, genName(name(), "bc"), params);
}
