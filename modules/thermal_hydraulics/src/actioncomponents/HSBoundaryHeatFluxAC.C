//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HSBoundaryHeatFluxAC.h"

// Registered under the classic Component's name; see FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryHeatFluxAC,
                           "HSBoundaryHeatFlux",
                           "add_bc");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryHeatFluxAC,
                           "HSBoundaryHeatFlux",
                           "add_postprocessor");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryHeatFluxAC,
                           "HSBoundaryHeatFlux",
                           "check_integrity");
registerActionComponentAliased("ThermalHydraulicsApp", HSBoundaryHeatFluxAC, "HSBoundaryHeatFlux");

InputParameters
HSBoundaryHeatFluxAC::validParams()
{
  InputParameters params = HSBoundaryAC::validParams();
  params.addRequiredParam<FunctionName>("q", "Heat flux [W/m^2]");
  params.addParam<FunctionName>("scale", 1.0, "Function by which to scale the boundary condition");
  params.addParam<bool>(
      "scale_heat_rate_pp",
      true,
      "If true, the scaling function is applied to the heat rate post-processor.");
  params.addClassDescription("Applies a specified heat flux to a heat structure boundary");
  return params;
}

HSBoundaryHeatFluxAC::HSBoundaryHeatFluxAC(const InputParameters & params)
  : ActionComponent(params), HSBoundaryAC(params)
{
  addRequiredTask("add_postprocessor");
}

void
HSBoundaryHeatFluxAC::addBCs()
{
  const std::string class_name = "ADHSHeatFluxBC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") = temperatureName();
  params.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params.set<FunctionName>("function") = getParam<FunctionName>("q");
  params.set<FunctionName>("scale") = getParam<FunctionName>("scale");
  getProblem().addBoundaryCondition(class_name, genName(name(), "bc"), params);
}

void
HSBoundaryHeatFluxAC::addPostprocessors()
{
  // Only cylindrical heat structures get a heat rate post-processor, matching classic
  // HSBoundaryHeatFlux (a plate's integral would be per unit depth)
  if (!cylindricalHS())
    return;

  const std::string class_name = "HeatRateHeatFlux";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params.set<FunctionName>("function") = getParam<FunctionName>("q");
  if (getParam<bool>("scale_heat_rate_pp"))
    params.set<FunctionName>("scale") = getParam<FunctionName>("scale");
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
  getProblem().addPostprocessor(class_name, genSafeName(name(), "integral"), params);
}
