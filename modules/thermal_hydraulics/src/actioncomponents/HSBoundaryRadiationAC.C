//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HSBoundaryRadiationAC.h"

// Registered under the classic Component's name; see FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryRadiationAC,
                           "HSBoundaryRadiation",
                           "add_bc");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryRadiationAC,
                           "HSBoundaryRadiation",
                           "add_postprocessor");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryRadiationAC,
                           "HSBoundaryRadiation",
                           "check_integrity");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HSBoundaryRadiationAC,
                               "HSBoundaryRadiation");

InputParameters
HSBoundaryRadiationAC::validParams()
{
  InputParameters params = HSBoundaryAC::validParams();
  params.addRequiredParam<MooseFunctorName>("emissivity", "Emissivity functor [-]");
  params.addParam<MooseFunctorName>("view_factor", 1.0, "View factor functor [-]");
  params.addRequiredParam<MooseFunctorName>("T_ambient", "Ambient temperature functor [K]");
  params.addParam<MooseFunctorName>(
      "scale", 1.0, "Functor by which to scale the boundary condition");
  params.addParam<bool>(
      "scale_heat_rate_pp",
      true,
      "If true, the scaling function is applied to the heat rate post-processor.");
  params.addClassDescription("Radiative heat transfer boundary condition for heat structure");
  return params;
}

HSBoundaryRadiationAC::HSBoundaryRadiationAC(const InputParameters & params)
  : ActionComponent(params), HSBoundaryAC(params)
{
  addRequiredTask("add_postprocessor");
}

void
HSBoundaryRadiationAC::addBCs()
{
  const std::string class_name = "ADRadiativeHeatFluxBC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") = temperatureName();
  params.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params.set<MooseFunctorName>("emissivity") = getParam<MooseFunctorName>("emissivity");
  params.set<MooseFunctorName>("T_ambient") = getParam<MooseFunctorName>("T_ambient");
  params.set<MooseFunctorName>("view_factor") = getParam<MooseFunctorName>("view_factor");
  params.set<MooseFunctorName>("scale") = getParam<MooseFunctorName>("scale");
  getProblem().addBoundaryCondition(class_name, genName(name(), "bc"), params);
}

void
HSBoundaryRadiationAC::addPostprocessors()
{
  // Only cylindrical heat structures get a heat rate post-processor, matching classic
  // HSBoundaryRadiation (a plate's integral would be per unit depth)
  if (!cylindricalHS())
    return;

  const std::string class_name = "HeatRateRadiation";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params.set<std::vector<VariableName>>("T") = {temperatureName()};
  params.set<MooseFunctorName>("T_ambient") = getParam<MooseFunctorName>("T_ambient");
  params.set<MooseFunctorName>("emissivity") = getParam<MooseFunctorName>("emissivity");
  params.set<MooseFunctorName>("view_factor") = getParam<MooseFunctorName>("view_factor");
  if (getParam<bool>("scale_heat_rate_pp"))
    params.set<MooseFunctorName>("scale") = getParam<MooseFunctorName>("scale");
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
  getProblem().addPostprocessor(class_name, genSafeName(name(), "integral"), params);
}
