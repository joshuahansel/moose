//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HSBoundaryAmbientConvectionAC.h"

// Registered under the classic Component's name; see FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryAmbientConvectionAC,
                           "HSBoundaryAmbientConvection",
                           "add_bc");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryAmbientConvectionAC,
                           "HSBoundaryAmbientConvection",
                           "add_postprocessor");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryAmbientConvectionAC,
                           "HSBoundaryAmbientConvection",
                           "check_integrity");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HSBoundaryAmbientConvectionAC,
                               "HSBoundaryAmbientConvection");

InputParameters
HSBoundaryAmbientConvectionAC::validParams()
{
  InputParameters params = HSBoundaryAC::validParams();
  params.addRequiredParam<MooseFunctorName>(
      "htc_ambient", "Ambient Convective heat transfer coefficient functor [W/(m^2-K)]");
  params.addRequiredParam<MooseFunctorName>("T_ambient", "Ambient temperature functor [K]");
  params.addParam<MooseFunctorName>(
      "scale", 1.0, "Functor by which to scale the boundary condition");
  params.addParam<bool>(
      "scale_heat_rate_pp",
      true,
      "If true, the 'scale' parameter is also applied to the heat rate post-processor.");
  params.addClassDescription("Applies a convective boundary condition to a heat structure");
  return params;
}

HSBoundaryAmbientConvectionAC::HSBoundaryAmbientConvectionAC(const InputParameters & params)
  : ActionComponent(params), HSBoundaryAC(params)
{
  addRequiredTask("add_postprocessor");
}

void
HSBoundaryAmbientConvectionAC::addBCs()
{
  const std::string class_name = "ADConvectionHeatTransferBC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") = temperatureName();
  params.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params.set<MooseFunctorName>("T_ambient") = getParam<MooseFunctorName>("T_ambient");
  params.set<MooseFunctorName>("htc_ambient") = getParam<MooseFunctorName>("htc_ambient");
  params.set<MooseFunctorName>("scale") = getParam<MooseFunctorName>("scale");
  getProblem().addBoundaryCondition(class_name, genName(name(), "bc"), params);
}

void
HSBoundaryAmbientConvectionAC::addPostprocessors()
{
  // Only cylindrical heat structures get a heat rate post-processor, matching classic
  // HSBoundaryAmbientConvection (a plate's integral would be per unit depth)
  if (!cylindricalHS())
    return;

  const std::string class_name = "HeatRateConvection";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params.set<std::vector<VariableName>>("T") = {temperatureName()};
  params.set<MooseFunctorName>("T_ambient") = getParam<MooseFunctorName>("T_ambient");
  params.set<MooseFunctorName>("htc") = getParam<MooseFunctorName>("htc_ambient");
  if (getParam<bool>("scale_heat_rate_pp"))
    params.set<MooseFunctorName>("scale") = getParam<MooseFunctorName>("scale");
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
  getProblem().addPostprocessor(class_name, genSafeName(name(), "integral"), params);
}
