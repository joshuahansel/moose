//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HSBoundaryExternalAppConvectionAC.h"

// Registered under the classic Component's name; see FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryExternalAppConvectionAC,
                           "HSBoundaryExternalAppConvection",
                           "add_aux_variable");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryExternalAppConvectionAC,
                           "HSBoundaryExternalAppConvection",
                           "add_bc");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryExternalAppConvectionAC,
                           "HSBoundaryExternalAppConvection",
                           "add_postprocessor");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryExternalAppConvectionAC,
                           "HSBoundaryExternalAppConvection",
                           "check_integrity");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HSBoundaryExternalAppConvectionAC,
                               "HSBoundaryExternalAppConvection");

InputParameters
HSBoundaryExternalAppConvectionAC::validParams()
{
  InputParameters params = HSBoundaryAC::validParams();
  params.addParam<VariableName>("T_ext", "T_ext", "Temperature from external application");
  params.addParam<VariableName>(
      "htc_ext", "htc_ext", "Heat transfer coefficient from external application");
  params.addParam<FunctionName>("scale", 1.0, "Function by which to scale the boundary condition");
  params.addParam<bool>(
      "scale_heat_rate_pp",
      true,
      "If true, the scaling function is applied to the heat rate post-processor.");
  params.addClassDescription("Heat structure boundary condition to perform convective heat "
                             "transfer with an external application");
  return params;
}

HSBoundaryExternalAppConvectionAC::HSBoundaryExternalAppConvectionAC(const InputParameters & params)
  : ActionComponent(params), HSBoundaryAC(params)
{
  addRequiredTask("add_aux_variable");
  addRequiredTask("add_postprocessor");
}

void
HSBoundaryExternalAppConvectionAC::addAuxVariables()
{
  requestAuxVariable(getParam<VariableName>("T_ext"), "LAGRANGE", "FIRST");
  requestAuxVariable(getParam<VariableName>("htc_ext"), "LAGRANGE", "FIRST");
}

void
HSBoundaryExternalAppConvectionAC::addBCs()
{
  const std::string class_name = "ADExternalAppConvectionHeatTransferBC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") = temperatureName();
  params.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params.set<std::vector<VariableName>>("T_ext") = {getParam<VariableName>("T_ext")};
  params.set<std::vector<VariableName>>("htc_ext") = {getParam<VariableName>("htc_ext")};
  params.set<FunctionName>("scale") = getParam<FunctionName>("scale");
  getProblem().addBoundaryCondition(class_name, genName(name(), "bc"), params);
}

void
HSBoundaryExternalAppConvectionAC::addPostprocessors()
{
  // Only cylindrical heat structures get a heat rate post-processor, matching classic
  // HSBoundaryExternalAppConvection (a plate's integral would be per unit depth)
  if (!cylindricalHS())
    return;

  // HeatRateConvection takes functors, so the transferred variables stand in for its ambient
  // temperature and heat transfer coefficient; this is equivalent to classic THM's
  // HeatRateExternalAppConvectionRZ without its explicit RZ circumference factor (see
  // HSBoundaryAC).
  const std::string class_name = "HeatRateConvection";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params.set<std::vector<VariableName>>("T") = {temperatureName()};
  params.set<MooseFunctorName>("T_ambient") = getParam<VariableName>("T_ext");
  params.set<MooseFunctorName>("htc") = getParam<VariableName>("htc_ext");
  if (getParam<bool>("scale_heat_rate_pp"))
    params.set<MooseFunctorName>("scale") = getParam<FunctionName>("scale");
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
  getProblem().addPostprocessor(class_name, genSafeName(name(), "integral"), params);
}
