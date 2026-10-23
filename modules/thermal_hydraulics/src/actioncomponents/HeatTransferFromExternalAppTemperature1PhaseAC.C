//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatTransferFromExternalAppTemperature1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"
#include "FlowModelSinglePhase.h"
#include "THMNames.h"
#include "THMVariableCoordinator.h"

// Registered under the classic Component's name ('HeatTransferFromExternalAppTemperature1Phase');
// see FlowChannel1PhaseAC.C for why.
registerHeatTransferFrom1PhaseACTasks("ThermalHydraulicsApp",
                                      HeatTransferFromExternalAppTemperature1PhaseAC,
                                      "HeatTransferFromExternalAppTemperature1Phase");

InputParameters
HeatTransferFromExternalAppTemperature1PhaseAC::validParams()
{
  InputParameters params = HeatTransferFrom1PhaseAC::validParams();
  params.addParam<FunctionName>("initial_T_wall", "Initial condition for wall temperature [K]");
  params.addParam<VariableName>("T_ext", "Name of the wall temperature variable");
  params.addParam<MooseEnum>("var_type",
                             MooseEnum("nodal elemental", "nodal"),
                             "The type of wall temperature variable (nodal, elemental).");
  params.addClassDescription("Heat transfer into 1-phase flow channel from temperature provided by "
                             "an external application");
  return params;
}

HeatTransferFromExternalAppTemperature1PhaseAC::HeatTransferFromExternalAppTemperature1PhaseAC(
    const InputParameters & params)
  : ActionComponent(params), HeatTransferFrom1PhaseAC(params)
{
  if (!isParamValid("T_ext") && !isParamValid("initial_T_wall") && !_app.isRestarting())
    paramError("initial_T_wall", "Missing initial condition for wall temperature.");
}

std::optional<VariableName>
HeatTransferFromExternalAppTemperature1PhaseAC::specifiedWallTemperatureName() const
{
  if (isParamValid("T_ext"))
    return getParam<VariableName>("T_ext");
  return std::nullopt;
}

void
HeatTransferFromExternalAppTemperature1PhaseAC::addHeatTransferVariables()
{
  // With 'T_ext' the variable is created elsewhere
  if (isParamValid("T_ext"))
    return;

  const bool nodal = getParam<MooseEnum>("var_type") == "nodal";
  THMVariableCoordinator::findOrCreate(_awh, _action_factory)
      .requestVariable(false,
                       wallTemperatureName(),
                       nodal ? "LAGRANGE" : "MONOMIAL",
                       nodal ? "FIRST" : "CONSTANT",
                       1.0,
                       closuresBlocks());
}

void
HeatTransferFromExternalAppTemperature1PhaseAC::addHeatTransferICs()
{
  if (isParamValid("T_ext") || !isParamValid("initial_T_wall"))
    return;

  const std::string class_name = "FunctionIC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<VariableName>("variable") = wallTemperatureName();
  params.set<std::vector<SubdomainName>>("block") = closuresBlocks();
  params.set<FunctionName>("function") = getParam<FunctionName>("initial_T_wall");
  getProblem().addInitialCondition(class_name, genName(name(), "T_wall_ic"), params);
}

void
HeatTransferFromExternalAppTemperature1PhaseAC::addHeatTransferKernels()
{
  const std::string class_name = "ADOneDEnergyWallHeating";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") = THM::RHOEA;
  params.set<std::vector<SubdomainName>>("block") = closuresBlocks();
  params.set<std::vector<VariableName>>("T_wall") = {wallTemperatureName()};
  params.set<MaterialPropertyName>("Hw") = wallHeatTransferCoefficientName();
  params.set<std::vector<VariableName>>("P_hf") = {heatedPerimeterName()};
  params.set<MaterialPropertyName>("T") = FlowModelSinglePhase::TEMPERATURE;
  getProblem().addKernel(class_name, genName(name(), "wall_heat_transfer"), params);
}
