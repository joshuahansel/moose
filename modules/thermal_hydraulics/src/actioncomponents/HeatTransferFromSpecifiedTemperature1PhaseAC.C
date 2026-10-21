//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatTransferFromSpecifiedTemperature1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"
#include "FlowModelSinglePhase.h"
#include "THMNames.h"
#include "THMVariableCoordinator.h"
#include "MooseUtils.h"

// Registered under the classic Component's name ('HeatTransferFromSpecifiedTemperature1Phase'); see
// FlowChannel1PhaseAC.C for why.
registerHeatTransferFrom1PhaseACTasks("ThermalHydraulicsApp",
                                      HeatTransferFromSpecifiedTemperature1PhaseAC,
                                      "HeatTransferFromSpecifiedTemperature1Phase");

InputParameters
HeatTransferFromSpecifiedTemperature1PhaseAC::validParams()
{
  InputParameters params = HeatTransferFrom1PhaseAC::validParams();
  params.addRequiredParam<FunctionName>("T_wall", "Specified wall temperature [K]");
  params.addParam<MooseEnum>("var_type",
                             MooseEnum("nodal elemental", "nodal"),
                             "The type of wall temperature variable (nodal, elemental).");
  params.addClassDescription("Applies a specified wall temperature to a FlowChannel1PhaseAC.");
  return params;
}

HeatTransferFromSpecifiedTemperature1PhaseAC::HeatTransferFromSpecifiedTemperature1PhaseAC(
    const InputParameters & params)
  : ActionComponent(params), HeatTransferFrom1PhaseAC(params)
{
}

void
HeatTransferFromSpecifiedTemperature1PhaseAC::addHeatTransferVariables()
{
  const bool nodal = getParam<MooseEnum>("var_type") == "nodal";

  THMVariableCoordinator::findOrCreate(_awh, _action_factory)
      .requestVariable(false,
                       THM::TEMPERATURE_WALL,
                       nodal ? "LAGRANGE" : "MONOMIAL",
                       nodal ? "FIRST" : "CONSTANT",
                       1.0,
                       connectedPhysics().getClosuresBlocks());
}

void
HeatTransferFromSpecifiedTemperature1PhaseAC::addHeatTransferICs()
{
  const std::string class_name = "FunctionIC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<VariableName>("variable") = THM::TEMPERATURE_WALL;
  params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
  params.set<FunctionName>("function") = getParam<FunctionName>("T_wall");
  getProblem().addInitialCondition(class_name, genName(name(), "T_wall_ic"), params);
}

void
HeatTransferFromSpecifiedTemperature1PhaseAC::addHeatTransferAuxKernels()
{
  const std::string class_name = "FunctionAux";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<AuxVariableName>("variable") = THM::TEMPERATURE_WALL;
  params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
  params.set<FunctionName>("function") = getParam<FunctionName>("T_wall");
  ExecFlagEnum execute_on(MooseUtils::getDefaultExecFlagEnum());
  execute_on = {EXEC_INITIAL, EXEC_LINEAR};
  params.set<ExecFlagEnum>("execute_on") = execute_on;
  getProblem().addAuxKernel(class_name, genName(name(), "T_wall_auxkernel"), params);
}

void
HeatTransferFromSpecifiedTemperature1PhaseAC::addHeatTransferKernels()
{
  const std::string class_name = "ADOneDEnergyWallHeating";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") = THM::RHOEA;
  params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
  params.set<std::vector<VariableName>>("T_wall") = {THM::TEMPERATURE_WALL};
  params.set<MaterialPropertyName>("Hw") = THM::HEAT_TRANSFER_COEFFICIENT_WALL;
  params.set<std::vector<VariableName>>("P_hf") = {THM::HEAT_FLUX_PERIMETER};
  params.set<MaterialPropertyName>("T") = FlowModelSinglePhase::TEMPERATURE;
  getProblem().addKernel(class_name, genName(name(), "wall_heat_transfer"), params);
}
