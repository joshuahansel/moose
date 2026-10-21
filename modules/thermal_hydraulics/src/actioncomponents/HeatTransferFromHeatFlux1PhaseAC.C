//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatTransferFromHeatFlux1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"
#include "THMNames.h"

// Registered under the classic Component's name ('HeatTransferFromHeatFlux1Phase'); see
// FlowChannel1PhaseAC.C for why.
registerHeatTransferFrom1PhaseACTasks("ThermalHydraulicsApp",
                                      HeatTransferFromHeatFlux1PhaseAC,
                                      "HeatTransferFromHeatFlux1Phase");

InputParameters
HeatTransferFromHeatFlux1PhaseAC::validParams()
{
  InputParameters params = HeatTransferFrom1PhaseAC::validParams();
  params.addRequiredParam<FunctionName>("q_wall", "Specified wall heat flux [W/m^2]");
  params.addClassDescription("Applies a specified wall heat flux to a FlowChannel1PhaseAC.");
  return params;
}

HeatTransferFromHeatFlux1PhaseAC::HeatTransferFromHeatFlux1PhaseAC(const InputParameters & params)
  : ActionComponent(params), HeatTransferFrom1PhaseAC(params)
{
}

void
HeatTransferFromHeatFlux1PhaseAC::addHeatTransferMaterials()
{
  const std::string class_name = "ADGenericFunctionMaterial";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
  params.set<std::vector<std::string>>("prop_names") = {THM::HEAT_FLUX_WALL};
  params.set<std::vector<FunctionName>>("prop_values") = {getParam<FunctionName>("q_wall")};
  getProblem().addMaterial(class_name, genName(name(), "q_wall_material"), params);
}

void
HeatTransferFromHeatFlux1PhaseAC::addHeatTransferKernels()
{
  const std::string class_name = "ADOneDEnergyWallHeatFlux";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") = THM::RHOEA;
  params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
  params.set<MaterialPropertyName>("q_wall") = THM::HEAT_FLUX_WALL;
  params.set<std::vector<VariableName>>("P_hf") = {THM::HEAT_FLUX_PERIMETER};
  getProblem().addKernel(class_name, genName(name(), "wall_heat"), params);
}
