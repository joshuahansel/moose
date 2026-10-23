//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatTransferFromExternalAppHeatFlux1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"
#include "FlowModelSinglePhase.h"
#include "THMNames.h"
#include "THMVariableCoordinator.h"

// Registered under the classic Component's name ('HeatTransferFromExternalAppHeatFlux1Phase'); see
// FlowChannel1PhaseAC.C for why.
registerHeatTransferFrom1PhaseACTasks("ThermalHydraulicsApp",
                                      HeatTransferFromExternalAppHeatFlux1PhaseAC,
                                      "HeatTransferFromExternalAppHeatFlux1Phase");

InputParameters
HeatTransferFromExternalAppHeatFlux1PhaseAC::validParams()
{
  InputParameters params = HeatTransferFrom1PhaseAC::validParams();
  params.addClassDescription("Heat transfer specified by heat flux provided by an external "
                             "application going into 1-phase flow channel.");
  return params;
}

HeatTransferFromExternalAppHeatFlux1PhaseAC::HeatTransferFromExternalAppHeatFlux1PhaseAC(
    const InputParameters & params)
  : ActionComponent(params), HeatTransferFrom1PhaseAC(params)
{
}

void
HeatTransferFromExternalAppHeatFlux1PhaseAC::addHeatTransferVariables()
{
  THMVariableCoordinator::findOrCreate(_awh, _action_factory)
      .requestVariable(false, wallHeatFluxName(), "MONOMIAL", "CONSTANT", 1.0, closuresBlocks());
}

void
HeatTransferFromExternalAppHeatFlux1PhaseAC::addHeatTransferICs()
{
  const std::string class_name = "ConstantIC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<VariableName>("variable") = wallHeatFluxName();
  params.set<std::vector<SubdomainName>>("block") = closuresBlocks();
  params.set<Real>("value") = 0;
  getProblem().addInitialCondition(class_name, genName(name(), "q_wall_ic"), params);
}

void
HeatTransferFromExternalAppHeatFlux1PhaseAC::addHeatTransferMaterials()
{
  const std::string class_name = "ADCoupledVariableValueMaterial";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<SubdomainName>>("block") = closuresBlocks();
  params.set<MaterialPropertyName>("prop_name") = wallHeatFluxName();
  params.set<std::vector<VariableName>>("coupled_variable") = {wallHeatFluxName()};
  getProblem().addMaterial(class_name, genName(name(), "q_wall_material"), params);
}

void
HeatTransferFromExternalAppHeatFlux1PhaseAC::addHeatTransferKernels()
{
  const std::string class_name = "ADOneDEnergyWallHeatFlux";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") = THM::RHOEA;
  params.set<std::vector<SubdomainName>>("block") = closuresBlocks();
  params.set<MaterialPropertyName>("q_wall") = wallHeatFluxName();
  params.set<std::vector<VariableName>>("P_hf") = {heatedPerimeterName()};
  getProblem().addKernel(class_name, genName(name(), "wall_heat"), params);
}
