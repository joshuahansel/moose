//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatSourceVolumetric1PhaseAC.h"
#include "FlowChannel1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"
#include "THMUtils.h"
#include "THMNames.h"

// Registered under the classic Component's name ('HeatSourceVolumetric1Phase'); see
// FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatSourceVolumetric1PhaseAC,
                           "HeatSourceVolumetric1Phase",
                           "add_kernel");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HeatSourceVolumetric1PhaseAC,
                               "HeatSourceVolumetric1Phase");

InputParameters
HeatSourceVolumetric1PhaseAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addRequiredParam<std::string>("flow_channel",
                                       "Flow channel name in which to apply heat source");
  params.addRequiredParam<MooseFunctorName>("q", "Volumetric heat source functor [W/m^3]");
  params.addClassDescription("Volumetric heat source applied on a flow channel");
  return params;
}

HeatSourceVolumetric1PhaseAC::HeatSourceVolumetric1PhaseAC(const InputParameters & params)
  : ActionComponent(params)
{
  THM::informExperimental(*this);
  addRequiredTask("add_kernel");
}

void
HeatSourceVolumetric1PhaseAC::actOnAdditionalTasks()
{
  if (_current_task != "add_kernel")
    return;

  // Look the flow channel up by type among all components, so a wrong name or type is reported
  // as in the classic component
  const auto & fch_name = getParam<std::string>("flow_channel");
  const FlowChannel1PhaseAC * fch = nullptr;
  for (const auto & action : _awh.getActionListByName("list_component"))
    if (action->name() == fch_name)
    {
      fch = dynamic_cast<const FlowChannel1PhaseAC *>(action);
      if (!fch)
        paramError(
            "flow_channel", "The component '", fch_name, "' is not of type 'FlowChannel1Phase'");
    }
  if (!fch)
    paramError("flow_channel", "The component '", fch_name, "' does not exist");

  const std::string class_name = "ADOneD3EqnEnergyHeatSource";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") = THM::RHOEA;
  params.set<std::vector<SubdomainName>>("block") = fch->getPhysics().getClosuresBlocks();
  params.set<MooseFunctorName>("q") = getParam<MooseFunctorName>("q");
  params.set<std::vector<VariableName>>("A") = {THM::AREA};
  getProblem().addKernel(class_name, genName(name(), "rhoE_heat_source"), params);
}
