//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FormLoss1PhaseBaseAC.h"
#include "FlowChannel1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"
#include "THMNames.h"
#include "THMUtils.h"

InputParameters
FormLoss1PhaseBaseAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addRequiredParam<std::string>("flow_channel",
                                       "Flow channel where form loss will be applied");
  return params;
}

FormLoss1PhaseBaseAC::FormLoss1PhaseBaseAC(const InputParameters & params) : ActionComponent(params)
{
  THM::informExperimental(*this);
  addRequiredTask("add_aux_variable");
  addRequiredTask("add_material");
  addRequiredTask("add_kernel");
}

const FlowChannel1PhasePhysics &
FormLoss1PhaseBaseAC::connectedPhysics() const
{
  const auto & fch_name = getParam<std::string>("flow_channel");
  for (const auto & action : _awh.getActionListByName("list_component"))
    if (action->name() == fch_name)
    {
      if (const auto * fch = dynamic_cast<const FlowChannel1PhaseAC *>(action))
        return fch->getPhysics();
      paramError(
          "flow_channel", "The component '", fch_name, "' is not of type 'FlowChannel1Phase'");
    }
  paramError("flow_channel", "The component '", fch_name, "' does not exist");
}

void
FormLoss1PhaseBaseAC::addMaterials()
{
  addFormLossMaterials();
}

void
FormLoss1PhaseBaseAC::actOnAdditionalTasks()
{
  if (_current_task == "add_aux_variable")
    addFormLossVariables();
  else if (_current_task == "add_kernel")
  {
    const std::string class_name = "ADOneD3EqnMomentumFormLoss";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<NonlinearVariableName>("variable") = THM::RHOUA;
    params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
    params.set<std::vector<VariableName>>("A") = {THM::AREA};
    params.set<MaterialPropertyName>("rho") = THM::DENSITY;
    params.set<MaterialPropertyName>("vel") = THM::VELOCITY;
    getProblem().addKernel(class_name, genName(name(), THM::RHOUA, "form_loss"), params);
  }
}
