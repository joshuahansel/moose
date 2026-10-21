//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatTransferFrom1PhaseAC.h"
#include "THMUtils.h"
#include "FlowChannel1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"
#include "ClosuresBase.h"
#include "ClosuresRegistry.h"
#include "THMNames.h"
#include "MooseUtils.h"

InputParameters
HeatTransferFrom1PhaseAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addRequiredParam<std::string>("flow_channel",
                                       "Name of the flow channel component to connect to");
  params.addParam<FunctionName>("Hw", "Convective heat transfer coefficient [W/(m^2-K)]");
  params.addParam<FunctionName>("P_hf",
                                "Heat flux perimeter [m]. If not provided, it is computed from the "
                                "flow channel area assuming a circular cross section.");
  return params;
}

HeatTransferFrom1PhaseAC::HeatTransferFrom1PhaseAC(const InputParameters & params)
  : ActionComponent(params),
    _P_hf_fn_name(isParamValid("P_hf") ? getParam<FunctionName>("P_hf")
                                       : FunctionName(genName(name(), "P_hf_fn"))),
    _Hw_name(THM::HEAT_TRANSFER_COEFFICIENT_WALL)
{
  THM::informExperimental(*this);
  addRequiredTask("THM:add_heat_transfers");
  addRequiredTask("add_function");
  addRequiredTask("add_aux_variable");
  addRequiredTask("add_ics_physics");
  addRequiredTask("add_material");
  addRequiredTask("add_kernel");
  addRequiredTask("add_aux_kernel");
  addRequiredTask("check_integrity");
}

const FlowChannel1PhasePhysics &
HeatTransferFrom1PhaseAC::connectedPhysics() const
{
  return _awh.getAction<FlowChannel1PhaseAC>(getParam<std::string>("flow_channel")).getPhysics();
}

void
HeatTransferFrom1PhaseAC::actOnAdditionalTasks()
{
  if (_current_task == "THM:add_heat_transfers")
    connectedPhysics().addHeatTransferName(name(), isTemperatureType());
  else if (_current_task == "add_function")
  {
    // Build the default (circular) heated perimeter function only if 'P_hf' was not provided;
    // mirrors HeatTransferBase::addHeatedPerimeter().
    if (!isParamValid("P_hf"))
    {
      const std::string class_name = "GeneralizedCircumference";
      InputParameters params = getFactory().getValidParams(class_name);
      params.set<FunctionName>("area_function") = connectedPhysics().areaFunctionName();
      getProblem().addFunction(class_name, _P_hf_fn_name, params);
    }
  }
  else if (_current_task == "add_aux_variable")
  {
    const std::string class_name = "MooseVariable";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
    params.set<MooseEnum>("family") = "MONOMIAL";
    params.set<MooseEnum>("order") = "CONSTANT";
    getProblem().addAuxVariable(class_name, THM::HEAT_FLUX_PERIMETER, params);

    addHeatTransferVariables();
  }
  else if (_current_task == "add_ics_physics")
  {
    const std::string class_name = "FunctionIC";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<VariableName>("variable") = THM::HEAT_FLUX_PERIMETER;
    params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
    params.set<FunctionName>("function") = _P_hf_fn_name;
    getProblem().addInitialCondition(class_name, genName(name(), "P_hf_ic"), params);

    addHeatTransferICs();
  }
  else if (_current_task == "add_kernel")
    addHeatTransferKernels();
  else if (_current_task == "add_aux_kernel")
  {
    const std::string class_name = "FunctionAux";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<AuxVariableName>("variable") = THM::HEAT_FLUX_PERIMETER;
    params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
    params.set<FunctionName>("function") = _P_hf_fn_name;
    ExecFlagEnum execute_on(MooseUtils::getDefaultExecFlagEnum());
    execute_on = {EXEC_TIMESTEP_BEGIN, EXEC_INITIAL};
    params.set<ExecFlagEnum>("execute_on") = execute_on;
    getProblem().addAuxKernel(class_name, genName(name(), "P_hf_auxkernel"), params);

    addHeatTransferAuxKernels();
  }
}

void
HeatTransferFrom1PhaseAC::addMaterials()
{
  const FlowChannel1PhasePhysics & physics = connectedPhysics();

  // The wall heat transfer coefficient material is supplied by the flow channel's closures objects,
  // reused here exactly as classic HeatTransferBase reuses its flow channel's closures.
  for (const auto & closures : physics.getClosuresObjects())
    closures->addMooseObjectsHeatTransfer(*this, physics);

  addHeatTransferMaterials();
}

void
HeatTransferFrom1PhaseAC::checkIntegrity()
{
  const FlowChannel1PhasePhysics & physics = connectedPhysics();

  ClosuresRegistry & registry = ClosuresRegistry::findOrCreate(_awh, _action_factory);
  for (const auto & closures : physics.getClosuresObjects())
    closures->checkHeatTransfer(*this, physics);
  registry.emitLoggedMessages();
}
