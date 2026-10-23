//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SimpleTurbine1PhaseAC.h"
#include "THMUtils.h"

// Registered under the classic Component's name ('SimpleTurbine1Phase'); see FlowChannel1PhaseAC.C
// for why. The tasks are the same as for VolumeJunction1PhaseAC.
#define registerSimpleTurbine1PhaseACTask(task)                                                    \
  registerMooseActionAliased(                                                                      \
      "ThermalHydraulicsApp", SimpleTurbine1PhaseAC, "SimpleTurbine1Phase", task)
registerSimpleTurbine1PhaseACTask("add_mesh_generator");
registerSimpleTurbine1PhaseACTask("add_variable");
registerSimpleTurbine1PhaseACTask("add_aux_variable");
registerSimpleTurbine1PhaseACTask("add_ic");
registerSimpleTurbine1PhaseACTask("add_user_object");
registerSimpleTurbine1PhaseACTask("add_bc");
registerSimpleTurbine1PhaseACTask("add_kernel");
registerSimpleTurbine1PhaseACTask("add_aux_kernel");
registerSimpleTurbine1PhaseACTask("add_material");
registerActionComponentAliased("ThermalHydraulicsApp",
                               SimpleTurbine1PhaseAC,
                               "SimpleTurbine1Phase");

InputParameters
SimpleTurbine1PhaseAC::validParams()
{
  InputParameters params = JunctionParallelChannels1PhaseAC::validParams();
  params.addRequiredParam<Real>("power", "Turbine power [W]");
  params.addRequiredParam<bool>("on", "Flag determining if turbine is operating or not [-]");
  params.declareControllable("power on");
  params.addClassDescription(
      "Simple turbine model that extracts prescribed power from the working fluid");
  return params;
}

SimpleTurbine1PhaseAC::SimpleTurbine1PhaseAC(const InputParameters & params)
  : ActionComponent(params),
    JunctionParallelChannels1PhaseAC(params),
    _on(getParam<bool>("on")),
    _power(getParam<Real>("power")),
    _W_dot_var_name("W_dot")
{
}

void
SimpleTurbine1PhaseAC::setJunctionUserObjectParameters(InputParameters & params) const
{
  params.set<bool>("on") = _on;
  params.set<Real>("W_dot") = _power;
}

void
SimpleTurbine1PhaseAC::connectJunctionUserObjectControllableParams(const std::string & class_name)
{
  THM::connectComponentControllableParam(*this, "power", class_name, _junction_uo_name, "W_dot");
  THM::connectComponentControllableParam(*this, "on", class_name, _junction_uo_name);
}

void
SimpleTurbine1PhaseAC::actOnAdditionalTasks()
{
  JunctionParallelChannels1PhaseAC::actOnAdditionalTasks();

  if (_current_task == "add_aux_variable")
    _coordinator.requestVariable(false, _W_dot_var_name, "MONOMIAL", "CONSTANT", 1.0, _blocks);
  else if (_current_task == "add_aux_kernel")
    addPowerAuxKernel();
}

void
SimpleTurbine1PhaseAC::addPowerAuxKernel()
{
  const std::string class_name = "SimpleTurbinePowerFieldAux";
  const std::string obj_name = genName(name(), "W_dot_aux");
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<AuxVariableName>("variable") = _W_dot_var_name;
  params.set<Real>("value") = _power;
  params.set<bool>("on") = _on;
  params.set<std::vector<SubdomainName>>("block") = _blocks;
  getProblem().addAuxKernel(class_name, obj_name, params);
  THM::connectComponentControllableParam(*this, "power", class_name, obj_name, "value");
  THM::connectComponentControllableParam(*this, "on", class_name, obj_name);
}
