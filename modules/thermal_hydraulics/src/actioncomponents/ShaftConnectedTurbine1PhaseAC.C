//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ShaftConnectedTurbine1PhaseAC.h"
#include "THMUtils.h"

// Registered under the classic Component's name ('ShaftConnectedTurbine1Phase'); see
// FlowChannel1PhaseAC.C for why. The tasks are the same as for VolumeJunction1PhaseAC.
#define registerShaftConnectedTurbine1PhaseACTask(task)                                            \
  registerMooseActionAliased(                                                                      \
      "ThermalHydraulicsApp", ShaftConnectedTurbine1PhaseAC, "ShaftConnectedTurbine1Phase", task)
registerShaftConnectedTurbine1PhaseACTask("add_mesh_generator");
registerShaftConnectedTurbine1PhaseACTask("add_variable");
registerShaftConnectedTurbine1PhaseACTask("add_aux_variable");
registerShaftConnectedTurbine1PhaseACTask("add_ic");
registerShaftConnectedTurbine1PhaseACTask("add_user_object");
registerShaftConnectedTurbine1PhaseACTask("add_bc");
registerShaftConnectedTurbine1PhaseACTask("add_kernel");
registerShaftConnectedTurbine1PhaseACTask("add_aux_kernel");
registerShaftConnectedTurbine1PhaseACTask("add_material");
registerActionComponentAliased("ThermalHydraulicsApp",
                               ShaftConnectedTurbine1PhaseAC,
                               "ShaftConnectedTurbine1Phase");

InputParameters
ShaftConnectedTurbine1PhaseAC::validParams()
{
  InputParameters params = ShaftConnectedJunction1PhaseAC::validParams();
  params.addRequiredParam<Real>("D_wheel",
                                "Any physical dimension of the turbine as a measure of the "
                                "turbine's size, usually the rotor diameter [m]");
  params.addRequiredParam<FunctionName>("head_coefficient",
                                        "Head coefficient vs flow coefficient function [-]");
  params.addRequiredParam<FunctionName>("power_coefficient",
                                        "Power coefficient vs flow coefficient function [-]");
  params.addClassDescription(
      "Turbine between two 1-phase flow channels whose speed is given by a connected shaft");
  return params;
}

ShaftConnectedTurbine1PhaseAC::ShaftConnectedTurbine1PhaseAC(const InputParameters & params)
  : ActionComponent(params), ShaftConnectedJunction1PhaseAC(params)
{
}

void
ShaftConnectedTurbine1PhaseAC::buildVolumeJunctionUserObject()
{
  const std::string class_name = "ADShaftConnectedTurbine1PhaseUserObject";
  InputParameters params = shaftConnectedUserObjectParameters(class_name);
  params.set<Real>("D_wheel") = getParam<Real>("D_wheel");
  for (const auto & param : {"head_coefficient", "power_coefficient"})
    params.set<FunctionName>(param) = getParam<FunctionName>(param);
  params.set<std::string>("turbine_name") = name();
  getProblem().addUserObject(class_name, _junction_uo_name, params);
  THM::connectComponentControllableParam(*this, "K", class_name, _junction_uo_name);
}

std::vector<std::pair<VariableName, Real>>
ShaftConnectedTurbine1PhaseAC::shaftConnectedVariables() const
{
  return {{"delta_p", 0.0}, {"power", 0.0}, {"driving_torque", 0.0}, {"flow_coeff", 0.0}};
}

std::vector<std::pair<std::string, VariableName>>
ShaftConnectedTurbine1PhaseAC::auxKernelQuantities() const
{
  return {{"delta_p", "delta_p"},
          {"flow_coefficient", "flow_coeff"},
          {"driving_torque", "driving_torque"},
          {"friction_torque", _friction_torque_var_name},
          {"moment_of_inertia", _moment_of_inertia_var_name},
          {"power", "power"}};
}
