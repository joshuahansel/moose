//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ShaftConnectedPump1PhaseAC.h"
#include "THMUtils.h"
#include "Numerics.h"

// Registered under the classic Component's name ('ShaftConnectedPump1Phase'); see
// FlowChannel1PhaseAC.C for why. The tasks are the same as for VolumeJunction1PhaseAC.
#define registerShaftConnectedPump1PhaseACTask(task)                                               \
  registerMooseActionAliased(                                                                      \
      "ThermalHydraulicsApp", ShaftConnectedPump1PhaseAC, "ShaftConnectedPump1Phase", task)
registerShaftConnectedPump1PhaseACTask("add_mesh_generator");
registerShaftConnectedPump1PhaseACTask("add_variable");
registerShaftConnectedPump1PhaseACTask("add_aux_variable");
registerShaftConnectedPump1PhaseACTask("add_ic");
registerShaftConnectedPump1PhaseACTask("add_user_object");
registerShaftConnectedPump1PhaseACTask("add_bc");
registerShaftConnectedPump1PhaseACTask("add_kernel");
registerShaftConnectedPump1PhaseACTask("add_aux_kernel");
registerShaftConnectedPump1PhaseACTask("add_material");
registerActionComponentAliased("ThermalHydraulicsApp",
                               ShaftConnectedPump1PhaseAC,
                               "ShaftConnectedPump1Phase");

InputParameters
ShaftConnectedPump1PhaseAC::validParams()
{
  InputParameters params = ShaftConnectedJunction1PhaseAC::validParams();
  params.addRequiredParam<Real>("volumetric_rated", "Rated pump volumetric flow rate [m^3/s]");
  params.addRequiredParam<Real>("head_rated", "Rated pump head [m]");
  params.addRequiredParam<Real>("torque_rated", "Rated pump torque [N-m]");
  params.addRequiredParam<Real>("density_rated", "Rated pump fluid density [kg/m^3]");
  params.addRequiredParam<FunctionName>("head", "Function to compute data for pump head [-]");
  params.addRequiredParam<FunctionName>("torque_hydraulic",
                                        "Function to compute data for pump torque [-]");
  params.addParam<Real>("transition_width",
                        1e-3,
                        "Transition width for sign of the frictional torque at 0 speed over rated "
                        "speed ratio.");
  params.addClassDescription(
      "Pump between two 1-phase flow channels whose speed is given by a connected shaft");
  return params;
}

ShaftConnectedPump1PhaseAC::ShaftConnectedPump1PhaseAC(const InputParameters & params)
  : ActionComponent(params), ShaftConnectedJunction1PhaseAC(params)
{
}

void
ShaftConnectedPump1PhaseAC::buildVolumeJunctionUserObject()
{
  const std::string class_name = "ADShaftConnectedPump1PhaseUserObject";
  InputParameters params = shaftConnectedUserObjectParameters(class_name);
  params.set<Real>("gravity_magnitude") = THM::gravity_const;
  for (const auto & param :
       {"volumetric_rated", "head_rated", "torque_rated", "density_rated", "transition_width"})
    params.set<Real>(param) = getParam<Real>(param);
  for (const auto & param : {"head", "torque_hydraulic"})
    params.set<FunctionName>(param) = getParam<FunctionName>(param);
  params.set<std::string>("pump_name") = name();
  getProblem().addUserObject(class_name, _junction_uo_name, params);
  THM::connectComponentControllableParam(*this, "K", class_name, _junction_uo_name);
}

std::vector<std::pair<VariableName, Real>>
ShaftConnectedPump1PhaseAC::shaftConnectedVariables() const
{
  return {{"head", 0.0}, {"hydraulic_torque", 0.0}};
}

std::vector<std::pair<std::string, VariableName>>
ShaftConnectedPump1PhaseAC::auxKernelQuantities() const
{
  return {{"pump_head", "head"},
          {"hydraulic_torque", "hydraulic_torque"},
          {"friction_torque", _friction_torque_var_name},
          {"moment_of_inertia", _moment_of_inertia_var_name}};
}
