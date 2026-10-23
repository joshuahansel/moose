//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ShaftConnectedMotorAC.h"
#include "ShaftAC.h"
#include "THMUtils.h"

// Registered under the classic Component's name ('ShaftConnectedMotor'); see FlowChannel1PhaseAC.C
// for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           ShaftConnectedMotorAC,
                           "ShaftConnectedMotor",
                           "add_mesh_generator");
registerMooseActionAliased("ThermalHydraulicsApp",
                           ShaftConnectedMotorAC,
                           "ShaftConnectedMotor",
                           "add_user_object");
registerActionComponentAliased("ThermalHydraulicsApp",
                               ShaftConnectedMotorAC,
                               "ShaftConnectedMotor");

InputParameters
ShaftConnectedMotorAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addRequiredParam<FunctionName>("torque", "Driving torque supplied by the motor [kg-m^2]");
  params.addRequiredParam<FunctionName>("inertia", "Moment of inertia from the motor [N-m]");
  params.addParam<bool>("ad", true, "Use AD version or not");
  params.addClassDescription("Motor to drive a shaft component");
  return params;
}

ShaftConnectedMotorAC::ShaftConnectedMotorAC(const InputParameters & params)
  : ActionComponent(params), _uo_name(genName(name(), "shaftconnected_uo"))
{
  THM::informExperimental(*this);
  addRequiredTask("add_mesh_generator");
  addRequiredTask("add_user_object");
}

void
ShaftConnectedMotorAC::addMeshGenerators()
{
  // The motor has no mesh. This is the first task on which every component exists, so the shaft has
  // already told this motor that it is connected, if it is, and an unconnected motor is reported
  // before the other tasks of the simulation can fail for another reason.
  if (!connectedToShaft())
    paramError("ad", "This component must be connected to a shaft.");
}

void
ShaftConnectedMotorAC::addUserObjects()
{
  const auto & shaft = _awh.getAction<ShaftAC>(_shaft_name);

  const std::string class_name =
      getParam<bool>("ad") ? "ADShaftConnectedMotorUserObject" : "ShaftConnectedMotorUserObject";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<FunctionName>("torque") = getParam<FunctionName>("torque");
  params.set<FunctionName>("inertia") = getParam<FunctionName>("inertia");
  params.set<std::vector<VariableName>>("shaft_speed") = {shaft.omegaVariableName()};
  getProblem().addUserObject(class_name, _uo_name, params);
}
