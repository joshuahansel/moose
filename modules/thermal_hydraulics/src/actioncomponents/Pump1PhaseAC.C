//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Pump1PhaseAC.h"
#include "THMUtils.h"
#include "Numerics.h"

// Registered under the classic Component's name ('Pump1Phase'); see FlowChannel1PhaseAC.C for why.
// The tasks are the same as for VolumeJunction1PhaseAC.
#define registerPump1PhaseACTask(task)                                                             \
  registerMooseActionAliased("ThermalHydraulicsApp", Pump1PhaseAC, "Pump1Phase", task)
registerPump1PhaseACTask("add_mesh_generator");
registerPump1PhaseACTask("add_variable");
registerPump1PhaseACTask("add_aux_variable");
registerPump1PhaseACTask("add_ic");
registerPump1PhaseACTask("add_user_object");
registerPump1PhaseACTask("add_bc");
registerPump1PhaseACTask("add_kernel");
registerPump1PhaseACTask("add_aux_kernel");
registerPump1PhaseACTask("add_material");
registerActionComponentAliased("ThermalHydraulicsApp", Pump1PhaseAC, "Pump1Phase");

InputParameters
Pump1PhaseAC::validParams()
{
  InputParameters params = VolumeJunction1PhaseAC::validParams();
  params.addRequiredParam<Real>("head", "Pump head [m]");
  params.makeParamRequired<Real>("A_ref");
  params.declareControllable("head");
  params.addClassDescription("Pump between two 1-phase flow channels that has a non-zero volume");
  return params;
}

Pump1PhaseAC::Pump1PhaseAC(const InputParameters & params)
  : ActionComponent(params), VolumeJunction1PhaseAC(params), _head(getParam<Real>("head"))
{
}

void
Pump1PhaseAC::buildVolumeJunctionUserObject()
{
  const std::string class_name = "ADPump1PhaseUserObject";
  InputParameters params = junctionUserObjectParameters(class_name);
  params.set<Real>("head") = _head;
  params.set<Real>("gravity_magnitude") = THM::gravity_const;
  params.set<bool>("apply_velocity_scaling") = _apply_velocity_scaling;
  getProblem().addUserObject(class_name, _junction_uo_name, params);
  THM::connectComponentControllableParam(*this, "head", class_name, _junction_uo_name);
  THM::connectComponentControllableParam(*this, "K", class_name, _junction_uo_name);
}
