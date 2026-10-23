//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "JunctionParallelChannels1PhaseAC.h"
#include "FlowChannel1PhaseAC.h"

// Registered under the classic Component's name ('JunctionParallelChannels1Phase'); see
// FlowChannel1PhaseAC.C for why. The tasks are the same as for VolumeJunction1PhaseAC.
#define registerJunctionParallelChannels1PhaseACTask(task)                                         \
  registerMooseActionAliased("ThermalHydraulicsApp",                                               \
                             JunctionParallelChannels1PhaseAC,                                     \
                             "JunctionParallelChannels1Phase",                                     \
                             task)
registerJunctionParallelChannels1PhaseACTask("add_mesh_generator");
registerJunctionParallelChannels1PhaseACTask("add_variable");
registerJunctionParallelChannels1PhaseACTask("add_aux_variable");
registerJunctionParallelChannels1PhaseACTask("add_ic");
registerJunctionParallelChannels1PhaseACTask("add_user_object");
registerJunctionParallelChannels1PhaseACTask("add_bc");
registerJunctionParallelChannels1PhaseACTask("add_kernel");
registerJunctionParallelChannels1PhaseACTask("add_aux_kernel");
registerJunctionParallelChannels1PhaseACTask("add_material");
registerActionComponentAliased("ThermalHydraulicsApp",
                               JunctionParallelChannels1PhaseAC,
                               "JunctionParallelChannels1Phase");

InputParameters
JunctionParallelChannels1PhaseAC::validParams()
{
  InputParameters params = VolumeJunction1PhaseAC::validParams();
  params.addClassDescription("Junction between 1-phase flow channels that are parallel");
  return params;
}

JunctionParallelChannels1PhaseAC::JunctionParallelChannels1PhaseAC(const InputParameters & params)
  : ActionComponent(params), VolumeJunction1PhaseAC(params)
{
}

void
JunctionParallelChannels1PhaseAC::buildVolumeJunctionUserObject()
{
  const std::string class_name = junctionUserObjectClass();
  InputParameters params = junctionUserObjectParameters(class_name);
  params.set<std::string>("component_name") = name();
  params.set<RealVectorValue>("dir_c0") = connectedDirection(0);
  params.set<bool>("apply_velocity_scaling") = _apply_velocity_scaling;
  setJunctionUserObjectParameters(params);
  getProblem().addUserObject(class_name, _junction_uo_name, params);
  THM::connectComponentControllableParam(*this, "K", class_name, _junction_uo_name);
  connectJunctionUserObjectControllableParams(class_name);
}
