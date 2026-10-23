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

RealVectorValue
JunctionParallelChannels1PhaseAC::firstDirection() const
{
  const auto & fch_name = _connected_component_names[0];
  const auto & fch = _awh.getAction<FlowChannel1PhaseAC>(fch_name);
  // The flow channel is generated along the x axis if it has no orientation
  return fch.isParamValid("orientation") ? fch.direction() : RealVectorValue(1, 0, 0);
}

void
JunctionParallelChannels1PhaseAC::buildVolumeJunctionUserObject()
{
  const std::string class_name = junctionUserObjectClass();
  InputParameters params = junctionUserObjectParameters(class_name);
  params.set<std::string>("component_name") = name();
  params.set<RealVectorValue>("dir_c0") = firstDirection();
  params.set<bool>("apply_velocity_scaling") = _apply_velocity_scaling;
  setJunctionUserObjectParameters(params);
  getProblem().addUserObject(class_name, _junction_uo_name, params);
  THM::connectComponentControllableParam(*this, "K", class_name, _junction_uo_name);
  connectJunctionUserObjectControllableParams(class_name);
}
