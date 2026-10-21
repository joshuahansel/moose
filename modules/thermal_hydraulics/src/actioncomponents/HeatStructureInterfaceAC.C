//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatStructureInterfaceAC.h"
#include "ActionFactory.h"

InputParameters
HeatStructureInterfaceAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addRequiredParam<FunctionName>("initial_T", "Initial temperature [K]");
  params.addParam<Real>(
      "scaling_factor_temperature", 1.0, "Scaling factor for solid temperature variable.");
  return params;
}

HeatStructureInterfaceAC::HeatStructureInterfaceAC(const InputParameters & params)
  : ActionComponent(params), _physics(nullptr)
{
}

void
HeatStructureInterfaceAC::createPhysics(const InputParameters & params, bool second_order)
{
  const std::string physics_type = "HeatStructurePhysics";
  InputParameters physics_params = _action_factory.getValidParams(physics_type);
  physics_params.blockFullpath() = params.blockFullpath();
  physics_params.set<bool>("_built_by_moose") = true;
  physics_params.set<std::string>("registered_identifier") = "(AutoBuilt)";

  physics_params.set<VariableName>("temperature_name") = "T_solid";
  physics_params.set<FunctionName>("initial_temperature") = getParam<FunctionName>("initial_T");
  physics_params.set<bool>("second_order") = second_order;
  physics_params.set<Real>("scaling_factor_temperature") =
      getParam<Real>("scaling_factor_temperature");
  // HeatConductionPhysicsBase::addPreconditioning() otherwise unconditionally injects '-pc_type
  // hypre -pc_hypre_type boomeramg' into the PETSc options, silently overriding whatever the
  // simulation's own [Preconditioning] block requests. Deferring keeps preconditioning entirely
  // under the user's/executioner's control, matching how FlowChannel1PhasePhysics never touches
  // preconditioning either.
  physics_params.set<MooseEnum>("preconditioning") = "defer";

  auto physics_action = _action_factory.create(physics_type, name() + "_physics", physics_params);
  _physics = dynamic_cast<HeatStructurePhysics *>(physics_action.get());
  if (!_physics)
    mooseError("Internal error: failed to create this component's HeatStructurePhysics");
  _awh.addActionBlock(physics_action);
}

void
HeatStructureInterfaceAC::addPhysics()
{
  _physics->addComponent(*this);
}
