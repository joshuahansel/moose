//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FlowChannelGasMixAC.h"

// Registered under the classic Component's name ('FlowChannelGasMix'); see FlowChannel1PhaseAC.C
// for why. The tasks are the same as for FlowChannel1PhaseAC.
#define registerFlowChannelGasMixACTask(task)                                                      \
  registerMooseActionAliased("ThermalHydraulicsApp", FlowChannelGasMixAC, "FlowChannelGasMix", task)
registerFlowChannelGasMixACTask("add_mesh_generator");
registerFlowChannelGasMixACTask("init_component_physics");
registerFlowChannelGasMixACTask("add_functor_material");
registerFlowChannelGasMixACTask("add_postprocessor");
registerFlowChannelGasMixACTask("add_convergence");
registerFlowChannelGasMixACTask("add_vector_postprocessor");
registerActionComponentAliased("ThermalHydraulicsApp", FlowChannelGasMixAC, "FlowChannelGasMix");

InputParameters
FlowChannelGasMixAC::validParams()
{
  InputParameters params = FlowChannel1PhaseAC::validParams();
  params.set<std::string>("_physics_type") = "FlowChannel1PhaseGasMixPhysics";

  params.addRequiredParam<FunctionName>("initial_mass_fraction",
                                        "Initial mass fraction of the secondary gas");
  params.addParam<Real>(
      "scaling_factor_xirhoA", 1.0, "Scaling factor for the secondary component mass equation");
  params.addParam<Real>("scaling_factor_rhoA", 1.0, "Scaling factor for the mixture mass equation");
  params.addParam<Real>("scaling_factor_rhouA", 1.0, "Scaling factor for the momentum equation");
  params.addParam<Real>("scaling_factor_rhoEA", 1.0, "Scaling factor for the energy equation");
  params.addParamNamesToGroup("initial_mass_fraction", "Variable initialization");
  params.addParamNamesToGroup(
      "scaling_factor_xirhoA scaling_factor_rhoA scaling_factor_rhouA scaling_factor_rhoEA",
      "Numerical scheme");

  // These are given one by one instead
  params.suppressParameter<std::vector<Real>>("scaling_factor_1phase");
  // Passive transport is not supported for a gas mixture
  params.suppressParameter<std::vector<FunctionName>>("initial_passives");
  params.suppressParameter<std::vector<VariableName>>("passives_names");
  params.suppressParameter<std::vector<Real>>("scaling_factor_passives");

  params.addClassDescription("Single-phase flow channel with a binary gas mixture");
  return params;
}

FlowChannelGasMixAC::FlowChannelGasMixAC(const InputParameters & params)
  : ActionComponent(params), FlowChannel1PhaseAC(params)
{
}

std::optional<ConvergenceName>
FlowChannelGasMixAC::nonlinearConvergenceName() const
{
  mooseError("The flow channel '",
             name(),
             "' is a gas mixture flow channel, which does not provide a nonlinear convergence.");
}
