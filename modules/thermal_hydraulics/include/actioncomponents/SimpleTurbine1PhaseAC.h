//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "JunctionParallelChannels1PhaseAC.h"

/**
 * Simple turbine model that extracts a prescribed power from the working fluid, between parallel
 * FlowChannel1PhaseAC boundaries, the native-ActionComponent counterpart of SimpleTurbine1Phase.
 * The power is also stored in the junction auxiliary variable 'W_dot'.
 */
class SimpleTurbine1PhaseAC : public JunctionParallelChannels1PhaseAC
{
public:
  static InputParameters validParams();
  SimpleTurbine1PhaseAC(const InputParameters & params);

protected:
  virtual void actOnAdditionalTasks() override;
  virtual std::string junctionUserObjectClass() const override
  {
    return "ADSimpleTurbine1PhaseUserObject";
  }
  virtual void setJunctionUserObjectParameters(InputParameters & params) const override;
  virtual void connectJunctionUserObjectControllableParams(const std::string & class_name) override;

  /// Adds the auxiliary kernel storing the turbine power in the 'W_dot' variable
  void addPowerAuxKernel();

  /// Whether the turbine is operating
  const bool _on;
  /// Turbine power
  const Real _power;
  /// Name of the variable holding the turbine power
  const VariableName _W_dot_var_name;
};
