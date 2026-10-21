//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ActionComponent.h"
#include "NamingInterface.h"

/**
 * Provides a total power, held in a scalar variable, for components such as
 * HeatSourceFromTotalPowerAC to use.
 */
class TotalPowerAC : public virtual ActionComponent, public NamingInterface
{
public:
  static InputParameters validParams();
  TotalPowerAC(const InputParameters & params);

  /// Name of the scalar variable holding the power
  const VariableName & getPowerVariableName() const { return _power_var_name; }

protected:
  virtual void actOnAdditionalTasks() override;

  /// Total power [W]
  const Real _power;
  /// Name of the scalar variable holding the power
  const VariableName _power_var_name;
};
