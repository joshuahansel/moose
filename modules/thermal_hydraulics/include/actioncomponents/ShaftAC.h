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
 * Component that connects the torque of turbomachinery ActionComponents (ShaftConnectableAC)
 * through the rotational speed of a shaft, the native-ActionComponent counterpart of Shaft. The
 * speed is a scalar solver variable with a time derivative scaled by the total moment of inertia
 * and a source from each connected component's torque.
 */
class ShaftAC : public virtual ActionComponent, public NamingInterface
{
public:
  static InputParameters validParams();
  ShaftAC(const InputParameters & params);

  /// Name of the shaft speed variable
  const VariableName & omegaVariableName() const { return _omega_var_name; }

protected:
  virtual void addSolverVariables() override;
  virtual void actOnAdditionalTasks() override;

private:
  /// Checks the connected components and the initial speed, and tells the components about the
  /// shaft
  void connectComponents();
  /// Names of the user objects of the connected components
  std::vector<UserObjectName> connectedUserObjectNames() const;

  /// Scaling factor of the shaft speed
  const Real _scaling_factor_omega;
  /// Name of the shaft speed variable
  const VariableName _omega_var_name;
  /// Names of the connected components
  const std::vector<std::string> & _connected_components;
};
