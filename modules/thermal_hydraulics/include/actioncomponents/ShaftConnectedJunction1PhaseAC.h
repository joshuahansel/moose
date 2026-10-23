//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "VolumeJunction1PhaseAC.h"
#include "ShaftConnectableAC.h"

/**
 * Base class for turbomachinery junctions that connect to a ShaftAC: ShaftConnectedPump1PhaseAC,
 * ShaftConnectedTurbine1PhaseAC and ShaftConnectedCompressor1PhaseAC. They take an 'inlet' and an
 * 'outlet' instead of 'connections', have the friction torque and moment of inertia of the rotating
 * parts as variables, and a single user object that serves as both the junction user object and
 * the shaft-connected user object.
 *
 * Derived classes provide the auxiliary variables and kernels of their own quantities.
 */
class ShaftConnectedJunction1PhaseAC : public VolumeJunction1PhaseAC, public ShaftConnectableAC
{
public:
  static InputParameters validParams();
  ShaftConnectedJunction1PhaseAC(const InputParameters & params);

  virtual UserObjectName shaftConnectedUserObjectName() const override { return _junction_uo_name; }

protected:
  virtual void addMeshGenerators() override;
  virtual void addAuxiliaryVariables() override;
  virtual void addInitialConditions() override;
  virtual void addAuxiliaryKernels() override;

  /// Requests a constant monomial junction auxiliary variable
  void addConstantJunctionVariable(const VariableName & var);
  /// Adds the auxiliary kernel for each pair of (quantity, variable name) with the class and the
  /// parameter naming the shaft-connected user object
  void addQuantityAuxKernels(const std::string & class_name,
                             const std::string & uo_param,
                             const std::vector<std::pair<std::string, VariableName>> & quantities);

  /// The auxiliary variables of the derived class and their initial values
  virtual std::vector<std::pair<VariableName, Real>> shaftConnectedVariables() const = 0;
  /// The class of the auxiliary kernels of the derived class
  virtual std::string auxKernelClass() const = 0;
  /// The name of the parameter of the auxiliary kernels that takes the user object
  virtual std::string auxKernelUserObjectParam() const = 0;
  /// The quantities that the auxiliary kernels compute, and the variable of each
  virtual std::vector<std::pair<std::string, VariableName>> auxKernelQuantities() const = 0;

  /// Parameters shared by the shaft-connected user objects: the junction parameters, the outlet
  /// direction, the shaft speed and the rated, friction and inertia parameters
  InputParameters shaftConnectedUserObjectParameters(const std::string & class_name) const;

  /// Direction of the outlet flow channel, which is the second connection
  RealVectorValue outletDirection() const { return connectedDirection(1).unit(); }

  /// Name of the shaft speed variable of the connected shaft
  VariableName omegaVariableName() const;

  /// Pump/turbine/compressor inertia constant
  const Real _inertia_const;
  /// Name of the friction torque variable
  const VariableName _friction_torque_var_name;
  /// Name of the moment of inertia variable
  const VariableName _moment_of_inertia_var_name;
};
