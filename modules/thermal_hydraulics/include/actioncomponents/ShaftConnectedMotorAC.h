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
#include "ShaftConnectableAC.h"

/**
 * Motor that drives a ShaftAC with a torque, the native-ActionComponent counterpart of
 * ShaftConnectedMotor. The torque and the moment of inertia are functions of the shaft speed.
 */
class ShaftConnectedMotorAC : public virtual ActionComponent,
                              public NamingInterface,
                              public ShaftConnectableAC
{
public:
  static InputParameters validParams();
  ShaftConnectedMotorAC(const InputParameters & params);

  virtual UserObjectName shaftConnectedUserObjectName() const override { return _uo_name; }

protected:
  virtual void addMeshGenerators() override;
  virtual void addUserObjects() override;

  /// Name of the user object computing the torque and the moment of inertia
  const UserObjectName _uo_name;
};
