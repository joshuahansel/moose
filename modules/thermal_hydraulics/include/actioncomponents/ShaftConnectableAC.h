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
 * Interface for ActionComponents that connect to a shaft, the native-ActionComponent counterpart of
 * ShaftConnectable. The component provides the name of the user object that gives the torque and
 * the moment of inertia it contributes to the shaft equation. The shaft sets the name of the shaft
 * when it is connected, which the component needs to couple to the speed of the shaft.
 */
class ShaftConnectableAC
{
public:
  /// Name of the user object computing the torque and the moment of inertia of the component
  virtual UserObjectName shaftConnectedUserObjectName() const = 0;

  /// Sets the name of the shaft this component is connected to
  void setShaftName(const std::string & shaft_name) const
  {
    _shaft_name = shaft_name;
    _connected_to_shaft = true;
  }

  /// Whether a shaft has this component connected
  bool connectedToShaft() const { return _connected_to_shaft; }

protected:
  /// Name of the shaft component
  mutable std::string _shaft_name;
  /// Flag indicating that a shaft has this component connected
  mutable bool _connected_to_shaft = false;
};
