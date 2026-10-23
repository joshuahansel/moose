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

/**
 * Junction between FlowChannel1PhaseAC boundaries whose flow channels are parallel, the
 * native-ActionComponent counterpart of JunctionParallelChannels1Phase. It differs from
 * VolumeJunction1PhaseAC only in the user object that computes the junction fluxes.
 */
class JunctionParallelChannels1PhaseAC : public VolumeJunction1PhaseAC
{
public:
  static InputParameters validParams();
  JunctionParallelChannels1PhaseAC(const InputParameters & params);

protected:
  virtual void buildVolumeJunctionUserObject() override;

  /// Direction of the flow channel connected first, which the others must be parallel to
  RealVectorValue firstDirection() const;

  /// Name of the class of the junction user object. Derived classes with other physics override it.
  virtual std::string junctionUserObjectClass() const
  {
    return "ADJunctionParallelChannels1PhaseUserObject";
  }
  /// Sets the parameters of the junction user object that are specific to a derived class
  virtual void setJunctionUserObjectParameters(InputParameters & /*params*/) const {}
  /// Connects the controllable parameters of a derived class to the junction user object
  virtual void connectJunctionUserObjectControllableParams(const std::string & /*class_name*/) {}
};
