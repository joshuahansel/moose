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
 * Pump between two FlowChannel1PhaseAC boundaries that has a non-zero volume, the
 * native-ActionComponent counterpart of Pump1Phase. It differs from VolumeJunction1PhaseAC in the
 * user object that computes the junction fluxes, which adds the pump head as a momentum source in
 * the direction of the outlet.
 */
class Pump1PhaseAC : public VolumeJunction1PhaseAC
{
public:
  static InputParameters validParams();
  Pump1PhaseAC(const InputParameters & params);

protected:
  virtual void buildVolumeJunctionUserObject() override;

  /// Pump head
  const Real _head;
};
