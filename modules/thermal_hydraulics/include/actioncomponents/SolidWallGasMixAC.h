//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FlowBoundary1PhaseAC.h"

/**
 * Wall boundary condition for a FlowChannelGasMixAC, the counterpart of SolidWallGasMix. Its
 * boundary condition kernels are those of the gas mixture equations, which have the additional
 * equation of the secondary gas.
 */
class SolidWallGasMixAC : public FlowBoundary1PhaseAC
{
public:
  static InputParameters validParams();
  SolidWallGasMixAC(const InputParameters & params);

protected:
  virtual void addUserObjects() override;
  virtual void addWeakBCs() override;
};
