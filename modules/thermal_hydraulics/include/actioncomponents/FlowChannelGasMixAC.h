//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FlowChannel1PhaseAC.h"

/**
 * Single-phase flow channel with a binary gas mixture as a native ActionComponent, the counterpart
 * of FlowChannelGasMix. It differs from FlowChannel1PhaseAC in its equations, which are those of
 * FlowChannel1PhaseGasMixPhysics, and in having no nonlinear convergence object, as in classic THM.
 */
class FlowChannelGasMixAC : public FlowChannel1PhaseAC
{
public:
  static InputParameters validParams();
  FlowChannelGasMixAC(const InputParameters & params);

  virtual std::optional<ConvergenceName> nonlinearConvergenceName() const override;

protected:
  virtual void addConvergence() override {}
};
