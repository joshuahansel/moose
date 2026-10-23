//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "JunctionOneToOne1PhaseAC.h"

/**
 * Gate valve between two FlowChannel1PhaseAC boundaries for 1-phase flow, the
 * native-ActionComponent counterpart of GateValve1Phase. The fraction of the flow area that is open
 * can be controlled.
 */
class GateValve1PhaseAC : public JunctionOneToOne1PhaseAC
{
public:
  static InputParameters validParams();
  GateValve1PhaseAC(const InputParameters & params);

protected:
  virtual void buildJunctionUserObject() override;
  virtual void addWeakBC(unsigned int i, const NonlinearVariableName & var) override;
};
