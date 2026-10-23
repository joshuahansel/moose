//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FormLoss1PhaseBaseAC.h"

/**
 * Prescribes a form loss over a FlowChannel1PhaseAC, given by a function, the
 * native-ActionComponent counterpart of FormLossFromFunction1Phase.
 */
class FormLossFromFunction1PhaseAC : public FormLoss1PhaseBaseAC
{
public:
  static InputParameters validParams();
  FormLossFromFunction1PhaseAC(const InputParameters & params);

protected:
  virtual void addFormLossMaterials() override;
};
