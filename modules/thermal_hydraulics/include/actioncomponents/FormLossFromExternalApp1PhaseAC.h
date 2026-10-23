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
 * Applies a form loss over a FlowChannel1PhaseAC that is computed by an external application, the
 * native-ActionComponent counterpart of FormLossFromExternalApp1Phase. The form loss coefficient
 * per unit length is the nodal variable 'K_prime', which the external application's transfer
 * writes into.
 */
class FormLossFromExternalApp1PhaseAC : public FormLoss1PhaseBaseAC
{
public:
  static InputParameters validParams();
  FormLossFromExternalApp1PhaseAC(const InputParameters & params);

protected:
  virtual void addFormLossVariables() override;
  virtual void addFormLossMaterials() override;

  /// Name of the variable holding the form loss coefficient per unit length
  const VariableName _K_prime_var_name;
};
