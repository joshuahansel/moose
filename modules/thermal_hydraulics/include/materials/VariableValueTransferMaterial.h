//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "VariableTransferMaterialBase.h"

/**
 * Creates an AD material property for a variable transferred from the boundary
 * of a 2D mesh onto a 1D mesh.
 */
class VariableValueTransferMaterial : public VariableTransferMaterialBase
{
public:
  VariableValueTransferMaterial(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  virtual void computeProperties() override;

  /// Basis function for transferred variable
  const VariablePhiValue & _phi;
  /// Material property for the transferred variable value
  ADMaterialProperty<Real> & _prop;
};
