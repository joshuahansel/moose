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
 * Creates an AD material property for a variable gradient transferred from the boundary
 * of a 2D mesh onto a 1D mesh.
 */
class VariableGradientTransferMaterial : public VariableTransferMaterialBase
{
public:
  VariableGradientTransferMaterial(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  virtual void computeProperties() override;

  /// Gradient of basis function for transferred variable
  const VariablePhiGradient & _grad_phi;
  /// Material property for the transferred variable gradient
  ADMaterialProperty<RealVectorValue> & _prop;
};
