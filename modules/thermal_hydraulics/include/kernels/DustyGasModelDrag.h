//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernel.h"

/**
 * Adds drag according to the Dusty Gas Model.
 */
class DustyGasModelDrag : public ADKernel
{
public:
  static InputParameters validParams();

  DustyGasModelDrag(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// Area
  const ADVariableValue & _A;
  /// Density
  const ADMaterialProperty<Real> & _rho;
  /// Velocity
  const ADMaterialProperty<Real> & _vel;
  /// Pressure
  const ADMaterialProperty<Real> & _p;
  /// Dynamic viscosity
  const ADMaterialProperty<Real> & _mu;
};
