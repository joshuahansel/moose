//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SideIntegralPostprocessor.h"

/**
 * Integrates a conduction heat flux over a boundary.
 */
class HeatRateConduction : public SideIntegralPostprocessor
{
public:
  static InputParameters validParams();

  HeatRateConduction(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  /// Temperature gradient
  const VariableGradient & _grad_T;
  /// Thermal conductivity
  const ADMaterialProperty<Real> & _k;
  /// Inward/outward factor
  const Real _direction_factor;
};
