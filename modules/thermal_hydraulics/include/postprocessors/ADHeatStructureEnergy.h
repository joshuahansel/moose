//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralPostprocessor.h"

/**
 * Computes the total energy for a heat structure.
 */
class ADHeatStructureEnergy : public ElementIntegralPostprocessor
{
public:
  ADHeatStructureEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  /// Number of units that heat structure is multiplied by
  const Real _n_units;

  /// Reference temperature
  const Real & _T_ref;

  /// Density of the heat structure
  const ADMaterialProperty<Real> & _rho;

  /// Isobaric specific heat capacity
  const ADMaterialProperty<Real> & _cp;

  /// Temperature variable
  MooseVariable * _T_var;
  /// Temperature variable value
  const VariableValue & _T;

  /// Factor by which to scale integral
  const Real _scale;

public:
  static InputParameters validParams();
};
