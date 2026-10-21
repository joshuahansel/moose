//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HeatStructure2DCouplerBCBase.h"

/**
 * Applies BC for a radiation coupler between two 2D heat structures in a coordinate system that
 * already accounts for the geometry (Cartesian, or the native RZ coordinate system)
 */
class HeatStructure2DRadiationCouplerBC : public HeatStructure2DCouplerBCBase
{
public:
  HeatStructure2DRadiationCouplerBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// Emissivity of this boundary
  const Real & _emissivity;
  /// Emissivity of the coupled boundary
  const Real & _coupled_emissivity;
  /// View factor of this boundary
  const Real & _view_factor;
  /// Area of this boundary
  const Real & _area;
  /// Area of the coupled boundary
  const Real & _coupled_area;
  /// Stefan-Boltzmann constant
  const Real & _sigma;
  /// Radiation resistance
  const Real _radiation_resistance;

public:
  static InputParameters validParams();
};
