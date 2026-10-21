//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADIntegratedBC.h"

class HSCoupler2D2DRadiationUserObject;

/**
 * Adds heat flux terms for HSCoupler2D2DRadiation in a coordinate system that already accounts for
 * the geometry (Cartesian, or the native RZ coordinate system)
 */
class HSCoupler2D2DRadiationBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  HSCoupler2D2DRadiationBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// HSCoupler2D2DRadiation user object
  const HSCoupler2D2DRadiationUserObject & _hs_coupler_2d2d_uo;
};
