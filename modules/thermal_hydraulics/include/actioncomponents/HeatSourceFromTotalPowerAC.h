//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HeatSourceBaseAC.h"

class TotalPowerAC;

/**
 * Heat generation in regions of a heat structure ActionComponent from the power of a
 * TotalPowerAC.
 *
 * Classic HeatSourceFromTotalPower picks RZ-specific objects for a cylindrical heat structure
 * (FunctionElementIntegralRZ, ADHeatStructureHeatSourceRZ), which multiply by the local
 * circumference. HeatStructureCylindricalAC uses MOOSE's native COORD_RZ coordinate system, which
 * applies that weighting itself, so the plain objects are used for both geometries.
 */
class HeatSourceFromTotalPowerAC : public HeatSourceBaseAC
{
public:
  static InputParameters validParams();
  HeatSourceFromTotalPowerAC(const InputParameters & params);

protected:
  virtual void actOnAdditionalTasks() override;

  /// The total power component, looked up by name lazily
  const TotalPowerAC & totalPower() const;

  /// Whether a power shape function was provided
  const bool _has_psf;
  /// Name of the power shape function (provided, or created here when it is not)
  const FunctionName _power_shape_fn;
  /// Name of the postprocessor integrating the power shape function over the source's subdomains
  const PostprocessorName _power_shape_integral_name;
};
