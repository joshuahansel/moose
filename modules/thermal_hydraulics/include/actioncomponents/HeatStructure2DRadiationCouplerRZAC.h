//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HeatStructure2DCouplerBaseAC.h"

/**
 * Couples boundaries of two 2D cylindrical heat structure ActionComponents via radiation, the
 * native-ActionComponent counterpart of HeatStructure2DRadiationCouplerRZ.
 */
class HeatStructure2DRadiationCouplerRZAC : public HeatStructure2DCouplerBaseAC
{
public:
  static InputParameters validParams();
  HeatStructure2DRadiationCouplerRZAC(const InputParameters & params);

protected:
  virtual void checkCoupling() override;
  virtual void addBCs() override;

  /// Emissivities of the primary and secondary sides
  const std::array<Real, 2> _emissivities;
};
