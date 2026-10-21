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

/**
 * Heat generation in regions of a heat structure ActionComponent from a power density variable.
 *
 * Classic HeatSourceFromPowerDensity picks CoupledForceRZ for a cylindrical heat structure.
 * HeatStructureCylindricalAC uses MOOSE's native COORD_RZ coordinate system, which applies the RZ
 * weighting itself, so CoupledForce is used for both geometries.
 */
class HeatSourceFromPowerDensityAC : public HeatSourceBaseAC
{
public:
  static InputParameters validParams();
  HeatSourceFromPowerDensityAC(const InputParameters & params);

protected:
  virtual void actOnAdditionalTasks() override;
};
