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
 * Couples boundaries of two 2D heat structure ActionComponents via a heat transfer coefficient, the
 * native-ActionComponent counterpart of HeatStructure2DCoupler.
 */
class HeatStructure2DCouplerAC : public HeatStructure2DCouplerBaseAC
{
public:
  static InputParameters validParams();
  HeatStructure2DCouplerAC(const InputParameters & params);

protected:
  virtual void checkCoupling() override;
  virtual void addBCs() override;
};
