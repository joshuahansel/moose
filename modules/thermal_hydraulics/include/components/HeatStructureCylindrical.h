//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HeatStructureCylindricalBase.h"

/**
 * Component to model cylindrical heat structure
 */
class HeatStructureCylindrical : public HeatStructureCylindricalBase
{
public:
  HeatStructureCylindrical(const InputParameters & params);

  virtual void check() const override;

  virtual std::vector<std::string> getRegionNames() const override;
  virtual std::vector<Real> getRegionWidths() const override;
  virtual std::vector<unsigned int> getRegionNumbersOfElements() const override;
  virtual Real getNumberOfUnits() const override;
  virtual Real getInnerRadius() const override;

public:
  static InputParameters validParams();
};
