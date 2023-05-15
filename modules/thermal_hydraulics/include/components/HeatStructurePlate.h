//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HeatStructureBase.h"
#include "HeatConductionModel.h"

/**
 * Component to model plate heat structure
 */
class HeatStructurePlate : public HeatStructureBase
{
public:
  HeatStructurePlate(const InputParameters & params);

  virtual void check() const override;

  virtual bool isCylindrical() const override { return false; }
  virtual std::vector<std::string> getRegionNames() const override;
  virtual std::vector<Real> getRegionWidths() const override;
  virtual std::vector<unsigned int> getRegionNumbersOfElements() const override;
  virtual Real getNumberOfUnits() const override;
  virtual Real getDepth() const override { return _depth; }

protected:
  virtual bool useCylindricalTransformation() const override { return false; }

  /// plate fuel depth
  const Real & _depth;

public:
  static InputParameters validParams();
};
