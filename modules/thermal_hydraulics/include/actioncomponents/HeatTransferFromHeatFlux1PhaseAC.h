//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HeatTransferFrom1PhaseAC.h"

/**
 * Applies a function-specified wall heat flux to a FlowChannel1PhaseAC.
 */
class HeatTransferFromHeatFlux1PhaseAC : public HeatTransferFrom1PhaseAC
{
public:
  static InputParameters validParams();
  HeatTransferFromHeatFlux1PhaseAC(const InputParameters & params);

protected:
  virtual bool isTemperatureType() const override { return false; }

  virtual void addHeatTransferMaterials() override;
  virtual void addHeatTransferKernels() override;
};
