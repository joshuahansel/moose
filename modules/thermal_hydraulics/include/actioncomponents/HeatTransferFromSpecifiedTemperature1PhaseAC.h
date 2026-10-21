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
 * Applies a function-specified wall temperature to a FlowChannel1PhaseAC.
 */
class HeatTransferFromSpecifiedTemperature1PhaseAC : public HeatTransferFrom1PhaseAC
{
public:
  static InputParameters validParams();
  HeatTransferFromSpecifiedTemperature1PhaseAC(const InputParameters & params);

protected:
  virtual bool isTemperatureType() const override { return true; }

  virtual void addHeatTransferVariables() override;
  virtual void addHeatTransferICs() override;
  virtual void addHeatTransferKernels() override;
  virtual void addHeatTransferAuxKernels() override;
};
