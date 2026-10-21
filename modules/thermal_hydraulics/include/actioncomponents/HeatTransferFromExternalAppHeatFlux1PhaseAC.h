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
 * Applies a wall heat flux provided by an external application to a FlowChannel1PhaseAC, the
 * native-ActionComponent counterpart of HeatTransferFromExternalAppHeatFlux1Phase.
 *
 * The heat flux is a constant monomial auxiliary variable named after the wall heat flux, which
 * the external application's transfer writes into (and which starts at zero).
 */
class HeatTransferFromExternalAppHeatFlux1PhaseAC : public HeatTransferFrom1PhaseAC
{
public:
  static InputParameters validParams();
  HeatTransferFromExternalAppHeatFlux1PhaseAC(const InputParameters & params);

protected:
  virtual bool isTemperatureType() const override { return false; }

  virtual void addHeatTransferVariables() override;
  virtual void addHeatTransferICs() override;
  virtual void addHeatTransferMaterials() override;
  virtual void addHeatTransferKernels() override;
};
