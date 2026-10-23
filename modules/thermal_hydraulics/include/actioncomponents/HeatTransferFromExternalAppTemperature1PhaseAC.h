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
 * Applies a wall temperature provided by an external application to a FlowChannel1PhaseAC, the
 * native-ActionComponent counterpart of HeatTransferFromExternalAppTemperature1Phase.
 *
 * The wall temperature is a variable that the external application's transfer writes into. By
 * default it is named after the wall temperature and created by this component (with 'var_type'
 * selecting nodal or elemental, and 'initial_T_wall' its initial value). If 'T_ext' is given, the
 * variable of that name, which is created elsewhere (for example in an AuxVariables block), is
 * used instead.
 */
class HeatTransferFromExternalAppTemperature1PhaseAC : public HeatTransferFrom1PhaseAC
{
public:
  static InputParameters validParams();
  HeatTransferFromExternalAppTemperature1PhaseAC(const InputParameters & params);

protected:
  virtual bool isTemperatureType() const override { return true; }

  virtual void addHeatTransferVariables() override;
  virtual void addHeatTransferICs() override;
  virtual void addHeatTransferKernels() override;

  /// Name of the wall temperature variable
  virtual std::optional<VariableName> specifiedWallTemperatureName() const override;
};
