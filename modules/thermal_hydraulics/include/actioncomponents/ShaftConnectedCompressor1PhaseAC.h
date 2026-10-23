//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ShaftConnectedJunction1PhaseAC.h"

/**
 * Compressor between two 1-phase flow channels whose speed is given by a connected shaft, the
 * native-ActionComponent counterpart of ShaftConnectedCompressor1Phase.
 */
class ShaftConnectedCompressor1PhaseAC : public ShaftConnectedJunction1PhaseAC
{
public:
  static InputParameters validParams();
  ShaftConnectedCompressor1PhaseAC(const InputParameters & params);

protected:
  virtual void buildVolumeJunctionUserObject() override;
  virtual std::vector<std::pair<VariableName, Real>> shaftConnectedVariables() const override;
  virtual std::string auxKernelClass() const override
  {
    return "ShaftConnectedCompressor1PhaseAux";
  }
  virtual std::string auxKernelUserObjectParam() const override { return "compressor_uo"; }
  virtual std::vector<std::pair<std::string, VariableName>> auxKernelQuantities() const override;

  virtual void actOnAdditionalTasks() override;
};
