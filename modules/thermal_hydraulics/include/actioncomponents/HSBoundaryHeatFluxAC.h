//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HSBoundaryAC.h"

/**
 * Applies a specified heat flux to a heat structure ActionComponent boundary.
 */
class HSBoundaryHeatFluxAC : public HSBoundaryAC
{
public:
  static InputParameters validParams();
  HSBoundaryHeatFluxAC(const InputParameters & params);

protected:
  virtual void addBCs() override;
  virtual void addPostprocessors() override;
};
