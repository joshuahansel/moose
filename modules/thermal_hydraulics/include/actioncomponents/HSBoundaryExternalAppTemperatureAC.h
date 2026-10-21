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
 * Heat structure ActionComponent boundary condition setting the temperature to values computed by
 * an external application (transferred into the 'T_ext' variable this component creates).
 */
class HSBoundaryExternalAppTemperatureAC : public HSBoundaryAC
{
public:
  static InputParameters validParams();
  HSBoundaryExternalAppTemperatureAC(const InputParameters & params);

protected:
  virtual void addAuxVariables() override;
  virtual void addBCs() override;
};
