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
 * Heat structure ActionComponent boundary condition performing convective heat transfer with an
 * external application, which transfers the ambient temperature and heat transfer coefficient into
 * the 'T_ext'/'htc_ext' variables this component creates.
 */
class HSBoundaryExternalAppConvectionAC : public HSBoundaryAC
{
public:
  static InputParameters validParams();
  HSBoundaryExternalAppConvectionAC(const InputParameters & params);

protected:
  virtual void addAuxVariables() override;
  virtual void addBCs() override;
  virtual void addPostprocessors() override;
};
