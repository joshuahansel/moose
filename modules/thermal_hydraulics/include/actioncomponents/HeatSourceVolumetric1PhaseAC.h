//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ActionComponent.h"
#include "NamingInterface.h"

/// Volumetric heat source applied on a single-phase flow channel ActionComponent
class HeatSourceVolumetric1PhaseAC : public virtual ActionComponent, public NamingInterface
{
public:
  static InputParameters validParams();
  HeatSourceVolumetric1PhaseAC(const InputParameters & params);

protected:
  virtual void actOnAdditionalTasks() override;
};
