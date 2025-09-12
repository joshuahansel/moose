//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Component3D.h"

/**
 * Base class for 2D generated heat structures
 */
class VisualizationComponent1Dto3D : public Component3D
{
public:
  static InputParameters validParams();

  VisualizationComponent1Dto3D(const InputParameters & params);

  bool usingSecondOrderMesh() const override;
};
