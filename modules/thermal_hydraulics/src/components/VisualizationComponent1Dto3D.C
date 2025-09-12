//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VisualizationComponent1Dto3D.h"

registerMooseObject("ThermalHydraulicsApp", VisualizationComponent1Dto3D);

InputParameters
VisualizationComponent1Dto3D::validParams()
{
  InputParameters params = Component3D::validParams();
  return params;
}

VisualizationComponent1Dto3D::VisualizationComponent1Dto3D(const InputParameters & params)
  : Component3D(params)
{
}

bool
VisualizationComponent1Dto3D::usingSecondOrderMesh() const
{
  return false;
}

