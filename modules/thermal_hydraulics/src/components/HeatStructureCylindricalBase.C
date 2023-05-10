//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatStructureCylindricalBase.h"
#include "MooseUtils.h"

InputParameters
HeatStructureCylindricalBase::validParams()
{
  InputParameters params = HeatStructureBase::validParams();
  params.addParam<bool>(
      "offset_mesh_by_inner_radius", false, "Offset the mesh by the inner radius?");
  return params;
}

HeatStructureCylindricalBase::HeatStructureCylindricalBase(const InputParameters & params)
  : HeatStructureBase(params)
{
}

void
HeatStructureCylindricalBase::setupMesh()
{
  const Real inner_radius = getInnerRadius();

  if (getParam<bool>("offset_mesh_by_inner_radius") || !_connected_to_flow_channel)
    _axial_offset = inner_radius;
  else if (!MooseUtils::absoluteFuzzyEqual(inner_radius, 0.0))
    mooseDeprecated(
        "Cylindrical heat structure meshes must now be offset by their inner radii. Set "
        "'offset_mesh_by_inner_radius = true', and re-gold any output files depending "
        "on heat structure mesh position.");

  HeatStructureBase::setupMesh();
}
