//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GrayLambertRadiationHeatFluxAux.h"
#include "GrayLambertSurfaceRadiationBase.h"

registerMooseObject("HeatTransferApp", GrayLambertRadiationHeatFluxAux);

InputParameters
GrayLambertRadiationHeatFluxAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Radiation heat flux from a GrayLambertSurfaceRadiationBase object.");
  params.addRequiredParam<UserObjectName>("surface_radiation_object",
                                          "GrayLambertSurfaceRadiationBase UO name");
  params.addRequiredParam<std::vector<BoundaryName>>("radiation_patch_boundary",
                                                     "Patch radiation boundaries");
  return params;
}

GrayLambertRadiationHeatFluxAux::GrayLambertRadiationHeatFluxAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _glsr_uo(getUserObject<GrayLambertSurfaceRadiationBase>("surface_radiation_object"))
{
  if (isNodal())
    mooseError("The variable may not be nodal because the heat fluxes are looked up by the current "
               "boundary ID, which is ambiguous for a given node.");

  const auto radiation_patch_boundary =
      getParam<std::vector<BoundaryName>>("radiation_patch_boundary");
  for (const auto & boundary : radiation_patch_boundary)
    _radiation_boundary_ids.push_back(_mesh.getBoundaryID(boundary));
  // std::set_intersection requires vectors to be sorted
  std::sort(_radiation_boundary_ids.begin(), _radiation_boundary_ids.end());
}

Real
GrayLambertRadiationHeatFluxAux::computeValue()
{
  return -_glsr_uo.getSurfaceHeatFluxDensity(getRadiationPatchBoundaryID());
}

BoundaryID
GrayLambertRadiationHeatFluxAux::getRadiationPatchBoundaryID() const
{
  // get the interior parent element and side
  const auto * higher_d_elem = _current_elem->interior_parent();
  if (!higher_d_elem)
    mooseError("Element ", _current_elem->id(), " has no interior parent");
  const auto higher_d_side = _mesh.getHigherDSide(_current_elem);

  // get the boundary IDs of the interior parent side
  std::vector<BoundaryID> boundary_ids;
  _mesh.getMesh().get_boundary_info().boundary_ids(higher_d_elem, higher_d_side, boundary_ids);

  // get only the boundary IDs corresponding to the radiation patch boundaries; there
  // should be exactly 1
  std::vector<BoundaryID> intersection_boundary_ids;
  std::sort(boundary_ids.begin(), boundary_ids.end());
  std::set_intersection(boundary_ids.begin(),
                        boundary_ids.end(),
                        _radiation_boundary_ids.begin(),
                        _radiation_boundary_ids.end(),
                        std::back_inserter(intersection_boundary_ids));
  if (intersection_boundary_ids.size() != 1)
    mooseError("Element ",
               _current_elem->id(),
               ", side ",
               higher_d_side,
               " has ",
               intersection_boundary_ids.size(),
               " radiation patch boundary IDs but must have exactly 1.");

  return intersection_boundary_ids[0];
}
