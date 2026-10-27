//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatStructureCylindricalAC.h"
#include "MooseMesh.h"
#include "MultiMooseEnum.h"

#include <numeric>

// Registered under the classic Component's name ('HeatStructureCylindrical') rather than this
// class's own name; see FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructureCylindricalAC,
                           "HeatStructureCylindrical",
                           "add_mesh_generator");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructureCylindricalAC,
                           "HeatStructureCylindrical",
                           "init_component_physics");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructureCylindricalAC,
                           "HeatStructureCylindrical",
                           "add_material");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructureCylindricalAC,
                           "HeatStructureCylindrical",
                           "setup_component");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructureCylindricalAC,
                           "HeatStructureCylindrical",
                           "check_integrity");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructureCylindricalAC,
                           "HeatStructureCylindrical",
                           "add_functor_material");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructureCylindricalAC,
                           "HeatStructureCylindrical",
                           "add_postprocessor");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructureCylindricalAC,
                           "HeatStructureCylindrical",
                           "add_convergence");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HeatStructureCylindricalAC,
                               "HeatStructureCylindrical");

InputParameters
HeatStructureCylindricalAC::validParams()
{
  InputParameters params = HeatStructureBaseAC::validParams();
  params.addParam<Real>("inner_radius", 0., "Inner radius of the heat structure [m]");
  params.addClassDescription("Cylindrical heat structure.");
  return params;
}

HeatStructureCylindricalAC::HeatStructureCylindricalAC(const InputParameters & params)
  : ActionComponent(params),
    HeatStructureBaseAC(params),
    _inner_radius(getParam<Real>("inner_radius"))
{
  addRequiredTask("setup_component");
}

Real
HeatStructureCylindricalAC::getUnitPerimeter(ExternalBoundaryType side) const
{
  switch (side)
  {
    case ExternalBoundaryType::INNER:
      return 2 * libMesh::pi * _inner_radius;
    case ExternalBoundaryType::OUTER:
      return 2 * libMesh::pi *
             (_inner_radius + std::accumulate(_widths.begin(), _widths.end(), 0.0));
    case ExternalBoundaryType::START:
    case ExternalBoundaryType::END:
      return std::numeric_limits<Real>::quiet_NaN();
  }
  mooseError(name(), ": Unknown external boundary type.");
}

Real
HeatStructureCylindricalAC::computeRadialBoundaryArea(Real length, Real y) const
{
  return length * 2 * libMesh::pi * (_inner_radius + y);
}

Real
HeatStructureCylindricalAC::computeAxialBoundaryArea(Real y_min, Real y_max) const
{
  return libMesh::pi * (std::pow(_inner_radius + y_max, 2) - std::pow(_inner_radius + y_min, 2));
}

void
HeatStructureCylindricalAC::setupComponent()
{
  // MooseMesh::setCoordSystem()/setGeneralAxisymmetricCoordAxes() must each be called exactly once
  // per simulation, covering every RZ block at once (see class documentation) - so only the first
  // HeatStructureCylindricalAC instance (in action-warehouse order) does anything here, making one
  // combined call on behalf of every instance; the rest are no-ops.
  const auto instances = _awh.getActions<HeatStructureCylindricalAC>();
  if (instances.empty() || instances[0] != this)
    return;

  std::vector<SubdomainName> all_blocks;
  std::vector<std::pair<Point, RealVectorValue>> all_axes;
  for (const auto * instance : instances)
  {
    for (const auto & block : instance->_blocks)
      all_blocks.push_back(block);
    // The heat structure is generated along the x axis, which is its direction if no orientation
    // is given
    const RealVectorValue direction =
        instance->isParamValid("orientation") ? instance->direction() : RealVectorValue(1, 0, 0);
    all_axes.insert(all_axes.end(),
                    instance->_blocks.size(),
                    std::make_pair(instance->translation(), direction));
  }

  // Every subdomain of the mesh is listed, with RZ for the cylindrical heat structure blocks and
  // XYZ for the others (e.g. flow channels): MooseMesh::setCoordSystem() applies a single
  // coordinate system to all of the subdomains if only one block is given, and requires every
  // subdomain to be listed otherwise.
  //
  // Note: MultiMooseEnum's single-string constructor sets the *valid options* list, not the
  // selected value - the two-argument form is required to actually select a value, and the valid
  // Moose::CoordinateSystemType strings are "XYZ"/"RZ"/"RSPHERICAL" (see Conversion.C), not the
  // C++ enumerator names themselves.
  auto & mesh = *_awh.getMesh();
  std::vector<SubdomainName> mesh_blocks;
  std::string coord_systems;
  for (const auto subdomain_id : mesh.meshSubdomains())
  {
    SubdomainName block = mesh.getSubdomainName(subdomain_id);
    if (block.empty())
      block = std::to_string(subdomain_id);
    const bool is_rz = std::find(all_blocks.begin(), all_blocks.end(), block) != all_blocks.end();
    mesh_blocks.push_back(block);
    coord_systems += (coord_systems.empty() ? "" : " ") + std::string(is_rz ? "RZ" : "XYZ");
  }
  mesh.setCoordSystem(mesh_blocks, MultiMooseEnum("XYZ RZ", coord_systems));
  mesh.setGeneralAxisymmetricCoordAxes(all_blocks, all_axes);
}
