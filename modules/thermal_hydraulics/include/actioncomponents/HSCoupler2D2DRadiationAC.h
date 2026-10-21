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

class HeatStructureBaseAC;
class MeshAlignment2D2D;

/**
 * Couples boundaries of several 2D cylindrical heat structure ActionComponents via radiation,
 * optionally with an environment, the native-ActionComponent counterpart of HSCoupler2D2DRadiation.
 *
 * The heat flux on each surface comes from the radiosity of an enclosure of the surfaces, solved
 * per quadrature point from the temperatures of the surfaces next to each other
 * (HSCoupler2D2DRadiationUserObject). The surfaces are paired with a MeshAlignment2D2D, and the
 * Jacobian needs an explicit coupling of all the paired elements, added with
 * AugmentSparsityBetweenElements (as for JunctionOneToOne1PhaseAC).
 *
 * Classic THM applies the heat flux of a cylindrical heat structure with an RZ-specific boundary
 * condition that multiplies by the local circumference. HeatStructureCylindricalAC uses MOOSE's
 * native COORD_RZ coordinate system, which applies that factor itself, so the plain boundary
 * condition is used.
 */
class HSCoupler2D2DRadiationAC : public virtual ActionComponent, public NamingInterface
{
public:
  static InputParameters validParams();
  HSCoupler2D2DRadiationAC(const InputParameters & params);

protected:
  virtual void actOnAdditionalTasks() override;
  virtual void addUserObjects() override;
  virtual void addRelationshipManagers(Moose::RelationshipManagerType when_type) override;

private:
  /// The i-th heat structure, which must be a cylindrical heat structure. The lookup is done
  /// lazily so that the heat structures may be declared after this component.
  const HeatStructureBaseAC & hs(unsigned int i) const;
  /// Checks the inputs
  void checkInputs();
  /// Builds the alignment between the boundaries if it has not been built yet
  MeshAlignment2D2D & meshAlignment();
  /// Adds the boundary conditions
  void addBCs();

  /// Names of the heat structures
  const std::vector<std::string> _hs_names;
  /// Names of the boundaries of the heat structures that are coupled
  const std::vector<BoundaryName> _hs_boundaries;
  /// Whether an environment surrounds all of the surfaces
  const bool _include_environment;
  /// Number of heat structures
  const std::size_t _n_hs;
  /// Number of surfaces: the heat structures, and the environment if it is included
  const std::size_t _n_surfaces;
  /// Name of the user object caching the temperature of each boundary by element
  const UserObjectName _temperature_uo_name;
  /// Name of the user object computing the heat fluxes
  const UserObjectName _coupler_uo_name;

  /// Alignment of the boundaries
  std::unique_ptr<MeshAlignment2D2D> _mesh_alignment;
  /// Map of element ID -> coupled element IDs, consumed by AugmentSparsityBetweenElements
  std::map<dof_id_type, std::vector<dof_id_type>> _elem_map;
  /// Whether the AugmentSparsityBetweenElements relationship manager has been added
  bool _rm_added;
};
