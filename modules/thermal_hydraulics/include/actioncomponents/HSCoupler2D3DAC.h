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

class HeatStructureCylindricalAC;
class HeatStructureFromFile3DAC;
class MeshAlignment2D3D;

/**
 * Couples a boundary of a 2D cylindrical heat structure ActionComponent to a boundary of a 3D heat
 * structure ActionComponent loaded from a file, using gap heat transfer, the native-ActionComponent
 * counterpart of HSCoupler2D3D.
 *
 * The heat flux across the gap (conduction, convection and optionally radiation) is computed on
 * the 3D boundary from the temperature of the 3D side and the temperature of the 2D side next to
 * it (HSCoupler2D3DUserObject), and applied to both boundaries. The boundaries are paired with a
 * MeshAlignment2D3D, and the Jacobian needs an explicit coupling of the paired elements, added with
 * AugmentSparsityBetweenElements (as for JunctionOneToOne1PhaseAC).
 */
class HSCoupler2D3DAC : public virtual ActionComponent, public NamingInterface
{
public:
  static InputParameters validParams();
  HSCoupler2D3DAC(const InputParameters & params);

protected:
  virtual void actOnAdditionalTasks() override;
  virtual void addUserObjects() override;
  virtual void addRelationshipManagers(Moose::RelationshipManagerType when_type) override;

private:
  /// The 2D heat structure, which must be a cylindrical heat structure (looked up lazily so that
  /// the heat structures may be declared after this component)
  const HeatStructureCylindricalAC & hs2D() const;
  /// The 3D heat structure, which must be loaded from a file
  const HeatStructureFromFile3DAC & hs3D() const;
  /// Checks the inputs
  void checkInputs();
  /// Builds the alignment between the boundaries if it has not been built yet
  MeshAlignment2D3D & meshAlignment();
  /// Adds the boundary conditions
  void addBCs();

  /// Name of the 2D heat structure
  const std::string _hs_name_2d;
  /// Name of the 3D heat structure
  const std::string _hs_name_3d;
  /// Boundary of the 2D heat structure
  const BoundaryName _boundary_2d;
  /// Boundary of the 3D heat structure
  const BoundaryName _boundary_3d;
  /// Name of the user object caching the temperature of the 2D boundary by element
  const UserObjectName _temperature_2d_uo_name;
  /// Name of the user object computing the heat fluxes
  const UserObjectName _coupler_uo_name;

  /// Alignment of the boundaries
  std::unique_ptr<MeshAlignment2D3D> _mesh_alignment;
  /// Map of element ID -> coupled element IDs, consumed by AugmentSparsityBetweenElements
  std::map<dof_id_type, std::vector<dof_id_type>> _elem_map;
  /// Whether the AugmentSparsityBetweenElements relationship manager has been added
  bool _rm_added;
};
