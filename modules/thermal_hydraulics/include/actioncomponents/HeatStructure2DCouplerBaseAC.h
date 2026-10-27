//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ActionComponentConvergenceInterface.h"
#include "ActionComponent.h"
#include "NamingInterface.h"
#include "HeatStructureBaseAC.h"

class MeshAlignment;

/**
 * Base class for ActionComponents that couple a boundary of one 2D heat structure ActionComponent
 * (the primary) to a boundary of another (the secondary), the native-ActionComponent counterpart of
 * HeatStructure2DCouplerBase.
 *
 * The two boundaries are paired with a MeshAlignment. The elements of the two heat structures are
 * not adjacent in the combined mesh, so the Jacobian needs an explicit cross-element coupling,
 * added with AugmentSparsityBetweenElements (as for JunctionOneToOne1PhaseAC).
 *
 * Classic THM applies the coupling to a cylindrical heat structure with RZ-specific boundary
 * conditions that multiply by the local circumference. HeatStructureCylindricalAC uses MOOSE's
 * native COORD_RZ coordinate system, which applies that factor itself, so the plain boundary
 * conditions are used for both geometries.
 */
class HeatStructure2DCouplerBaseAC : public virtual ActionComponent,
                                     public NamingInterface,
                                     public NoNonlinearConvergenceAC
{
public:
  static InputParameters validParams();
  HeatStructure2DCouplerBaseAC(const InputParameters & params);

protected:
  virtual void actOnAdditionalTasks() override;
  virtual void addRelationshipManagers(Moose::RelationshipManagerType when_type) override;

  /// Adds the boundary conditions, once the checks on the heat structures have been done
  virtual void addBCs() = 0;
  /// Additional checks of the derived class, done after the checks that the heat structures and
  /// their boundaries exist but before the check that the heat structures have the same type, so
  /// that a derived class that restricts the type reports its own message
  virtual void checkCoupling() {}

  /// The i-th heat structure (0 is the primary and 1 the secondary), which must be a 2D heat
  /// structure. The lookup is done lazily so that the heat structures may be declared after this
  /// component.
  const HeatStructureBaseAC & hs(unsigned int i) const;
  /// Builds the alignment between the two boundaries if it has not been built yet
  MeshAlignment & meshAlignment();
  /// Whether the i-th heat structure is cylindrical
  bool isCylindrical(unsigned int i) const;

  /// Names of the two coupled heat structures
  const std::array<std::string, 2> _hs_names;
  /// Names of the boundaries of the two heat structures that are coupled
  const std::array<BoundaryName, 2> _hs_boundaries;
  /// Area of each boundary
  std::array<Real, 2> _areas;
  /// Fraction of each boundary that couples to the other, which is not one for the smaller of the
  /// two boundaries if the areas are not equal
  std::array<Real, 2> _coupling_area_fractions;

private:
  /// Checks the heat structures and their boundaries, and sets the areas and area fractions
  void checkHeatStructures();

  /// Alignment of the two boundaries
  std::unique_ptr<MeshAlignment> _mesh_alignment;
  /// Map of element ID -> coupled element IDs, consumed by AugmentSparsityBetweenElements
  std::map<dof_id_type, std::vector<dof_id_type>> _elem_map;
  /// Whether the AugmentSparsityBetweenElements relationship manager has been added
  bool _rm_added;
};
