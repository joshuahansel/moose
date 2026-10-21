//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HeatTransferFrom1PhaseAC.h"
#include "HeatStructureBaseAC.h"

class MeshAlignment;

/**
 * Connects a FlowChannel1PhaseAC to a side of a 2D heat structure ActionComponent, the
 * native-ActionComponent counterpart of HeatTransferFromHeatStructure1Phase.
 *
 * The wall temperature seen by the flow channel is the temperature of the heat structure, and the
 * heat flux computed from it is applied to both: as a source in the flow channel energy equation
 * and as a boundary condition on the heat structure. The two meshes are coupled through a
 * MeshAlignment, which pairs the flow channel elements with the elements next to the heat structure
 * boundary and requires the two to be offset by the same translation. The elements are not
 * adjacent in the combined mesh, so the Jacobian needs an explicit cross-element coupling, added
 * with AugmentSparsityBetweenElements (as for JunctionOneToOne1PhaseAC).
 *
 * Classic THM integrates the heat flux boundary condition of a cylindrical heat structure on a
 * Cartesian mesh and multiplies by the circumference. HeatStructureCylindricalAC uses MOOSE's
 * native COORD_RZ coordinate system, which applies that factor during integration, so the boundary
 * condition is told that the coordinate system is not cylindrical for both geometries.
 */
class HeatTransferFromHeatStructure1PhaseAC : public HeatTransferFrom1PhaseAC
{
public:
  static InputParameters validParams();
  HeatTransferFromHeatStructure1PhaseAC(const InputParameters & params);

protected:
  virtual bool isTemperatureType() const override { return true; }

  virtual void actOnAdditionalTasks() override;
  virtual void addUserObjects() override;
  virtual void addRelationshipManagers(Moose::RelationshipManagerType when_type) override;

  virtual void addHeatTransferVariables() override;
  virtual void addHeatTransferICs() override;
  virtual void addHeatTransferMaterials() override;
  virtual void addHeatTransferKernels() override;
  virtual void addHeatTransferAuxKernels() override;

private:
  /// Looks up the heat structure by name, which must be a 2D heat structure
  const HeatStructureBaseAC & hs() const;
  /// Looks up the component named by \p param, which must be of type \p T (named \p type_name in
  /// the error message)
  template <typename T>
  const T & component(const std::string & param, const std::string & type_name) const;

  /// Checks the heat structure and the side of it that is connected, and determines the name and
  /// type of the connected boundary
  void checkHeatStructure();
  /// Adds the boundary condition on the heat structure
  void addHeatStructureBC();
  /// Builds the alignment between the flow channel and the heat structure boundary if it has not
  /// been built yet
  MeshAlignment & meshAlignment();

  /// Name of the heat structure boundary that is connected
  BoundaryName _hs_boundary;
  /// Type of the heat structure boundary that is connected
  HeatStructureBaseAC::ExternalBoundaryType _hs_boundary_type;
  /// Name of the user object computing the heat flux
  const UserObjectName _heat_flux_uo_name;
  /// Name of the material property holding the heat structure temperature on the flow channel
  const MaterialPropertyName _T_wall_coupled_name;

  /// Alignment of the flow channel and the heat structure boundary
  std::unique_ptr<MeshAlignment> _mesh_alignment;
  /// Map of element ID -> coupled element IDs, consumed by AugmentSparsityBetweenElements
  std::map<dof_id_type, std::vector<dof_id_type>> _elem_map;
  /// Whether the AugmentSparsityBetweenElements relationship manager has been added
  bool _rm_added;
};

template <typename T>
const T &
HeatTransferFromHeatStructure1PhaseAC::component(const std::string & param,
                                                 const std::string & type_name) const
{
  const auto & comp_name = getParam<std::string>(param);
  for (const auto & action : _awh.getActionListByName("list_component"))
    if (action->name() == comp_name)
    {
      if (const auto * comp = dynamic_cast<const T *>(action))
        return *comp;
      paramError(param, "The component '", comp_name, "' is not of type '", type_name, "'");
    }
  paramError(param, "The component '", comp_name, "' does not exist");
}
