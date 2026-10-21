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
#include "MooseEnum.h"

class FlowChannel1PhaseAC;
class HeatStructureFromFile3DAC;
class MeshAlignment1D3D;

/**
 * Connects one or more FlowChannel1PhaseAC to a boundary of a HeatStructureFromFile3DAC, the
 * native-ActionComponent counterpart of HeatTransferFromHeatStructure3D1Phase.
 *
 * The wall temperature seen by each flow channel is the average of the heat structure temperature
 * over the layer of the boundary next to each flow channel element (NearestPointLayeredSideAverage,
 * with one point per flow channel), and the heat transfer to the heat structure is computed per
 * boundary face from the flow channel element it is paired with, through a MeshAlignment1D3D. The
 * flow channels must all be aligned with the same coordinate axis and have the same number of
 * elements and the same length.
 *
 * The elements of the flow channels are not adjacent to the boundary elements of the heat
 * structure, so the Jacobian needs an explicit cross-element coupling, added with
 * AugmentSparsityBetweenElements (as for HeatTransferFromHeatStructure1PhaseAC).
 */
class HeatTransferFromHeatStructure3D1PhaseAC : public HeatTransferFrom1PhaseAC
{
public:
  static InputParameters validParams();
  HeatTransferFromHeatStructure3D1PhaseAC(const InputParameters & params);

protected:
  virtual bool isTemperatureType() const override { return true; }

  virtual void actOnAdditionalTasks() override;
  virtual void addUserObjects() override;
  virtual void addRelationshipManagers(Moose::RelationshipManagerType when_type) override;
  virtual std::vector<const FlowChannel1PhasePhysics *> connectedPhysicsList() const override;

  virtual void addHeatTransferVariables() override;
  virtual void addHeatTransferICs() override;
  virtual void addHeatTransferKernels() override;
  virtual void addHeatTransferAuxKernels() override;

private:
  /// The flow channel components, which must all be flow channels
  std::vector<const FlowChannel1PhaseAC *> flowChannels() const;
  /// The heat structure, which must be a heat structure loaded from a file
  const HeatStructureFromFile3DAC & hs() const;
  /// Checks the flow channels and the heat structure, and determines the axis the flow channels
  /// are aligned with
  void checkInputs();
  /// Builds the alignment between the flow channels and the heat structure boundary if it has not
  /// been built yet
  MeshAlignment1D3D & meshAlignment();

  /// Names of the flow channels
  const std::vector<std::string> & _flow_channel_names;
  /// Name of the heat structure boundary
  const BoundaryName & _boundary;
  /// Axis the flow channels are aligned with
  MooseEnum _layered_average_uo_direction;
  /// Number of layers (elements) along the flow channels
  unsigned int _num_layers;
  /// Name of the user object averaging the heat structure temperature over the layers
  const UserObjectName _T_wall_avg_uo_name;
  /// Name of the user object coupling the flow channels and the heat structure
  const UserObjectName _heat_transfer_uo_name;

  /// Alignment of the flow channels and the heat structure boundary
  std::unique_ptr<MeshAlignment1D3D> _mesh_alignment;
  /// Map of element ID -> coupled element IDs, consumed by AugmentSparsityBetweenElements
  std::map<dof_id_type, std::vector<dof_id_type>> _elem_map;
  /// Whether the AugmentSparsityBetweenElements relationship manager has been added
  bool _rm_added;
};
