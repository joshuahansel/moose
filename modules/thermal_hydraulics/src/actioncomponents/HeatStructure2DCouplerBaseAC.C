//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatStructure2DCouplerBaseAC.h"
#include "HeatStructureCylindricalAC.h"
#include "HeatStructurePlateAC.h"
#include "MeshAlignment.h"
#include "AugmentSparsityBetweenElements.h"
#include "RelationshipManager.h"
#include "MooseMesh.h"
#include "THMUtils.h"

InputParameters
HeatStructure2DCouplerBaseAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addRequiredParam<std::string>("primary_heat_structure",
                                       "The first heat structure to couple");
  params.addRequiredParam<std::string>("secondary_heat_structure",
                                       "The second heat structure to couple");
  params.addRequiredParam<BoundaryName>("primary_boundary",
                                        "The boundary of the first heat structure to couple");
  params.addRequiredParam<BoundaryName>("secondary_boundary",
                                        "The boundary of the second heat structure to couple");
  return params;
}

HeatStructure2DCouplerBaseAC::HeatStructure2DCouplerBaseAC(const InputParameters & params)
  : ActionComponent(params),
    _hs_names({getParam<std::string>("primary_heat_structure"),
               getParam<std::string>("secondary_heat_structure")}),
    _hs_boundaries(
        {getParam<BoundaryName>("primary_boundary"), getParam<BoundaryName>("secondary_boundary")}),
    _areas({0.0, 0.0}),
    _coupling_area_fractions({1.0, 1.0}),
    _rm_added(false)
{
  THM::informExperimental(*this);
  addRequiredTask("add_bc");
}

const HeatStructureBaseAC &
HeatStructure2DCouplerBaseAC::hs(unsigned int i) const
{
  const auto param = i == 0 ? "primary_heat_structure" : "secondary_heat_structure";
  for (const auto & action : _awh.getActionListByName("list_component"))
    if (action->name() == _hs_names[i])
    {
      if (const auto * hs = dynamic_cast<const HeatStructureBaseAC *>(action))
        return *hs;
      paramError(param, "The component '", _hs_names[i], "' is not of type 'HeatStructureBase'");
    }
  paramError(param, "The component '", _hs_names[i], "' does not exist");
}

bool
HeatStructure2DCouplerBaseAC::isCylindrical(unsigned int i) const
{
  return dynamic_cast<const HeatStructureCylindricalAC *>(&hs(i)) != nullptr;
}

void
HeatStructure2DCouplerBaseAC::checkHeatStructures()
{
  for (const auto i : make_range(2u))
    if (!hs(i).hasBoundary(_hs_boundaries[i]))
      paramError(i == 0 ? "primary_boundary" : "secondary_boundary",
                 "The heat structure '",
                 _hs_names[i],
                 "' does not have the boundary '",
                 _hs_boundaries[i],
                 "'.");

  checkCoupling();

  if (isCylindrical(0) != isCylindrical(1))
    paramError("primary_heat_structure", "The coupled heat structures must have the same type.");

  for (const auto i : make_range(2u))
    _areas[i] = hs(i).getBoundaryArea(_hs_boundaries[i]);

  // Compute the coupling area fractions in case areas are not equal
  if (_areas[0] > _areas[1])
    _coupling_area_fractions = {_areas[1] / _areas[0], 1.0};
  else
    _coupling_area_fractions = {1.0, _areas[0] / _areas[1]};
}

MeshAlignment &
HeatStructure2DCouplerBaseAC::meshAlignment()
{
  if (_mesh_alignment)
    return *_mesh_alignment;

  auto & mesh = getProblem().mesh();
  _mesh_alignment = std::make_unique<MeshAlignment>(mesh);

  std::array<std::vector<std::tuple<dof_id_type, unsigned short int>>, 2> boundary_infos;
  for (const auto i : make_range(2u))
  {
    const auto boundary_id = mesh.getBoundaryID(_hs_boundaries[i]);
    for (const auto & [elem_id, side, bid] : mesh.buildSideList())
      if (bid == boundary_id)
        boundary_infos[i].emplace_back(elem_id, side);
  }
  _mesh_alignment->initialize(boundary_infos[0], boundary_infos[1]);

  // The coupled elements are not adjacent in the mesh, so their coupling must be declared
  if (_mesh_alignment->meshesAreAligned())
    for (const auto & elem_id : _mesh_alignment->getPrimaryElemIDs())
      if (_mesh_alignment->hasCoupledElemID(elem_id))
      {
        const auto coupled_elem_id = _mesh_alignment->getCoupledElemID(elem_id);
        _elem_map[elem_id].push_back(coupled_elem_id);
        _elem_map[coupled_elem_id].push_back(elem_id);
      }

  return *_mesh_alignment;
}

void
HeatStructure2DCouplerBaseAC::actOnAdditionalTasks()
{
  if (_current_task != "add_bc")
    return;

  checkHeatStructures();
  addBCs();
}

void
HeatStructure2DCouplerBaseAC::addRelationshipManagers(Moose::RelationshipManagerType when_type)
{
  // This virtual is called once per RelationshipManagerType category (geometric, algebraic,
  // coupling); add the (single) RM only once, on the first call. _elem_map is populated later, when
  // the mesh alignment is built - the RM holds a reference to it, consulted only once mesh
  // distribution actually needs it, well after the mesh (and thus _elem_map) exist.
  if (_rm_added || when_type != Moose::RelationshipManagerType::GEOMETRIC)
    return;
  _rm_added = true;

  const std::string class_name = "AugmentSparsityBetweenElements";
  InputParameters params = _factory.getValidParams(class_name);
  params.set<Moose::RelationshipManagerType>("rm_type") =
      Moose::RelationshipManagerType::COUPLING | Moose::RelationshipManagerType::ALGEBRAIC |
      Moose::RelationshipManagerType::GEOMETRIC;
  params.set<std::string>("for_whom") = name();
  params.set<MooseMesh *>("mesh") = _awh.getMesh().get();
  params.set<std::map<dof_id_type, std::vector<dof_id_type>> *>("_elem_map") = &_elem_map;
  auto rm =
      _factory.create<RelationshipManager>(class_name, genName(name(), "sparsity_rm"), params);
  if (!_app.addRelationshipManager(rm))
    _factory.releaseSharedObjects(*rm);
}
