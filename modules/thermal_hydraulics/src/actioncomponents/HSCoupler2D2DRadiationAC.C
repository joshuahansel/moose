//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HSCoupler2D2DRadiationAC.h"
#include "HeatStructureCylindricalAC.h"
#include "MeshAlignment2D2D.h"
#include "AugmentSparsityBetweenElements.h"
#include "RelationshipManager.h"
#include "MooseMesh.h"
#include "MooseUtils.h"
#include "THMUtils.h"

#include <numeric>

// Registered under the classic Component's name ('HSCoupler2D2DRadiation'); see
// FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSCoupler2D2DRadiationAC,
                           "HSCoupler2D2DRadiation",
                           "add_user_object");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSCoupler2D2DRadiationAC,
                           "HSCoupler2D2DRadiation",
                           "add_bc");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HSCoupler2D2DRadiationAC,
                               "HSCoupler2D2DRadiation");

InputParameters
HSCoupler2D2DRadiationAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addRequiredParam<std::vector<std::string>>("heat_structures",
                                                    "The heat structures to couple");
  params.addRequiredParam<std::vector<BoundaryName>>(
      "boundaries", "The boundaries of the heat structures to couple");
  params.addRequiredRangeCheckedParam<std::vector<Real>>(
      "emissivities",
      "emissivities > 0 & emissivities < 1",
      "Emissivities of each heat structure surface");
  params.addRequiredRangeCheckedParam<std::vector<std::vector<Real>>>(
      "view_factors",
      "view_factors >= 0 & view_factors <= 1",
      "The view factors between each surface, as a matrix. The row/column ordering corresponds to "
      "the ordering in 'heat_structures', with an additional row and column if "
      "'include_environment' is set to 'true'. Each row must sum to one.");
  params.addRequiredParam<bool>(
      "include_environment",
      "Whether or not to include an environment surrounding all of the surfaces. If the heat "
      "structure surfaces themselves form an enclosure, then set this parameter to 'false'.");
  params.addParam<Real>(
      "T_environment",
      "Environment temperature [K]. Only set if 'include_environment' is set to true.");
  params.addClassDescription("Couples boundaries of multiple 2D heat structures via radiation");
  return params;
}

HSCoupler2D2DRadiationAC::HSCoupler2D2DRadiationAC(const InputParameters & params)
  : ActionComponent(params),
    _hs_names(getParam<std::vector<std::string>>("heat_structures")),
    _hs_boundaries(getParam<std::vector<BoundaryName>>("boundaries")),
    _include_environment(getParam<bool>("include_environment")),
    _n_hs(_hs_names.size()),
    _n_surfaces(_include_environment ? _n_hs + 1 : _n_hs),
    _temperature_uo_name(genName(name(), "temperature_uo")),
    _coupler_uo_name(genName(name(), "uo")),
    _rm_added(false)
{
  THM::informExperimental(*this);
  addRequiredTask("add_user_object");
  addRequiredTask("add_bc");
}

const HeatStructureBaseAC &
HSCoupler2D2DRadiationAC::hs(unsigned int i) const
{
  for (const auto & action : _awh.getActionListByName("list_component"))
    if (action->name() == _hs_names[i])
    {
      if (dynamic_cast<const HeatStructureCylindricalAC *>(action))
        return *dynamic_cast<const HeatStructureBaseAC *>(action);
      paramError("heat_structures",
                 "The component '",
                 _hs_names[i],
                 "' is not of type 'HeatStructureCylindricalBase'");
    }
  paramError("heat_structures", "The component '", _hs_names[i], "' does not exist");
}

void
HSCoupler2D2DRadiationAC::checkInputs()
{
  if (_hs_boundaries.size() != _n_hs)
    paramError("boundaries",
               "The number of entries in parameter 'boundaries' (",
               _hs_boundaries.size(),
               ") must equal the number of entries of parameter 'heat_structures' (",
               _n_hs,
               ")");
  if (getParam<std::vector<Real>>("emissivities").size() != _n_hs)
    paramError("emissivities",
               "The number of entries in parameter 'emissivities' (",
               getParam<std::vector<Real>>("emissivities").size(),
               ") must equal the number of entries of parameter 'heat_structures' (",
               _n_hs,
               ")");

  // check view factor matrix dimensions
  const auto & view_factors = getParam<std::vector<std::vector<Real>>>("view_factors");
  bool correct_size = view_factors.size() == _n_surfaces;
  if (correct_size)
    for (const auto & row : view_factors)
      if (row.size() != _n_surfaces)
        correct_size = false;
  if (!correct_size)
    paramError("view_factors",
               "The parameter 'view_factors' must be a square matrix of size ",
               _n_surfaces,
               ". For example, a size 2 matrix is provided as '0.2 0.8; 0.7 0.3'. The row/column "
               "ordering corresponds to the ordering in 'heat_structures', with an additional "
               "row/column if 'include_environment' is set to 'true'.");

  // check that all view factor matrix rows sum to one
  for (const auto & row : view_factors)
    if (!MooseUtils::absoluteFuzzyEqual(std::accumulate(row.begin(), row.end(), 0.0), 1.0))
      paramError("view_factors", "All rows in 'view_factors' must sum to one.");

  for (const auto i : index_range(_hs_names))
    if (!hs(i).hasBoundary(_hs_boundaries[i]))
      paramError("boundaries",
                 "The heat structure '",
                 _hs_names[i],
                 "' does not have the boundary '",
                 _hs_boundaries[i],
                 "'.");

  if (_n_hs > 1 && !meshAlignment().meshesAreAligned())
    paramError("heat_structures", "The meshes of the heat structures are not aligned.");
}

MeshAlignment2D2D &
HSCoupler2D2DRadiationAC::meshAlignment()
{
  if (_mesh_alignment)
    return *_mesh_alignment;

  auto & mesh = getProblem().mesh();
  _mesh_alignment = std::make_unique<MeshAlignment2D2D>(mesh);

  // With a single heat structure there is nothing to align
  if (_n_hs < 2)
    return *_mesh_alignment;

  std::vector<std::vector<std::tuple<dof_id_type, unsigned short int>>> boundary_infos(_n_hs);
  for (const auto i : index_range(_hs_names))
  {
    const auto boundary_id = mesh.getBoundaryID(_hs_boundaries[i]);
    for (const auto & [elem_id, side, bid] : mesh.buildSideList())
      if (bid == boundary_id)
        boundary_infos[i].emplace_back(elem_id, side);
  }
  const auto & first = hs(0);
  RealVectorValue axis(1, 0, 0);
  if (first.isParamValid("orientation"))
    axis = first.direction();
  _mesh_alignment->initialize(boundary_infos, first.translation(), axis);

  // The coupled elements are not adjacent in the mesh, so their coupling must be declared
  if (_mesh_alignment->meshesAreAligned())
    for (const auto primary_elem_id : _mesh_alignment->getPrimaryElemIDs())
    {
      std::vector<dof_id_type> elem_ids = {primary_elem_id};
      const auto & coupled_elem_ids = _mesh_alignment->getCoupledSecondaryElemIDs(primary_elem_id);
      elem_ids.insert(elem_ids.end(), coupled_elem_ids.begin(), coupled_elem_ids.end());
      for (const auto i : index_range(elem_ids))
        for (const auto j : make_range(i + 1, elem_ids.size()))
        {
          _elem_map[elem_ids[i]].push_back(elem_ids[j]);
          _elem_map[elem_ids[j]].push_back(elem_ids[i]);
        }
    }

  return *_mesh_alignment;
}

void
HSCoupler2D2DRadiationAC::actOnAdditionalTasks()
{
  if (_current_task == "add_bc")
    addBCs();
}

void
HSCoupler2D2DRadiationAC::addUserObjects()
{
  // The inputs are checked on the first task that uses the heat structures
  checkInputs();

  const auto temperature_name = hs(0).getPhysics().getParam<VariableName>("temperature_name");

  // side user object on the boundaries to cache temperature values by element ID
  {
    const std::string class_name = "StoreVariableByElemIDSideUserObject";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::vector<BoundaryName>>("boundary") = _hs_boundaries;
    params.set<std::vector<VariableName>>("variable") = {temperature_name};
    params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};
    // This user object needs to execute before the one below
    params.set<int>("execution_order_group") = -1;
    getProblem().addUserObject(class_name, _temperature_uo_name, params);
  }

  // side user object to compute the heat fluxes across each boundary
  {
    const std::string class_name = "HSCoupler2D2DRadiationUserObject";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::vector<BoundaryName>>("boundary") = {_hs_boundaries[0]};
    params.applySpecificParameters(
        parameters(), {"emissivities", "view_factors", "include_environment", "T_environment"});
    params.set<UserObjectName>("temperature_uo") = _temperature_uo_name;
    params.set<MeshAlignment2D2D *>("mesh_alignment") = &meshAlignment();
    params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};
    getProblem().addUserObject(class_name, _coupler_uo_name, params);
  }
}

void
HSCoupler2D2DRadiationAC::addBCs()
{
  for (const auto i : index_range(_hs_names))
  {
    const std::string class_name = "HSCoupler2D2DRadiationBC";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<NonlinearVariableName>("variable") =
        hs(i).getPhysics().getParam<VariableName>("temperature_name");
    params.set<std::vector<BoundaryName>>("boundary") = {_hs_boundaries[i]};
    params.set<UserObjectName>("hs_coupler_2d2d_uo") = _coupler_uo_name;
    getProblem().addBoundaryCondition(class_name, genName(name(), class_name, i), params);
  }
}

void
HSCoupler2D2DRadiationAC::addRelationshipManagers(Moose::RelationshipManagerType when_type)
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
