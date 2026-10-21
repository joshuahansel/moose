//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HSCoupler2D3DAC.h"
#include "HeatStructureCylindricalAC.h"
#include "HeatStructureFromFile3DAC.h"
#include "MeshAlignment2D3D.h"
#include "AugmentSparsityBetweenElements.h"
#include "RelationshipManager.h"
#include "MooseMesh.h"
#include "THMUtils.h"

// Registered under the classic Component's name ('HSCoupler2D3D'); see FlowChannel1PhaseAC.C for
// why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSCoupler2D3DAC,
                           "HSCoupler2D3D",
                           "add_user_object");
registerMooseActionAliased("ThermalHydraulicsApp", HSCoupler2D3DAC, "HSCoupler2D3D", "add_bc");
registerActionComponentAliased("ThermalHydraulicsApp", HSCoupler2D3DAC, "HSCoupler2D3D");

InputParameters
HSCoupler2D3DAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addRequiredParam<std::string>("heat_structure_2d", "The 2D heat structure to couple");
  params.addRequiredParam<std::string>("heat_structure_3d", "The 3D heat structure to couple");
  params.addRequiredParam<BoundaryName>("boundary_2d",
                                        "The boundary of the 2D heat structure to couple");
  params.addRequiredParam<BoundaryName>("boundary_3d",
                                        "The boundary of the 3D heat structure to couple");

  params.addParam<bool>("include_radiation", true, "Include radiation component of heat flux");
  params.addParam<FunctionName>(
      "emissivity_2d",
      "Emissivity of the 2D heat structure boundary as a function of temperature [K]");
  params.addParam<FunctionName>(
      "emissivity_3d",
      "Emissivity of the 3D heat structure boundary as a function of temperature [K]");
  params.addRequiredParam<FunctionName>("gap_thickness",
                                        "Gap thickness [m] as a function of temperature [K]");
  params.addRequiredParam<FunctionName>(
      "gap_thermal_conductivity",
      "Gap thermal conductivity [W/(m-K)] as a function of temperature [K]");
  params.addParam<FunctionName>(
      "gap_htc", 0, "Gap heat transfer coefficient [W/(m^2-K)] as a function of temperature [K]");
  params.addRangeCheckedParam<Real>("symmetry_factor",
                                    1.0,
                                    "symmetry_factor>=1.0",
                                    "Azimuthal symmetry correction factor (>= 1.0). Equal to 2*pi "
                                    "divided by the azimuthal angle covered by the 3D mesh.");

  params.addClassDescription("Couples a 2D heat structure boundary to a 3D heat structure boundary "
                             "using gap heat transfer.");
  return params;
}

HSCoupler2D3DAC::HSCoupler2D3DAC(const InputParameters & params)
  : ActionComponent(params),
    _hs_name_2d(getParam<std::string>("heat_structure_2d")),
    _hs_name_3d(getParam<std::string>("heat_structure_3d")),
    _boundary_2d(getParam<BoundaryName>("boundary_2d")),
    _boundary_3d(getParam<BoundaryName>("boundary_3d")),
    _temperature_2d_uo_name(genName(name(), "2d_uo")),
    _coupler_uo_name(genName(name(), "3d_uo")),
    _rm_added(false)
{
  THM::informExperimental(*this);
  addRequiredTask("add_user_object");
  addRequiredTask("add_bc");
}

const HeatStructureCylindricalAC &
HSCoupler2D3DAC::hs2D() const
{
  for (const auto & action : _awh.getActionListByName("list_component"))
    if (action->name() == _hs_name_2d)
      if (const auto * hs = dynamic_cast<const HeatStructureCylindricalAC *>(action))
        return *hs;
  paramError("heat_structure_2d",
             "There is no 2D cylindrical heat structure with the name '",
             _hs_name_2d,
             "'.");
}

const HeatStructureFromFile3DAC &
HSCoupler2D3DAC::hs3D() const
{
  for (const auto & action : _awh.getActionListByName("list_component"))
    if (action->name() == _hs_name_3d)
      if (const auto * hs = dynamic_cast<const HeatStructureFromFile3DAC *>(action))
        return *hs;
  paramError(
      "heat_structure_3d", "There is no 3D heat structure with the name '", _hs_name_3d, "'.");
}

void
HSCoupler2D3DAC::checkInputs()
{
  if (getParam<bool>("include_radiation"))
  {
    if (!(isParamValid("emissivity_2d") && isParamValid("emissivity_3d")))
      paramError("include_radiation",
                 "If 'include_radiation' is 'true', then 'emissivity_2d' and 'emissivity_3d' are "
                 "required.");
  }
  else if (isParamValid("emissivity_2d") || isParamValid("emissivity_3d"))
    paramError("include_radiation",
               "If 'include_radiation' is 'false', then neither 'emissivity_2d' nor "
               "'emissivity_3d' can be specified.");

  if (!hs2D().hasBoundary(_boundary_2d))
    paramError("boundary_2d",
               "The heat structure '",
               _hs_name_2d,
               "' does not have the boundary '",
               _boundary_2d,
               "'.");
  if (!hs3D().hasBoundary(_boundary_3d))
    paramError("boundary_3d",
               "The heat structure '",
               _hs_name_3d,
               "' does not have the boundary '",
               _boundary_3d,
               "'.");

  auto & alignment = meshAlignment();
  if (!alignment.meshesAreAligned())
    paramError("heat_structure_2d", "The meshes of the heat structures are not aligned.");

  const unsigned int needed_ad_container_size = 4 * alignment.getMaxCouplingSize() + 6;
  if (MOOSE_AD_MAX_DOFS_PER_ELEM < needed_ad_container_size)
    mooseError(name(),
               ": MOOSE must be configured with a larger AD container size (>= ",
               needed_ad_container_size,
               "). See HSCoupler2D3D's documentation for more information.");
}

MeshAlignment2D3D &
HSCoupler2D3DAC::meshAlignment()
{
  if (_mesh_alignment)
    return *_mesh_alignment;

  auto & mesh = getProblem().mesh();
  _mesh_alignment = std::make_unique<MeshAlignment2D3D>(mesh);

  std::array<std::vector<std::tuple<dof_id_type, unsigned short int>>, 2> boundary_infos;
  const std::array<BoundaryName, 2> boundaries = {_boundary_2d, _boundary_3d};
  for (const auto i : make_range(2u))
  {
    const auto boundary_id = mesh.getBoundaryID(boundaries[i]);
    for (const auto & [elem_id, side, bid] : mesh.buildSideList())
      if (bid == boundary_id)
        boundary_infos[i].emplace_back(elem_id, side);
  }

  const auto & hs_2d = hs2D();
  RealVectorValue axis(1, 0, 0);
  if (hs_2d.isParamValid("orientation"))
    axis = hs_2d.direction();
  _mesh_alignment->initialize(boundary_infos[0], boundary_infos[1], hs_2d.translation(), axis);

  // The coupled elements are not adjacent in the mesh, so their coupling must be declared
  if (_mesh_alignment->meshesAreAligned())
    for (const auto elem_id : _mesh_alignment->getSecondaryElemIDs())
      if (_mesh_alignment->hasCoupledPrimaryElemID(elem_id))
      {
        const auto coupled_elem_id = _mesh_alignment->getCoupledPrimaryElemID(elem_id);
        _elem_map[elem_id].push_back(coupled_elem_id);
        _elem_map[coupled_elem_id].push_back(elem_id);
      }

  return *_mesh_alignment;
}

void
HSCoupler2D3DAC::actOnAdditionalTasks()
{
  if (_current_task == "add_bc")
    addBCs();
}

void
HSCoupler2D3DAC::addUserObjects()
{
  // The inputs are checked on the first task that uses the heat structures
  checkInputs();

  // side user object on the 2D boundary to cache temperature values by element ID
  {
    const std::string class_name = "StoreVariableByElemIDSideUserObject";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::vector<BoundaryName>>("boundary") = {_boundary_2d};
    params.set<std::vector<VariableName>>("variable") = {
        hs2D().getPhysics().getParam<VariableName>("temperature_name")};
    params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};
    // This user object needs to execute before the one on the 3D boundary
    params.set<int>("execution_order_group") = -1;
    getProblem().addUserObject(class_name, _temperature_2d_uo_name, params);
  }

  // The radius of the 2D heat structure boundary
  const auto & hs_2d = hs2D();
  const auto radius_2d = hs_2d.innerRadius() + hs_2d.totalWidth();

  // side user object on the 3D boundary to compute the heat fluxes across each 3D boundary
  {
    const std::string class_name = "HSCoupler2D3DUserObject";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::vector<BoundaryName>>("boundary") = {_boundary_3d};
    params.set<std::vector<VariableName>>("temperature") = {
        hs3D().getPhysics().getParam<VariableName>("temperature_name")};
    params.set<Real>("radius_2d") = radius_2d;
    params.set<Real>("symmetry_factor") = getParam<Real>("symmetry_factor");
    if (getParam<bool>("include_radiation"))
    {
      params.set<FunctionName>("emissivity_2d") = getParam<FunctionName>("emissivity_2d");
      params.set<FunctionName>("emissivity_3d") = getParam<FunctionName>("emissivity_3d");
    }
    params.set<FunctionName>("gap_thickness") = getParam<FunctionName>("gap_thickness");
    params.set<FunctionName>("gap_thermal_conductivity") =
        getParam<FunctionName>("gap_thermal_conductivity");
    params.set<FunctionName>("gap_htc") = getParam<FunctionName>("gap_htc");
    params.set<UserObjectName>("temperature_2d_uo") = _temperature_2d_uo_name;
    params.set<MeshAlignment2D3D *>("mesh_alignment") = &meshAlignment();
    params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};
    getProblem().addUserObject(class_name, _coupler_uo_name, params);
  }
}

void
HSCoupler2D3DAC::addBCs()
{
  const std::array<std::pair<BoundaryName, const char *>, 2> boundaries = {
      std::make_pair(_boundary_2d, "2d"), std::make_pair(_boundary_3d, "3d")};
  const std::array<VariableName, 2> variables = {
      hs2D().getPhysics().getParam<VariableName>("temperature_name"),
      hs3D().getPhysics().getParam<VariableName>("temperature_name")};
  for (const auto i : make_range(2u))
  {
    const std::string class_name = "HSCoupler2D3DBC";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<NonlinearVariableName>("variable") = variables[i];
    params.set<std::vector<BoundaryName>>("boundary") = {boundaries[i].first};
    params.set<UserObjectName>("hs_coupler_2d3d_uo") = _coupler_uo_name;
    getProblem().addBoundaryCondition(
        class_name, genName(name(), class_name, boundaries[i].second), params);
  }
}

void
HSCoupler2D3DAC::addRelationshipManagers(Moose::RelationshipManagerType when_type)
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
