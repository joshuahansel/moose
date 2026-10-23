//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatTransferFromHeatStructure1PhaseAC.h"
#include "FlowChannel1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"
#include "FlowModelSinglePhase.h"
#include "MeshAlignment.h"
#include "AugmentSparsityBetweenElements.h"
#include "RelationshipManager.h"
#include "MooseMesh.h"
#include "THMNames.h"
#include "THMVariableCoordinator.h"
#include "MooseUtils.h"

// Registered under the classic Component's name ('HeatTransferFromHeatStructure1Phase'); see
// FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatTransferFromHeatStructure1PhaseAC,
                           "HeatTransferFromHeatStructure1Phase",
                           "add_user_object");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatTransferFromHeatStructure1PhaseAC,
                           "HeatTransferFromHeatStructure1Phase",
                           "add_bc");
registerHeatTransferFrom1PhaseACTasks("ThermalHydraulicsApp",
                                      HeatTransferFromHeatStructure1PhaseAC,
                                      "HeatTransferFromHeatStructure1Phase");

InputParameters
HeatTransferFromHeatStructure1PhaseAC::validParams()
{
  InputParameters params = HeatTransferFrom1PhaseAC::validParams();
  params.addRequiredParam<std::string>("hs", "Heat structure name");
  params.addParam<MooseEnum>("hs_side",
                             MooseEnum("inner outer start end"),
                             "Heat structure side. Exactly one of 'hs_side' and 'hs_boundary' must "
                             "be provided.");
  params.addParam<BoundaryName>("hs_boundary",
                                "Name of the coupled heat structure boundary. Exactly one of "
                                "'hs_side' and 'hs_boundary' must be provided.");
  params.addParam<MooseFunctorName>("scale", 1.0, "Functor by which to scale the heat flux");
  params.addClassDescription("Connects a 1-phase flow channel and a heat structure");
  return params;
}

HeatTransferFromHeatStructure1PhaseAC::HeatTransferFromHeatStructure1PhaseAC(
    const InputParameters & params)
  : ActionComponent(params),
    HeatTransferFrom1PhaseAC(params),
    _hs_boundary_type(HeatStructureBaseAC::ExternalBoundaryType::OUTER),
    _heat_flux_uo_name(genName(name(), "heat_flux_uo")),
    _rm_added(false)
{
  if (isParamValid("hs_side") == isParamValid("hs_boundary"))
    paramError("hs_side", "Exactly one of 'hs_side' and 'hs_boundary' must be provided.");

  addRequiredTask("add_user_object");
  addRequiredTask("add_bc");
}

const HeatStructureBaseAC &
HeatTransferFromHeatStructure1PhaseAC::hs() const
{
  return component<HeatStructureBaseAC>("hs", "HeatStructureBase");
}

void
HeatTransferFromHeatStructure1PhaseAC::checkHeatStructure()
{
  // Check the flow channel here as well, so a wrong name or type is reported as in the classic
  // component
  component<FlowChannel1PhaseAC>("flow_channel", "FlowChannelBase");

  const auto & heat_structure = hs();
  const auto & hs_name = getParam<std::string>("hs");

  if (isParamValid("hs_side"))
  {
    const std::string side = getParam<MooseEnum>("hs_side");
    _hs_boundary = genName(hs_name, side);
  }
  else
    _hs_boundary = getParam<BoundaryName>("hs_boundary");

  if (!heat_structure.hasExternalBoundary(_hs_boundary))
    paramError(isParamValid("hs_side") ? "hs_side" : "hs_boundary",
               "The heat structure boundary '",
               _hs_boundary,
               "' does not exist or is not an external boundary of the heat structure '",
               hs_name,
               "'.");
  _hs_boundary_type = heat_structure.getExternalBoundaryType(_hs_boundary);

  const Real P_hs = heat_structure.getUnitPerimeter(_hs_boundary_type);
  if (MooseUtils::absoluteFuzzyEqual(P_hs, 0.))
  {
    if (isParamValid("hs_side"))
      paramError("hs_side",
                 "'hs_side' parameter is set to '",
                 getParam<MooseEnum>("hs_side"),
                 "', but this side of the heat structure '",
                 hs_name,
                 "' has radius of zero.");
    else
      paramError("hs_boundary",
                 "The specified boundary '",
                 _hs_boundary,
                 "' of the heat structure '",
                 hs_name,
                 "' has a radius of zero.");
  }
  if (std::isnan(P_hs))
    paramError(isParamValid("hs_side") ? "hs_side" : "hs_boundary",
               "The specified boundary '",
               _hs_boundary,
               "' of the heat structure '",
               hs_name,
               "' is either START and END boundary, which may not be used.");
}

MeshAlignment &
HeatTransferFromHeatStructure1PhaseAC::meshAlignment()
{
  if (_mesh_alignment)
    return *_mesh_alignment;

  auto & mesh = getProblem().mesh();
  _mesh_alignment = std::make_unique<MeshAlignment>(mesh, true);

  // Elements of the flow channel
  const auto flow_channel_subdomain_ids =
      mesh.getSubdomainIDs(connectedPhysics().getClosuresBlocks());
  const std::set<SubdomainID> flow_channel_subdomains(flow_channel_subdomain_ids.begin(),
                                                      flow_channel_subdomain_ids.end());
  std::vector<dof_id_type> flow_channel_elem_ids;
  for (const auto * elem : mesh.getMesh().active_element_ptr_range())
    if (flow_channel_subdomains.count(elem->subdomain_id()))
      flow_channel_elem_ids.push_back(elem->id());

  // Element sides on the heat structure boundary
  const auto hs_boundary_id = mesh.getBoundaryID(_hs_boundary);
  std::vector<std::tuple<dof_id_type, unsigned short int>> hs_boundary_info;
  for (const auto & [elem_id, side, boundary_id] : mesh.buildSideList())
    if (boundary_id == hs_boundary_id)
      hs_boundary_info.emplace_back(elem_id, side);

  _mesh_alignment->initialize(flow_channel_elem_ids, hs_boundary_info);

  if (!_mesh_alignment->meshesAreAligned())
    paramError("flow_channel",
               "The centers of the elements of flow channel '",
               getParam<std::string>("flow_channel"),
               "' do not align with the centers of the specified heat structure side.");

  // The coupled elements are not adjacent in the mesh, so their coupling must be declared
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
HeatTransferFromHeatStructure1PhaseAC::actOnAdditionalTasks()
{
  // The inputs that refer to other components are checked before the base class uses them
  if (_current_task == "THM:add_heat_transfers")
    checkHeatStructure();

  HeatTransferFrom1PhaseAC::actOnAdditionalTasks();

  if (_current_task == "add_bc")
    addHeatStructureBC();
}

void
HeatTransferFromHeatStructure1PhaseAC::addHeatTransferVariables()
{
  THMVariableCoordinator::findOrCreate(_awh, _action_factory)
      .requestVariable(false,
                       wallTemperatureName(),
                       "LAGRANGE",
                       "FIRST",
                       1.0,
                       connectedPhysics().getClosuresBlocks());
}

void
HeatTransferFromHeatStructure1PhaseAC::addHeatTransferICs()
{
  // The heat structure has no initial temperature when the simulation is restarted
  if (!hs().isParamValid("initial_T"))
    return;

  const std::string class_name = "FunctionIC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<VariableName>("variable") = wallTemperatureName();
  params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
  params.set<FunctionName>("function") = hs().getParam<FunctionName>("initial_T");
  getProblem().addInitialCondition(class_name, genName(name(), "T_wall_ic"), params);
}

void
HeatTransferFromHeatStructure1PhaseAC::addUserObjects()
{
  ExecFlagEnum execute_on(MooseUtils::getDefaultExecFlagEnum());
  execute_on = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};

  const std::string class_name = "ADHeatFluxFromHeatStructure3EqnUserObject";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
  params.set<MeshAlignment *>("_mesh_alignment") = &meshAlignment();
  params.set<MaterialPropertyName>("T_wall") = coupledWallTemperatureName();
  params.set<std::vector<VariableName>>("P_hf") = {heatedPerimeterName()};
  params.set<MaterialPropertyName>("Hw") = wallHeatTransferCoefficientName();
  params.set<MaterialPropertyName>("T") = FlowModelSinglePhase::TEMPERATURE;
  params.set<MooseFunctorName>("scale") = getParam<MooseFunctorName>("scale");
  params.set<ExecFlagEnum>("execute_on") = execute_on;
  getProblem().addUserObject(class_name, _heat_flux_uo_name, params);
}

void
HeatTransferFromHeatStructure1PhaseAC::addHeatTransferKernels()
{
  const std::string class_name = "ADOneD3EqnEnergyHeatFlux";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
  params.set<NonlinearVariableName>("variable") = THM::RHOEA;
  params.set<UserObjectName>("q_uo") = _heat_flux_uo_name;
  getProblem().addKernel(class_name, genName(name(), "heat_flux_kernel"), params);
}

void
HeatTransferFromHeatStructure1PhaseAC::addHeatStructureBC()
{
  const auto & heat_structure = hs();

  const std::string class_name = "ADHeatFlux3EqnBC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<BoundaryName>>("boundary") = {_hs_boundary};
  params.set<NonlinearVariableName>("variable") =
      heat_structure.getPhysics().getParam<VariableName>("temperature_name");
  params.set<UserObjectName>("q_uo") = _heat_flux_uo_name;
  params.set<Real>("P_hs_unit") = heat_structure.getUnitPerimeter(_hs_boundary_type);
  params.set<unsigned int>("n_unit") = heat_structure.numberOfUnits();
  // See the class documentation
  params.set<bool>("hs_coord_system_is_cylindrical") = false;
  getProblem().addBoundaryCondition(class_name, genName(name(), "heat_flux_bc"), params);
}

void
HeatTransferFromHeatStructure1PhaseAC::addHeatTransferMaterials()
{
  // Transfer the temperature of the solid onto the flow channel
  const std::string class_name = "MeshAlignmentVariableTransferMaterial";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
  params.set<MaterialPropertyName>("property_name") = coupledWallTemperatureName();
  params.set<std::string>("paired_variable") =
      hs().getPhysics().getParam<VariableName>("temperature_name");
  params.set<MeshAlignment *>("_mesh_alignment") = &meshAlignment();
  getProblem().addMaterial(class_name, genName(name(), "T_wall_transfer_mat"), params);
}

void
HeatTransferFromHeatStructure1PhaseAC::addHeatTransferAuxKernels()
{
  // Transfer the temperature of the solid onto the flow channel as an auxiliary variable for
  // visualization
  const std::string class_name = "VariableValueTransferAux";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<AuxVariableName>("variable") = wallTemperatureName();
  params.set<std::vector<BoundaryName>>("boundary") = {getParam<std::string>("flow_channel")};
  params.set<BoundaryName>("paired_boundary") = _hs_boundary;
  params.set<std::vector<VariableName>>("paired_variable") = {
      hs().getPhysics().getParam<VariableName>("temperature_name")};
  getProblem().addAuxKernel(class_name, genName(name(), "T_wall_transfer"), params);
}

void
HeatTransferFromHeatStructure1PhaseAC::addRelationshipManagers(
    Moose::RelationshipManagerType when_type)
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
