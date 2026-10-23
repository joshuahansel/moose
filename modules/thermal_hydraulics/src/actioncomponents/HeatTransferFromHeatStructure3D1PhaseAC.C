//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatTransferFromHeatStructure3D1PhaseAC.h"
#include "FlowChannel1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"
#include "HeatStructureFromFile3DAC.h"
#include "FlowModelSinglePhase.h"
#include "MeshAlignment1D3D.h"
#include "AugmentSparsityBetweenElements.h"
#include "RelationshipManager.h"
#include "MooseMesh.h"
#include "THMNames.h"
#include "THMVariableCoordinator.h"
#include "MooseUtils.h"

// Registered under the classic Component's name ('HeatTransferFromHeatStructure3D1Phase'); see
// FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatTransferFromHeatStructure3D1PhaseAC,
                           "HeatTransferFromHeatStructure3D1Phase",
                           "add_user_object");
registerHeatTransferFrom1PhaseACTasks("ThermalHydraulicsApp",
                                      HeatTransferFromHeatStructure3D1PhaseAC,
                                      "HeatTransferFromHeatStructure3D1Phase");

InputParameters
HeatTransferFromHeatStructure3D1PhaseAC::validParams()
{
  InputParameters params = HeatTransferFrom1PhaseAC::validParams();
  // The flow channels are given by 'flow_channels' instead
  params.suppressParameter<std::string>("flow_channel");
  params.makeParamNotRequired<std::string>("flow_channel");
  params.set<std::string>("flow_channel", true) = "";
  params.addRequiredParam<std::vector<std::string>>(
      "flow_channels", "List of flow channel component names to connect to");
  params.addRequiredParam<BoundaryName>(
      "boundary", "The name of the heat structure boundary this heat transfer is applied on.");
  params.addRequiredParam<std::string>("hs", "Heat structure name");
  params.makeParamRequired<FunctionName>("P_hf");
  params.addClassDescription("Connects multiple 1-phase flow channels and a 3D heat structure");
  return params;
}

HeatTransferFromHeatStructure3D1PhaseAC::HeatTransferFromHeatStructure3D1PhaseAC(
    const InputParameters & params)
  : ActionComponent(params),
    HeatTransferFrom1PhaseAC(params),
    _flow_channel_names(getParam<std::vector<std::string>>("flow_channels")),
    _boundary(getParam<BoundaryName>("boundary")),
    _layered_average_uo_direction(MooseEnum("x y z")),
    _num_layers(0),
    _T_wall_avg_uo_name(genName(name(), "T_wall_avg_uo")),
    _heat_transfer_uo_name(genName(name(), "heat_flux_uo")),
    _rm_added(false)
{
  addRequiredTask("add_user_object");
}

std::vector<const FlowChannel1PhaseAC *>
HeatTransferFromHeatStructure3D1PhaseAC::flowChannels() const
{
  std::vector<const FlowChannel1PhaseAC *> flow_channels;
  for (const auto & fch_name : _flow_channel_names)
  {
    const FlowChannel1PhaseAC * fch = nullptr;
    for (const auto & action : _awh.getActionListByName("list_component"))
      if (action->name() == fch_name)
      {
        fch = dynamic_cast<const FlowChannel1PhaseAC *>(action);
        if (!fch)
          paramError("flow_channels",
                     "The component '",
                     fch_name,
                     "' is not of type 'FlowChannel1Phase'.");
      }
    if (!fch)
      paramError("flow_channels", "The component '", fch_name, "' does not exist.");
    flow_channels.push_back(fch);
  }
  return flow_channels;
}

std::vector<const FlowChannel1PhasePhysics *>
HeatTransferFromHeatStructure3D1PhaseAC::connectedPhysicsList() const
{
  std::vector<const FlowChannel1PhasePhysics *> physics;
  for (const auto * fch : flowChannels())
    physics.push_back(&fch->getPhysics());
  return physics;
}

const HeatStructureFromFile3DAC &
HeatTransferFromHeatStructure3D1PhaseAC::hs() const
{
  const auto & hs_name = getParam<std::string>("hs");
  for (const auto & action : _awh.getActionListByName("list_component"))
    if (action->name() == hs_name)
    {
      if (const auto * hs = dynamic_cast<const HeatStructureFromFile3DAC *>(action))
        return *hs;
      paramError("hs", "The component '", hs_name, "' is not a HeatStructureFromFile3D component.");
    }
  paramError("hs", "The component '", hs_name, "' does not exist.");
}

void
HeatTransferFromHeatStructure3D1PhaseAC::checkInputs()
{
  const auto flow_channels = flowChannels();

  // The axis that each flow channel is aligned with. The orientation defaults to the x axis, the
  // direction a flow channel is generated along.
  std::vector<unsigned int> axes;
  for (const auto i : index_range(flow_channels))
  {
    const auto * fch = flow_channels[i];
    RealVectorValue direction(1, 0, 0);
    if (fch->isParamValid("orientation"))
      direction = fch->getParam<RealVectorValue>("orientation");
    direction = direction.unit();
    const Real x = std::abs(direction * RealVectorValue(1, 0, 0));
    const Real y = std::abs(direction * RealVectorValue(0, 1, 0));
    const Real z = std::abs(direction * RealVectorValue(0, 0, 1));
    if (x == 1 && y == 0 && z == 0)
      axes.push_back(0);
    else if (x == 0 && y == 1 && z == 0)
      axes.push_back(1);
    else if (x == 0 && y == 0 && z == 1)
      axes.push_back(2);
    else
      paramError("flow_channels",
                 "The flow channel '",
                 _flow_channel_names[i],
                 "' must be aligned with the x-, y-, or z- axis.");
  }

  // All flow channels must be aligned with the same axis and have the same number of elements and
  // the same length
  for (const auto i : make_range(std::size_t(1), flow_channels.size()))
  {
    if (axes[i] != axes[0])
      paramError("flow_channels",
                 "Flow channel '",
                 _flow_channel_names[i],
                 "' has a different axis alignment (",
                 axes[i],
                 "). Make sure all flow channels are aligned with the same axis.");
  }
  for (const auto i : make_range(std::size_t(1), flow_channels.size()))
    if (flow_channels[i]->numElems() != flow_channels[0]->numElems())
      paramError("flow_channels",
                 "Flow channel '",
                 _flow_channel_names[i],
                 "' has ",
                 flow_channels[i]->numElems(),
                 " elements which is inconsistent with the rest of the flow channels. Make sure "
                 "all flow channels have the same number of elements.");
  for (const auto i : make_range(std::size_t(1), flow_channels.size()))
    if (flow_channels[i]->length() != flow_channels[0]->length())
      paramError("flow_channels",
                 "Flow channel '",
                 _flow_channel_names[i],
                 "' has length equal to ",
                 flow_channels[i]->length(),
                 " which is inconsistent with the rest of the flow channels. Make sure all flow "
                 "channels have the length.");
  _layered_average_uo_direction = axes[0];
  _num_layers = flow_channels[0]->numElems();

  if (!hs().hasBoundary(_boundary))
    paramError("boundary",
               "The boundary '",
               _boundary,
               "' does not exist on the component '",
               getParam<std::string>("hs"),
               "'.");
}

MeshAlignment1D3D &
HeatTransferFromHeatStructure3D1PhaseAC::meshAlignment()
{
  if (_mesh_alignment)
    return *_mesh_alignment;

  auto & mesh = getProblem().mesh();
  _mesh_alignment = std::make_unique<MeshAlignment1D3D>(mesh);

  // Elements of the flow channels
  const auto flow_channel_subdomain_ids = mesh.getSubdomainIDs(closuresBlocks());
  const std::set<SubdomainID> flow_channel_subdomains(flow_channel_subdomain_ids.begin(),
                                                      flow_channel_subdomain_ids.end());
  std::vector<dof_id_type> flow_channel_elem_ids;
  for (const auto * elem : mesh.getMesh().active_element_ptr_range())
    if (flow_channel_subdomains.count(elem->subdomain_id()))
      flow_channel_elem_ids.push_back(elem->id());

  // Element sides on the heat structure boundary
  const auto hs_boundary_id = mesh.getBoundaryID(_boundary);
  std::vector<std::tuple<dof_id_type, unsigned short int>> hs_boundary_info;
  for (const auto & [elem_id, side, boundary_id] : mesh.buildSideList())
    if (boundary_id == hs_boundary_id)
      hs_boundary_info.emplace_back(elem_id, side);

  _mesh_alignment->initialize(flow_channel_elem_ids, hs_boundary_info);

  // The coupled elements are not adjacent in the mesh, so their coupling must be declared
  for (const auto fc_elem_id : flow_channel_elem_ids)
    if (_mesh_alignment->hasCoupledSecondaryElemIDs(fc_elem_id))
      for (const auto hs_elem_id : _mesh_alignment->getCoupledSecondaryElemIDs(fc_elem_id))
      {
        _elem_map[fc_elem_id].push_back(hs_elem_id);
        _elem_map[hs_elem_id].push_back(fc_elem_id);
      }

  return *_mesh_alignment;
}

void
HeatTransferFromHeatStructure3D1PhaseAC::actOnAdditionalTasks()
{
  // The inputs that refer to other components are checked before the base class uses them
  if (_current_task == "THM:add_heat_transfers")
    checkInputs();

  HeatTransferFrom1PhaseAC::actOnAdditionalTasks();
}

void
HeatTransferFromHeatStructure3D1PhaseAC::addHeatTransferVariables()
{
  // The wall temperature is constant over each flow channel element, as it is the average of the
  // heat structure temperature over the layer next to it. 'T_wall' holds the same quantity under
  // the name that closures that read the wall temperature expect.
  THMVariableCoordinator::findOrCreate(_awh, _action_factory)
      .requestVariable(false, wallTemperatureName(), "MONOMIAL", "CONSTANT", 1.0, closuresBlocks());
}

void
HeatTransferFromHeatStructure3D1PhaseAC::addHeatTransferICs()
{
  // The heat structure has no initial temperature when the simulation is restarted
  if (!hs().isParamValid("initial_T"))
    return;

  const std::string class_name = "FunctionIC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<VariableName>("variable") = wallTemperatureName();
  params.set<std::vector<SubdomainName>>("block") = closuresBlocks();
  params.set<FunctionName>("function") = hs().getParam<FunctionName>("initial_T");
  getProblem().addInitialCondition(class_name, genName(name(), "T_wall_ic"), params);
}

void
HeatTransferFromHeatStructure3D1PhaseAC::addUserObjects()
{
  ExecFlagEnum execute_on(MooseUtils::getDefaultExecFlagEnum());
  execute_on = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};

  // The position of each flow channel selects the column of the heat structure boundary its
  // layers are averaged over
  std::vector<Point> fch_positions;
  for (const auto * fch : flowChannels())
    fch_positions.push_back(fch->translation());

  {
    const std::string class_name = "NearestPointLayeredSideAverage";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::vector<VariableName>>("variable") = {
        hs().getPhysics().getParam<VariableName>("temperature_name")};
    params.set<std::vector<BoundaryName>>("boundary") = {_boundary};
    params.set<ExecFlagEnum>("execute_on") = execute_on;
    params.set<MooseEnum>("direction") = _layered_average_uo_direction;
    params.set<unsigned int>("num_layers") = _num_layers;
    params.set<std::vector<Point>>("points") = fch_positions;
    getProblem().addUserObject(class_name, _T_wall_avg_uo_name, params);
  }
  {
    const std::string class_name = "ADHeatTransferFromHeatStructure3D1PhaseUserObject";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::vector<SubdomainName>>("block") = closuresBlocks();
    params.set<MeshAlignment1D3D *>("_mesh_alignment") = &meshAlignment();
    params.set<std::vector<VariableName>>("P_hf") = {heatedPerimeterName()};
    params.set<MaterialPropertyName>("Hw") = wallHeatTransferCoefficientName();
    params.set<MaterialPropertyName>("T") = FlowModelSinglePhase::TEMPERATURE;
    params.set<ExecFlagEnum>("execute_on") = execute_on;
    getProblem().addUserObject(class_name, _heat_transfer_uo_name, params);
  }
}

void
HeatTransferFromHeatStructure3D1PhaseAC::addHeatTransferAuxKernels()
{
  ExecFlagEnum execute_on(MooseUtils::getDefaultExecFlagEnum());
  execute_on = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};

  const std::string class_name = "SpatialUserObjectAux";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<SubdomainName>>("block") = closuresBlocks();
  params.set<AuxVariableName>("variable") = wallTemperatureName();
  params.set<ExecFlagEnum>("execute_on") = execute_on;
  params.set<UserObjectName>("user_object") = _T_wall_avg_uo_name;
  getProblem().addAuxKernel(class_name, genName(name(), "T_wall_transfer"), params);
}

void
HeatTransferFromHeatStructure3D1PhaseAC::addHeatTransferKernels()
{
  {
    const std::string class_name = "ADOneD3EqnEnergyHeatFluxFromHeatStructure3D";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::vector<SubdomainName>>("block") = closuresBlocks();
    params.set<std::vector<VariableName>>("P_hf") = {heatedPerimeterName()};
    params.set<MaterialPropertyName>("Hw") = wallHeatTransferCoefficientName();
    params.set<UserObjectName>("user_object") = _T_wall_avg_uo_name;
    params.set<MaterialPropertyName>("T") = FlowModelSinglePhase::TEMPERATURE;
    params.set<NonlinearVariableName>("variable") = THM::RHOEA;
    getProblem().addKernel(class_name, genName(name(), "heat_flux_kernel"), params);
  }
  {
    const std::string class_name = "ADConvectionHeatTransfer3DBC";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::vector<BoundaryName>>("boundary") = {_boundary};
    params.set<NonlinearVariableName>("variable") =
        hs().getPhysics().getParam<VariableName>("temperature_name");
    params.set<UserObjectName>("ht_uo") = _heat_transfer_uo_name;
    getProblem().addBoundaryCondition(class_name, genName(name(), "heat_flux_bc"), params);
  }
}

void
HeatTransferFromHeatStructure3D1PhaseAC::addRelationshipManagers(
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
