//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FlowChannel1PhaseAC.h"
#include "SlopeReconstruction1DInterface.h"
#include "FlowBoundary1PhaseAC.h"
#include "THMUtils.h"
#include "FlowChannel1PhasePhysics.h"
#include "ActionFactory.h"
#include "ActionWarehouse.h"
#include "FEProblemBase.h"

#include <numeric>

// Registered under the classic Component's name ('FlowChannel1Phase') rather than this class's own
// name, so an input file can move a flow channel from [Components] to [ActionComponents] without
// also renaming 'type'. The C++ class itself keeps its 'AC' suffix to stay distinguishable from
// the classic FlowChannel1Phase Component in the source code.
registerMooseActionAliased("ThermalHydraulicsApp",
                           FlowChannel1PhaseAC,
                           "FlowChannel1Phase",
                           "add_mesh_generator");
// FlowChannel1PhaseAC creates and attaches its own FlowChannel1PhasePhysics; see the constructor
registerMooseActionAliased("ThermalHydraulicsApp",
                           FlowChannel1PhaseAC,
                           "FlowChannel1Phase",
                           "init_component_physics");
registerMooseActionAliased("ThermalHydraulicsApp",
                           FlowChannel1PhaseAC,
                           "FlowChannel1Phase",
                           "add_functor_material");
registerMooseActionAliased("ThermalHydraulicsApp",
                           FlowChannel1PhaseAC,
                           "FlowChannel1Phase",
                           "add_postprocessor");
registerMooseActionAliased("ThermalHydraulicsApp",
                           FlowChannel1PhaseAC,
                           "FlowChannel1Phase",
                           "add_convergence");
registerMooseActionAliased("ThermalHydraulicsApp",
                           FlowChannel1PhaseAC,
                           "FlowChannel1Phase",
                           "add_vector_postprocessor");
registerActionComponentAliased("ThermalHydraulicsApp", FlowChannel1PhaseAC, "FlowChannel1Phase");

InputParameters
FlowChannel1PhaseAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params += ComponentMeshTransformHelper::validParams();
  params += GravityInterface::validParams();
  params.addClassDescription(
      "A 1D single-phase flow channel, with its mesh generated through the MeshGenerator system.");

  // The physics this component creates; derived components with other equations change it
  params.addPrivateParam<std::string>("_physics_type", "FlowChannel1PhasePhysics");

  params.addParam<Real>("p_ref", 101.325e3, "Reference pressure [Pa]");
  params.addParam<Real>("T_ref", 273.15, "Reference temperature [K]");
  params.addParam<Real>("vel_ref", 1.0, "Reference velocity [m/s]");
  params.addParam<Real>("p_rel_step_tol", 1e-5, "Pressure relative step tolerance");
  params.addParam<Real>("T_rel_step_tol", 1e-5, "Temperature relative step tolerance");
  params.addParam<Real>("vel_rel_step_tol", 1e-5, "Velocity relative step tolerance");
  params.addParam<Real>("mass_res_tol", 1e-5, "Mass equation normalized residual tolerance");
  params.addParam<Real>(
      "momentum_res_tol", 1e-5, "Momentum equation normalized residual tolerance");
  params.addParam<Real>("energy_res_tol", 1e-5, "Energy equation normalized residual tolerance");
  params.addParam<bool>("enable_heat_conduction", false, "Enable heat conduction?");
  params.addParam<std::vector<VariableName>>(
      "vpp_vars", {}, "Variables to add in an ElementValueSampler");
  params.addParam<bool>(
      "create_flux_vpp",
      false,
      "If true, create a VectorPostprocessor with the the mass, momentum, and energy side fluxes");
  params.addParam<bool>("name_multiple_ht_by_index",
                        true,
                        "If true, when there are multiple heat transfer components connected to "
                        "this flow channel, use their index for naming related quantities; "
                        "otherwise, use the name of the heat transfer component.");
  params.addRequiredRangeCheckedParam<std::vector<Real>>(
      "length", "length > 0", "Length of each axial section of the flow channel [m]");
  params.addRequiredRangeCheckedParam<std::vector<unsigned int>>(
      "n_elems", "n_elems > 0", "Number of elements of each axial section of the flow channel");
  params.addParam<std::vector<std::string>>(
      "axial_region_names",
      {},
      "Names of the axial sections. This is required if there is more than 1 axial section.");

  // These are forwarded as-is to the FlowChannel1PhasePhysics this component creates (see the
  // constructor) - kept here, rather than requiring a separate [Physics] block, so that an input
  // file can move a flow channel from [Components]/FlowChannel1Phase to
  // [ActionComponents]/FlowChannel1PhaseAC with no other changes.
  params.addRequiredParam<UserObjectName>("fp",
                                          "The name of the user object that defines fluid "
                                          "properties");
  params.addRequiredParam<FunctionName>(
      "A", "Area of the flow channel, can be a constant or a function");
  params.addRequiredParam<FunctionName>("initial_p", "Initial pressure [Pa]");
  params.addRequiredParam<FunctionName>("initial_T", "Initial temperature [K]");
  params.addRequiredParam<FunctionName>("initial_vel", "Initial velocity [m/s]");
  params.addParam<std::vector<Real>>(
      "scaling_factor_1phase",
      std::vector<Real>(3, 1.0),
      "Scaling factors for each single phase variable (rhoA, rhouA, rhoEA)");
  params.addParam<std::vector<FunctionName>>(
      "initial_passives",
      {},
      "Initial passive transport variable values in the flow channel, if any (units are "
      "[amount/m^3], where 'amount' may be mass (kg) or a number (molecules, moles, etc.))");
  params.addParam<std::vector<VariableName>>(
      "passives_names",
      {},
      "Names for each passive transport variable [amount/m^3]. Note that the conserved (solution) "
      "variables will be an amount per unit volume multiplied by the channel cross-sectional area, "
      "yielding an amount per unit length; these solution variable names will append '_times_area' "
      "to the names given in this parameter.");
  params.addParam<std::vector<Real>>("scaling_factor_passives",
                                     "Scaling factor for each passive transport variable");
  params.addParam<std::vector<std::string>>(
      "closures",
      {},
      "Closures object(s) providing wall friction. This is optional since closure relations can "
      "be supplied directly by Materials as well.");
  params.addParam<FunctionName>("f", "Wall friction factor [-]");
  params.addParam<FunctionName>("D_h", "Hydraulic diameter [m]");
  params.addParam<Real>("roughness", 0.0, "Roughness [m]");
  params.addParam<Real>("PoD", 1, "Pitch-to-diameter ratio for parallel bundle heat transfer [-]");
  params.addParam<MooseEnum>(
      "rdg_slope_reconstruction",
      SlopeReconstruction1DInterface<true>::getSlopeReconstructionMooseEnum("None"),
      "Slope reconstruction type for rDG spatial discretization");

  // Use THM's classic naming for the flow channel's orientation, instead of
  // ComponentMeshTransformHelper's generic 'direction'. InputParameters::renameParam() redirects
  // ComponentMeshTransformHelper's internal queryParam("direction") calls to the 'orientation'
  // value, so no change is needed in the mixin itself.
  params.renameParam("direction",
                     "orientation",
                     "Direction to orient the flow channel with, assuming it is initially "
                     "oriented along the X-axis (1, 0, 0).");

  return params;
}

FlowChannel1PhaseAC::FlowChannel1PhaseAC(const InputParameters & params)
  : ActionComponent(params),
    ComponentMeshTransformHelper(params),
    GravityInterface(params),
    _lengths(getParam<std::vector<Real>>("length")),
    _n_elems_sections(getParam<std::vector<unsigned int>>("n_elems")),
    _axial_region_names(getParam<std::vector<std::string>>("axial_region_names")),
    _length(std::accumulate(_lengths.begin(), _lengths.end(), 0.0)),
    _n_elems(std::accumulate(_n_elems_sections.begin(), _n_elems_sections.end(), 0u)),
    _physics(nullptr)
{
  if (_lengths.size() != _n_elems_sections.size())
    paramError("n_elems",
               "The parameters 'length' and 'n_elems' must have the same number of entries.");
  if (_lengths.size() > 1 && _axial_region_names.size() != _lengths.size())
    paramError(
        "axial_region_names",
        "If there is more than 1 axial section, then the parameter 'axial_region_names' must "
        "be provided, with as many entries as 'length'.");
  if (!_axial_region_names.empty() && _axial_region_names.size() != _lengths.size())
    paramError("axial_region_names", "The number of entries must be equal to that of 'length'.");

  THM::informExperimental(*this);
  _dimension = 1;
  // ComponentMeshTransformHelper adds its own required task
  addRequiredTask("add_mesh_generator");
  addRequiredTask("init_component_physics");
  addRequiredTask("add_functor_material");
  addRequiredTask("add_postprocessor");
  addRequiredTask("add_convergence");
  if (!getParam<std::vector<VariableName>>("vpp_vars").empty() || getParam<bool>("create_flux_vpp"))
    addRequiredTask("add_vector_postprocessor");

  // Build and register the physics defining this flow channel's governing equations, mirroring
  // how AddActionComponentAction itself builds this component - so the user does not need to
  // declare a separate [Physics] block (see class documentation).
  const std::string physics_type = getParam<std::string>("_physics_type");
  InputParameters physics_params = _action_factory.getValidParams(physics_type);
  physics_params.blockFullpath() = params.blockFullpath();
  physics_params.set<bool>("_built_by_moose") = true;
  physics_params.set<std::string>("registered_identifier") = "(AutoBuilt)";

  physics_params.set<UserObjectName>("fp") = getParam<UserObjectName>("fp");
  physics_params.set<FunctionName>("A") = getParam<FunctionName>("A");
  physics_params.set<FunctionName>("initial_p") = getParam<FunctionName>("initial_p");
  physics_params.set<FunctionName>("initial_T") = getParam<FunctionName>("initial_T");
  physics_params.set<FunctionName>("initial_vel") = getParam<FunctionName>("initial_vel");
  // A flow channel with other equations gives the scaling factors of its variables one by one
  if (params.have_parameter<Real>("scaling_factor_rhoA"))
    physics_params.set<std::vector<Real>>("scaling_factor_1phase") = {
        getParam<Real>("scaling_factor_rhoA"),
        getParam<Real>("scaling_factor_rhouA"),
        getParam<Real>("scaling_factor_rhoEA")};
  else
    physics_params.set<std::vector<Real>>("scaling_factor_1phase") =
        getParam<std::vector<Real>>("scaling_factor_1phase");
  if (params.have_parameter<Real>("scaling_factor_xirhoA"))
    physics_params.set<Real>("scaling_factor_xirhoA") = getParam<Real>("scaling_factor_xirhoA");
  if (params.have_parameter<FunctionName>("initial_mass_fraction"))
    physics_params.set<FunctionName>("initial_mass_fraction") =
        getParam<FunctionName>("initial_mass_fraction");
  physics_params.set<RealVectorValue>("gravity_vector") = _gravity_vector;
  physics_params.set<std::vector<FunctionName>>("initial_passives") =
      getParam<std::vector<FunctionName>>("initial_passives");
  physics_params.set<std::vector<VariableName>>("passives_names") =
      getParam<std::vector<VariableName>>("passives_names");
  if (isParamValid("scaling_factor_passives"))
    physics_params.set<std::vector<Real>>("scaling_factor_passives") =
        getParam<std::vector<Real>>("scaling_factor_passives");
  physics_params.set<std::vector<std::string>>("closures") =
      getParam<std::vector<std::string>>("closures");
  if (isParamValid("f"))
    physics_params.set<FunctionName>("f") = getParam<FunctionName>("f");
  if (isParamValid("D_h"))
    physics_params.set<FunctionName>("D_h") = getParam<FunctionName>("D_h");
  physics_params.set<Real>("roughness") = getParam<Real>("roughness");
  physics_params.set<Real>("PoD") = getParam<Real>("PoD");
  physics_params.set<bool>("enable_heat_conduction") = getParam<bool>("enable_heat_conduction");
  physics_params.set<bool>("name_multiple_ht_by_index") =
      getParam<bool>("name_multiple_ht_by_index");
  physics_params.set<MooseEnum>("rdg_slope_reconstruction") =
      getParam<MooseEnum>("rdg_slope_reconstruction");

  auto physics_action = _action_factory.create(physics_type, name() + "_physics", physics_params);
  _physics = dynamic_cast<FlowChannel1PhasePhysics *>(physics_action.get());
  if (!_physics)
    mooseError("Internal error: failed to create this component's FlowChannel1PhasePhysics");
  _awh.addActionBlock(physics_action);
}

void
FlowChannel1PhaseAC::addMeshGenerators()
{
  InputParameters params = _factory.getValidParams("GeneratedMeshGenerator");
  params.set<MooseEnum>("dim") = "1";
  params.set<Real>("xmax") = _length;
  params.set<unsigned int>("nx") = _n_elems;
  params.set<std::string>("boundary_name_prefix") = name();
  params.set<SubdomainName>("subdomain_name") = name();
  _app.getMeshGeneratorSystem().addMeshGenerator(
      "GeneratedMeshGenerator", name() + "_base", params);
  _mg_names.push_back(name() + "_base");
  _blocks.push_back(name());

  // Axial sections each have their own block, named like classic THM's '<name>:<section>'. The
  // elements are uniform in each section, so the sections are split out of the uniform mesh above
  // by their element ranges: the size of an element differs between sections, so the mesh is built
  // by moving the nodes of the uniform mesh to the section boundaries.
  if (_lengths.size() > 1)
  {
    // The uniform mesh has nx elements of the same size; map the node at element index k to the
    // position along the flow channel that the sections give
    // Built from the last section backwards as nested if(x < end_of_section, section, ...)
    const Real dx_uniform = _length / _n_elems;
    std::vector<std::string> section_functions;
    std::vector<Real> section_ends;
    {
      Real x_start = 0;
      unsigned int k_start = 0;
      for (const auto s : index_range(_lengths))
      {
        const Real dx = _lengths[s] / _n_elems_sections[s];
        // Within the section, x_new = x_start + (x - k_start * dx_uniform) / dx_uniform * dx
        section_functions.push_back(Moose::stringifyExact(x_start) + " + (x - " +
                                    Moose::stringifyExact(k_start * dx_uniform) + ") / " +
                                    Moose::stringifyExact(dx_uniform) + " * " +
                                    Moose::stringifyExact(dx));
        k_start += _n_elems_sections[s];
        section_ends.push_back(k_start * dx_uniform);
        x_start += _lengths[s];
      }
    }
    std::string x_function = section_functions.back();
    for (int s = int(_lengths.size()) - 2; s >= 0; --s)
      x_function = "if(x < " + Moose::stringifyExact(section_ends[s]) + ", " +
                   section_functions[s] + ", " + x_function + ")";

    InputParameters nodes_params = _factory.getValidParams("ParsedNodeTransformGenerator");
    nodes_params.set<MeshGeneratorName>("input") = _mg_names.back();
    nodes_params.set<ParsedFunctionExpression>("x_function") = x_function;
    _app.getMeshGeneratorSystem().addMeshGenerator(
        "ParsedNodeTransformGenerator", name() + "_sections", nodes_params);
    _mg_names.push_back(name() + "_sections");

    // Subdomains of the sections
    _blocks.clear();
    Real x_lo = 0;
    for (const auto s : index_range(_lengths))
    {
      const std::string step_name = name() + "_" + _axial_region_names[s] + "_block";
      InputParameters block_params = _factory.getValidParams("ParsedSubdomainMeshGenerator");
      block_params.set<MeshGeneratorName>("input") = _mg_names.back();
      block_params.set<std::string>("combinatorial_geometry") =
          "x > " + Moose::stringifyExact(x_lo) + " & x < " +
          Moose::stringifyExact(x_lo + _lengths[s]);
      // Block ids are made unique by the combiner of the components, so the ones here only need to
      // be different between the sections
      block_params.set<subdomain_id_type>("block_id") = s + 1;
      block_params.set<SubdomainName>("block_name") = name() + ":" + _axial_region_names[s];
      block_params.set<std::vector<SubdomainName>>("excluded_subdomains") = {};
      _app.getMeshGeneratorSystem().addMeshGenerator(
          "ParsedSubdomainMeshGenerator", step_name, block_params);
      _mg_names.push_back(step_name);
      _blocks.push_back(name() + ":" + _axial_region_names[s]);
      x_lo += _lengths[s];
    }
  }

  // Rename the generated boundaries to match THM's <name>:in / <name>:out convention, used by
  // classic THM Components (Component1D) for inlet/outlet/junction boundary references.
  InputParameters rename_params = _factory.getValidParams("RenameBoundaryGenerator");
  rename_params.set<MeshGeneratorName>("input") = _mg_names.back();
  rename_params.set<std::vector<BoundaryName>>("old_boundary") = {name() + "_left",
                                                                  name() + "_right"};
  rename_params.set<std::vector<BoundaryName>>("new_boundary") = {name() + ":in", name() + ":out"};
  _app.getMeshGeneratorSystem().addMeshGenerator(
      "RenameBoundaryGenerator", name() + "_renamed", rename_params);
  _mg_names.push_back(name() + "_renamed");

  addShapeMeshGenerators();

  // Classic THM boundary components create a sideset named after themselves at the end of the flow
  // channel they connect to, which inputs use to refer to the boundary (e.g. in postprocessors that
  // look up the boundary's UserObject by this name). Create the same sidesets for the boundary
  // components attached to this flow channel.
  for (const auto * boundary : _awh.getActions<FlowBoundary1PhaseAC>())
    if (boundary->connectedComponentName() == name())
    {
      const std::string step_name = name() + "_" + boundary->name() + "_sideset";
      InputParameters sideset_params = _factory.getValidParams("ParsedGenerateSideset");
      sideset_params.set<MeshGeneratorName>("input") = _mg_names.back();
      sideset_params.set<std::vector<BoundaryName>>("included_boundaries") = {boundary->input()};
      sideset_params.set<std::string>("combinatorial_geometry") = "1";
      sideset_params.set<BoundaryName>("new_sideset_name") = boundary->name();
      _app.getMeshGeneratorSystem().addMeshGenerator(
          "ParsedGenerateSideset", step_name, sideset_params);
      _mg_names.push_back(step_name);
    }

  // Add a nodeset with all the nodes of the flow channel, named after the component as in classic
  // THM's Component1D. Heat transfer components use it to transfer values onto the flow channel.
  InputParameters nodeset_params = _factory.getValidParams("ParsedGenerateNodeset");
  nodeset_params.set<MeshGeneratorName>("input") = _mg_names.back();
  nodeset_params.set<BoundaryName>("new_nodeset_name") = name();
  nodeset_params.set<std::vector<SubdomainName>>("included_subdomains") = _blocks;
  nodeset_params.set<std::string>("expression") = "1";
  // The expression does not depend on the distribution of the mesh
  nodeset_params.set<bool>("allow_distributed_meshes") = true;
  _app.getMeshGeneratorSystem().addMeshGenerator(
      "ParsedGenerateNodeset", name() + "_nodeset", nodeset_params);
  _mg_names.push_back(name() + "_nodeset");

  _top_mg_name = _mg_names.back();

  ComponentMeshTransformHelper::addMeshGenerators();
}

void
FlowChannel1PhaseAC::addPhysics()
{
  _physics->addComponent(*this);
}

const FlowChannel1PhaseAC &
FlowChannel1PhaseAC::connectedFlowChannel(const ActionWarehouse & awh,
                                          const std::string & name,
                                          const std::string & from)
{
  for (const auto & action : awh.getActionListByName("list_component"))
    if (action->name() == name)
    {
      if (const auto * fch = dynamic_cast<const FlowChannel1PhaseAC *>(action))
        return *fch;
      ::mooseError(
          from, ": Trying to connect to the component '", name, "', which is not a flow channel.");
    }
  ::mooseError(from,
               ": Trying to connect to a component '",
               name,
               "', but there is no such component in the simulation. Please check your spelling.");
}

const FlowChannel1PhasePhysics &
FlowChannel1PhaseAC::getPhysics() const
{
  mooseAssert(_physics, "The physics has not been created yet");
  return *_physics;
}

std::optional<ConvergenceName>
FlowChannel1PhaseAC::nonlinearConvergenceName() const
{
  return ConvergenceName(name() + ":nlconv");
}

void
FlowChannel1PhaseAC::actOnAdditionalTasks()
{
  if (_current_task == "add_functor_material" || _current_task == "add_postprocessor" ||
      _current_task == "add_convergence")
    addConvergence();
  else if (_current_task == "add_vector_postprocessor")
    addVectorPostprocessors();
}

std::string
FlowChannel1PhaseAC::sortBy() const
{
  // The coordinate of the dominant direction of the flow channel, along which values are sorted
  const RealVectorValue dir = isParamValid("orientation") ? direction() : RealVectorValue(1, 0, 0);
  const std::vector<Real> abs_dir = {std::abs(dir(0)), std::abs(dir(1)), std::abs(dir(2))};
  const auto dominant = std::max_element(abs_dir.begin(), abs_dir.end()) - abs_dir.begin();
  return std::string(1, "xyz"[dominant]);
}

void
FlowChannel1PhaseAC::addVectorPostprocessors()
{
  const auto & vpp_vars = getParam<std::vector<VariableName>>("vpp_vars");
  if (!vpp_vars.empty())
  {
    const std::string class_name = "ElementValueSampler";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<std::vector<SubdomainName>>("block") = _blocks;
    params.set<std::vector<VariableName>>("variable") = vpp_vars;
    params.set<std::string>("sort_by") = sortBy();
    params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
    getProblem().addVectorPostprocessor(class_name, name() + "_vars_vpp", params);
  }
  if (getParam<bool>("create_flux_vpp"))
  {
    const std::string class_name = "NumericalFlux3EqnInternalValues";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<std::vector<SubdomainName>>("block") = _blocks;
    params.set<UserObjectName>("numerical_flux") = getPhysics().numericalFluxUserObjectName();
    params.set<std::vector<VariableName>>("A_linear") = {THM::AREA_LINEAR};
    params.set<std::string>("sort_by") = sortBy();
    params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
    getProblem().addVectorPostprocessor(class_name, name() + "_flux_vpp", params);
  }
}

void
FlowChannel1PhaseAC::addConvergence()
{
  const std::vector<std::pair<std::string, Real>> var_norm_pairs{
      {THM::PRESSURE, getParam<Real>("p_ref")},
      {THM::TEMPERATURE, getParam<Real>("T_ref")},
      {THM::VELOCITY, getParam<Real>("vel_ref")}};
  const std::vector<std::pair<std::string, std::string>> var_eq_pairs{
      {THM::RHOA, "mass"}, {THM::RHOUA, "momentum"}, {THM::RHOEA, "energy"}};
  const auto pp_name = [this](const std::string & base) { return name() + ":" + base; };

  if (_current_task == "add_functor_material")
  {
    // The functor of the flow channel's state
    {
      const std::string class_name = "FlowModel1PhaseFunctorMaterial";
      InputParameters params = _factory.getValidParams(class_name);
      params.set<std::vector<SubdomainName>>("block") = _blocks;
      params.set<UserObjectName>("fluid_properties") = getPhysics().fluidPropertiesName();
      getProblem().addFunctorMaterial(class_name, pp_name("fm1phase_fmat"), params);
    }
    for (const auto & [var, norm] : var_norm_pairs)
      THM::addNonlinearStepFunctorMaterial(getProblem(),
                                           _factory,
                                           pp_name(var + "_change_fmat"),
                                           _blocks,
                                           THM::functorMaterialPropertyName<false>(var),
                                           var + "_change",
                                           false);
  }
  else if (_current_task == "add_postprocessor")
  {
    for (const auto & [var, norm] : var_norm_pairs)
      THM::addMaximumFunctorPostprocessor(
          getProblem(), _factory, var + "_change", pp_name(var + "_rel_step"), norm, _blocks);

    // The flow channel is straight and uniform, so the middle point is that of the line segment
    // and the minimum element size is the element size
    const Point start = translation();
    const Point end = start + _length * (isParamValid("orientation") ? direction().unit()
                                                                     : RealVectorValue(1, 0, 0));
    for (const auto & [var, eq] : var_eq_pairs)
    {
      const std::string class_name = "Normalized1PhaseResidualNorm";
      InputParameters params = _factory.getValidParams(class_name);
      params.set<VariableName>("variable") = var;
      params.set<std::vector<SubdomainName>>("block") = _blocks;
      params.set<MooseEnum>("norm_type") = "l_inf";
      params.set<Real>("p_ref") = getParam<Real>("p_ref");
      params.set<Real>("T_ref") = getParam<Real>("T_ref");
      params.set<Real>("vel_ref") = getParam<Real>("vel_ref");
      params.set<FunctionName>("A") = getParam<FunctionName>("A");
      params.set<Point>("point") = 0.5 * (start + end);
      params.set<UserObjectName>("fluid_properties") = getPhysics().fluidPropertiesName();
      params.set<Real>("min_elem_size") = _length / _n_elems;
      params.set<ExecFlagEnum>("execute_on") = EXEC_NONLINEAR_CONVERGENCE;
      params.set<std::vector<OutputName>>("outputs") = {"none"};
      getProblem().addPostprocessor(class_name, pp_name(eq + "_res"), params);
    }
  }
  else if (_current_task == "add_convergence")
    THM::addMultiPostprocessorConvergence(
        getProblem(),
        _factory,
        *nonlinearConvergenceName(),
        {pp_name(THM::PRESSURE + "_rel_step"),
         pp_name(THM::TEMPERATURE + "_rel_step"),
         pp_name(THM::VELOCITY + "_rel_step"),
         pp_name("mass_res"),
         pp_name("momentum_res"),
         pp_name("energy_res")},
        {"step: p", "step: T", "step: vel", "res: mass", "res: momentum", "res: energy"},
        {getParam<Real>("p_rel_step_tol"),
         getParam<Real>("T_rel_step_tol"),
         getParam<Real>("vel_rel_step_tol"),
         getParam<Real>("mass_res_tol"),
         getParam<Real>("momentum_res_tol"),
         getParam<Real>("energy_res_tol")});
}
