//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FlowChannel1PhaseGasMixPhysics.h"
#include "THMNames.h"
#include "THMUtils.h"
#include "VaporMixtureFluidProperties.h"
#include "FEProblemBase.h"

registerPhysicsBaseTasks("ThermalHydraulicsApp", FlowChannel1PhaseGasMixPhysics);
registerMooseAction("ThermalHydraulicsApp",
                    FlowChannel1PhaseGasMixPhysics,
                    "add_variables_physics");
registerMooseAction("ThermalHydraulicsApp", FlowChannel1PhaseGasMixPhysics, "add_aux_variable");
registerMooseAction("ThermalHydraulicsApp", FlowChannel1PhaseGasMixPhysics, "add_ics_physics");
registerMooseAction("ThermalHydraulicsApp", FlowChannel1PhaseGasMixPhysics, "add_kernel");
registerMooseAction("ThermalHydraulicsApp", FlowChannel1PhaseGasMixPhysics, "add_dg_kernel");
registerMooseAction("ThermalHydraulicsApp", FlowChannel1PhaseGasMixPhysics, "add_aux_kernel");
registerMooseAction("ThermalHydraulicsApp", FlowChannel1PhaseGasMixPhysics, "add_material");
registerMooseAction("ThermalHydraulicsApp", FlowChannel1PhaseGasMixPhysics, "add_user_object");
registerMooseAction("ThermalHydraulicsApp", FlowChannel1PhaseGasMixPhysics, "check_integrity");

InputParameters
FlowChannel1PhaseGasMixPhysics::validParams()
{
  InputParameters params = FlowChannel1PhasePhysics::validParams();
  params.addClassDescription("Physics for a 1D flow channel ActionComponent with a binary gas "
                             "mixture, without wall friction or wall heat transfer.");

  params.addRequiredParam<FunctionName>("initial_mass_fraction",
                                        "Initial mass fraction of the secondary gas");
  params.addParam<Real>(
      "scaling_factor_xirhoA", 1.0, "Scaling factor for the secondary component mass equation");
  return params;
}

FlowChannel1PhaseGasMixPhysics::FlowChannel1PhaseGasMixPhysics(const InputParameters & params)
  : FlowChannel1PhasePhysics(params)
{
  if (!_passives_times_area_names.empty())
    paramError("passives_names", "Passive transport is not supported for a gas mixture.");
}

std::vector<VariableName>
FlowChannel1PhaseGasMixPhysics::solutionVariableNames() const
{
  return {THM::XIRHOA, THM::RHOA, THM::RHOUA, THM::RHOEA};
}

void
FlowChannel1PhaseGasMixPhysics::checkFluidProperties() const
{
  if (!dynamic_cast<const VaporMixtureFluidProperties *>(&getProblem().getUserObjectBase(_fp_name)))
    paramError("fp",
               "The supplied fluid properties object must be of type "
               "'VaporMixtureFluidProperties'.");
}

void
FlowChannel1PhaseGasMixPhysics::addSolverVariables()
{
  FlowChannel1PhasePhysics::addSolverVariables();
  addFlowVariable(
      true, THM::XIRHOA, "MONOMIAL", "CONSTANT", getParam<Real>("scaling_factor_xirhoA"));
}

void
FlowChannel1PhaseGasMixPhysics::addAuxiliaryVariables()
{
  FlowChannel1PhasePhysics::addAuxiliaryVariables();
  addFlowVariable(false, THM::MASS_FRACTION);
}

void
FlowChannel1PhaseGasMixPhysics::addInitialConditions()
{
  FlowChannel1PhasePhysics::addInitialConditions();

  addFunctionIC(THM::MASS_FRACTION, getParam<FunctionName>("initial_mass_fraction"));

  const std::string class_name = "VariableProductIC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<VariableName>("variable") = THM::XIRHOA;
  assignBlocks(params, _blocks);
  params.set<std::vector<VariableName>>("values") = {THM::MASS_FRACTION, THM::RHOA};
  getProblem().addInitialCondition(class_name, genName(name(), "xirhoA_ic"), params);
}

void
FlowChannel1PhaseGasMixPhysics::addRhoEAInitialCondition()
{
  const std::string class_name = "FlowModelGasMixIC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<VariableName>("variable") = THM::RHOEA;
  assignBlocks(params, _blocks);
  params.set<MooseEnum>("quantity") = "rhoEA";
  params.set<FunctionName>("mass_fraction") = getParam<FunctionName>("initial_mass_fraction");
  params.set<FunctionName>("pressure") = _initial_p_fn;
  params.set<FunctionName>("temperature") = _initial_T_fn;
  params.set<FunctionName>("velocity") = _initial_vel_fn;
  params.set<std::vector<VariableName>>("area") = {THM::AREA};
  params.set<UserObjectName>("fluid_properties") = _fp_name;
  getProblem().addInitialCondition(class_name, genName(name(), "rhoEA_ic"), params);
}

void
FlowChannel1PhaseGasMixPhysics::addDensityInitialCondition()
{
  const std::string class_name = "FlowModelGasMixIC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<VariableName>("variable") = THM::DENSITY;
  assignBlocks(params, _blocks);
  params.set<MooseEnum>("quantity") = "rho";
  params.set<FunctionName>("mass_fraction") = getParam<FunctionName>("initial_mass_fraction");
  params.set<FunctionName>("pressure") = _initial_p_fn;
  params.set<FunctionName>("temperature") = _initial_T_fn;
  params.set<FunctionName>("velocity") = _initial_vel_fn;
  params.set<std::vector<VariableName>>("area") = {THM::AREA};
  params.set<UserObjectName>("fluid_properties") = _fp_name;
  getProblem().addInitialCondition(class_name, genName(name(), "rho_ic"), params);
}

void
FlowChannel1PhaseGasMixPhysics::addFEKernels()
{
  FlowChannel1PhasePhysics::addFEKernels();
  addTimeDerivativeKernelIfTransient(THM::XIRHOA);
}

void
FlowChannel1PhaseGasMixPhysics::addDGKernels()
{
  FlowChannel1PhasePhysics::addDGKernels();

  {
    const std::string class_name = "MassDiffusionSpeciesGasMixDGKernel";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<NonlinearVariableName>("variable") = THM::XIRHOA;
    assignBlocks(params, _blocks);
    params.set<std::vector<VariableName>>("A_linear") = {THM::AREA_LINEAR};
    params.set<MaterialPropertyName>("density") = THM::DENSITY;
    params.set<MaterialPropertyName>("diffusion_coefficient") = THM::MASS_DIFFUSION_COEFFICIENT;
    params.set<MaterialPropertyName>("mass_fraction") = THM::MASS_FRACTION;
    params.set<MaterialPropertyName>("direction") = THM::DIRECTION;
    params.set<bool>("implicit") = THM::implicitTimeIntegrationFlag(_app);
    getProblem().addDGKernel(class_name, genName(name(), "mass_diffusion_species"), params);
  }
  {
    const std::string class_name = "MassDiffusionEnergyGasMixDGKernel";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<NonlinearVariableName>("variable") = THM::RHOEA;
    assignBlocks(params, _blocks);
    params.set<std::vector<VariableName>>("A_linear") = {THM::AREA_LINEAR};
    params.set<MaterialPropertyName>("density") = THM::DENSITY;
    params.set<MaterialPropertyName>("pressure") = THM::PRESSURE;
    params.set<MaterialPropertyName>("temperature") = THM::TEMPERATURE;
    params.set<MaterialPropertyName>("velocity") = THM::VELOCITY;
    params.set<MaterialPropertyName>("diffusion_coefficient") = THM::MASS_DIFFUSION_COEFFICIENT;
    params.set<MaterialPropertyName>("mass_fraction") = THM::MASS_FRACTION;
    params.set<MaterialPropertyName>("direction") = THM::DIRECTION;
    params.set<UserObjectName>("fluid_properties") = _fp_name;
    params.set<bool>("implicit") = THM::implicitTimeIntegrationFlag(_app);
    getProblem().addDGKernel(class_name, genName(name(), "mass_diffusion_energy"), params);
  }
}

void
FlowChannel1PhaseGasMixPhysics::addAdvectionDGKernels()
{
  for (const auto & var : solutionVariableNames())
  {
    const std::string class_name = "NumericalFluxGasMixDGKernel";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<NonlinearVariableName>("variable") = var;
    assignBlocks(params, _blocks);
    params.set<std::vector<VariableName>>("A_linear") = {THM::AREA_LINEAR};
    params.set<std::vector<VariableName>>("xirhoA") = {THM::XIRHOA};
    params.set<std::vector<VariableName>>("rhoA") = {THM::RHOA};
    params.set<std::vector<VariableName>>("rhouA") = {THM::RHOUA};
    params.set<std::vector<VariableName>>("rhoEA") = {THM::RHOEA};
    params.set<UserObjectName>("numerical_flux") = _numerical_flux_name;
    params.set<bool>("implicit") = THM::implicitTimeIntegrationFlag(_app);
    getProblem().addDGKernel(class_name, genName(name(), "advection:" + var), params);
  }
}

void
FlowChannel1PhaseGasMixPhysics::addFluidPropertiesMaterial()
{
  const std::string class_name = "FluidPropertiesGasMixMaterial";
  InputParameters params = getFactory().getValidParams(class_name);
  assignBlocks(params, _blocks);
  params.set<std::vector<VariableName>>("xirhoA") = {THM::XIRHOA};
  params.set<std::vector<VariableName>>("rhoA") = {THM::RHOA};
  params.set<std::vector<VariableName>>("rhouA") = {THM::RHOUA};
  params.set<std::vector<VariableName>>("rhoEA") = {THM::RHOEA};
  params.set<std::vector<VariableName>>("area") = {THM::AREA};
  params.set<UserObjectName>("fluid_properties") = _fp_name;
  getProblem().addMaterial(class_name, genName(name(), "fp_mat"), params);
}

void
FlowChannel1PhaseGasMixPhysics::addSlopeReconstructionMaterial()
{
  const std::string class_name = "SlopeReconstructionGasMixMaterial";
  InputParameters params = getFactory().getValidParams(class_name);
  assignBlocks(params, _blocks);
  params.set<MooseEnum>("scheme") = getParam<MooseEnum>("rdg_slope_reconstruction");
  params.set<std::vector<VariableName>>("A_elem") = {THM::AREA};
  params.set<std::vector<VariableName>>("A_linear") = {THM::AREA_LINEAR};
  params.set<std::vector<VariableName>>("xirhoA") = {THM::XIRHOA};
  params.set<std::vector<VariableName>>("rhoA") = {THM::RHOA};
  params.set<std::vector<VariableName>>("rhouA") = {THM::RHOUA};
  params.set<std::vector<VariableName>>("rhoEA") = {THM::RHOEA};
  params.set<UserObjectName>("fluid_properties") = _fp_name;
  params.set<bool>("implicit") = THM::implicitTimeIntegrationFlag(_app);
  getProblem().addMaterial(class_name, genName(name(), "slope_mat"), params);
}

void
FlowChannel1PhaseGasMixPhysics::addNumericalFluxUserObject()
{
  const std::string class_name = "NumericalFluxGasMixHLLC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<UserObjectName>("fluid_properties") = _fp_name;
  params.set<MooseEnum>("emit_on_nan") = "none";
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};
  addUserObject(class_name, _numerical_flux_name, params);
}

void
FlowChannel1PhaseGasMixPhysics::addPressureTemperatureAuxKernels()
{
  for (const auto & [var, quantity, suffix] :
       std::vector<std::tuple<std::string, std::string, std::string>>{
           {THM::PRESSURE, "p", "p_aux"}, {THM::TEMPERATURE, "T", "T_aux"}})
  {
    const std::string class_name = "FlowModelGasMixAux";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<AuxVariableName>("variable") = var;
    assignBlocks(params, _blocks);
    params.set<MooseEnum>("quantity") = quantity;
    params.set<std::vector<VariableName>>("xirhoA") = {THM::XIRHOA};
    params.set<std::vector<VariableName>>("rhoA") = {THM::RHOA};
    params.set<std::vector<VariableName>>("rhouA") = {THM::RHOUA};
    params.set<std::vector<VariableName>>("rhoEA") = {THM::RHOEA};
    params.set<std::vector<VariableName>>("area") = {THM::AREA};
    params.set<UserObjectName>("fluid_properties") = _fp_name;
    params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
    getProblem().addAuxKernel(class_name, genName(name(), suffix), params);
  }
}

void
FlowChannel1PhaseGasMixPhysics::addAuxiliaryKernels()
{
  FlowChannel1PhasePhysics::addAuxiliaryKernels();

  const std::string class_name = "QuotientAux";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<AuxVariableName>("variable") = THM::MASS_FRACTION;
  assignBlocks(params, _blocks);
  params.set<std::vector<VariableName>>("numerator") = {THM::XIRHOA};
  params.set<std::vector<VariableName>>("denominator") = {THM::RHOA};
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
  getProblem().addAuxKernel(class_name, genName(name(), "xi_aux"), params);
}
