//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ShaftConnectedCompressor1PhaseAC.h"
#include "THMUtils.h"

// Registered under the classic Component's name ('ShaftConnectedCompressor1Phase'); see
// FlowChannel1PhaseAC.C for why. The tasks are the same as for VolumeJunction1PhaseAC.
#define registerShaftConnectedCompressor1PhaseACTask(task)                                         \
  registerMooseActionAliased("ThermalHydraulicsApp",                                               \
                             ShaftConnectedCompressor1PhaseAC,                                     \
                             "ShaftConnectedCompressor1Phase",                                     \
                             task)
registerShaftConnectedCompressor1PhaseACTask("add_mesh_generator");
registerShaftConnectedCompressor1PhaseACTask("add_variable");
registerShaftConnectedCompressor1PhaseACTask("add_aux_variable");
registerShaftConnectedCompressor1PhaseACTask("add_ic");
registerShaftConnectedCompressor1PhaseACTask("add_user_object");
registerShaftConnectedCompressor1PhaseACTask("add_bc");
registerShaftConnectedCompressor1PhaseACTask("add_kernel");
registerShaftConnectedCompressor1PhaseACTask("add_aux_kernel");
registerShaftConnectedCompressor1PhaseACTask("add_material");
registerShaftConnectedCompressor1PhaseACTask("add_postprocessor");
registerActionComponentAliased("ThermalHydraulicsApp",
                               ShaftConnectedCompressor1PhaseAC,
                               "ShaftConnectedCompressor1Phase");

InputParameters
ShaftConnectedCompressor1PhaseAC::validParams()
{
  InputParameters params = ShaftConnectedJunction1PhaseAC::validParams();
  params.addParam<bool>("treat_as_turbine", false, "Treat the compressor as a turbine?");
  params.addRequiredParam<Real>("mdot_rated", "Rated compressor mass flow rate [kg/s]");
  params.addRequiredParam<Real>("rho0_rated", "Rated compressor stagnation fluid density [kg/m^3]");
  params.addRequiredParam<Real>("c0_rated", "Rated compressor stagnation sound speed [m/s]");
  params.addRequiredParam<std::vector<Real>>(
      "speeds",
      "Relative corrected speeds. Order of speeds needs correspond to the "
      "orders of `Rp_functions` and `eff_functions` [-]");
  params.addRequiredParam<std::vector<FunctionName>>(
      "Rp_functions",
      "Functions of pressure ratio versus relative corrected flow. Each function is for a "
      "different, constant relative corrected speed. The order of function names should correspond "
      "to the order of speeds in the `speeds` parameter [-]");
  params.addRequiredParam<std::vector<FunctionName>>(
      "eff_functions",
      "Functions of adiabatic efficiency versus relative corrected flow. Each function is for a "
      "different, constant relative corrected speed. The order of function names should correspond "
      "to the order of speeds in the `speeds` parameter [-]");
  params.addParam<Real>("min_pressure_ratio", 0.0, "Minimum pressure ratio");
  params.addParam<Real>("max_pressure_ratio", 50.0, "Maximum pressure ratio");
  params.addClassDescription(
      "Compressor between two 1-phase flow channels whose speed is given by a connected shaft");
  return params;
}

ShaftConnectedCompressor1PhaseAC::ShaftConnectedCompressor1PhaseAC(const InputParameters & params)
  : ActionComponent(params), ShaftConnectedJunction1PhaseAC(params)
{
  addRequiredTask("add_postprocessor");
}

void
ShaftConnectedCompressor1PhaseAC::buildVolumeJunctionUserObject()
{
  const std::string class_name = "ADShaftConnectedCompressor1PhaseUserObject";
  InputParameters params = shaftConnectedUserObjectParameters(class_name);
  params.set<bool>("treat_as_turbine") = getParam<bool>("treat_as_turbine");
  for (const auto & param :
       {"mdot_rated", "rho0_rated", "c0_rated", "min_pressure_ratio", "max_pressure_ratio"})
    params.set<Real>(param) = getParam<Real>(param);
  params.set<std::vector<Real>>("speeds") = getParam<std::vector<Real>>("speeds");
  for (const auto & param : {"Rp_functions", "eff_functions"})
    params.set<std::vector<FunctionName>>(param) = getParam<std::vector<FunctionName>>(param);
  params.set<std::string>("compressor_name") = name();
  getProblem().addUserObject(class_name, _junction_uo_name, params);
  THM::connectComponentControllableParam(*this, "K", class_name, _junction_uo_name);
}

std::vector<std::pair<VariableName, Real>>
ShaftConnectedCompressor1PhaseAC::shaftConnectedVariables() const
{
  return {{"delta_p", 0.0}, {"isentropic_torque", 0.0}, {"dissipation_torque", 0.0}};
}

std::vector<std::pair<std::string, VariableName>>
ShaftConnectedCompressor1PhaseAC::auxKernelQuantities() const
{
  return {{"delta_p", "delta_p"},
          {"isentropic_torque", "isentropic_torque"},
          {"dissipation_torque", "dissipation_torque"},
          {"friction_torque", _friction_torque_var_name},
          {"moment_of_inertia", _moment_of_inertia_var_name}};
}

void
ShaftConnectedCompressor1PhaseAC::actOnAdditionalTasks()
{
  ShaftConnectedJunction1PhaseAC::actOnAdditionalTasks();

  if (_current_task == "add_postprocessor")
    for (const auto & quantity :
         {"pressure_ratio", "efficiency", "rel_corrected_flow", "rel_corrected_speed"})
    {
      const std::string class_name = "ShaftConnectedCompressor1PhasePostprocessor";
      InputParameters params = _factory.getValidParams(class_name);
      params.set<MooseEnum>("quantity") = quantity;
      params.set<UserObjectName>("compressor_uo") = _junction_uo_name;
      params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
      getProblem().addPostprocessor(class_name, genName(name(), quantity), params);
    }
}
