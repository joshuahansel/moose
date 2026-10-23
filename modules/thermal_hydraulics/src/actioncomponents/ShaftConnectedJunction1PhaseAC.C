//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ShaftConnectedJunction1PhaseAC.h"
#include "ShaftAC.h"
#include "THMNames.h"
#include "THMUtils.h"
#include "THMVariableCoordinator.h"

InputParameters
ShaftConnectedJunction1PhaseAC::validParams()
{
  InputParameters params = VolumeJunction1PhaseAC::validParams();
  params.makeParamRequired<Real>("A_ref");
  params.addRequiredParam<BoundaryName>("inlet", "Inlet of the component");
  params.addRequiredParam<BoundaryName>("outlet", "Outlet of the component");
  // The connections are the inlet and the outlet
  params.suppressParameter<std::vector<BoundaryName>>("connections");
  params.addRequiredParam<Real>("omega_rated", "Rated speed [rad/s]");
  params.addRequiredParam<Real>("speed_cr_fr", "Speed threshold for friction [-]");
  params.addRequiredParam<Real>("tau_fr_const", "Friction constant [N-m]");
  params.addRequiredParam<std::vector<Real>>("tau_fr_coeff", "Friction coefficients [N-m]");
  params.addRequiredParam<Real>("speed_cr_I", "Speed threshold for inertia [-]");
  params.addRequiredParam<Real>("inertia_const", "Inertia constant [kg-m^2]");
  params.addRequiredParam<std::vector<Real>>("inertia_coeff", "Inertia coefficients [kg-m^2]");
  return params;
}

ShaftConnectedJunction1PhaseAC::ShaftConnectedJunction1PhaseAC(const InputParameters & params)
  : ActionComponent(params),
    VolumeJunction1PhaseAC(params),
    _inertia_const(getParam<Real>("inertia_const")),
    _friction_torque_var_name("friction_torque"),
    _moment_of_inertia_var_name(genName(name(), "moment_of_inertia"))
{
  // The inlet is the first connection and the outlet the second
  _connections = {getParam<BoundaryName>("inlet"), getParam<BoundaryName>("outlet")};
  setupConnections();

  for (const auto & param : {"tau_fr_coeff", "inertia_coeff"})
    if (getParam<std::vector<Real>>(param).size() != 4)
      paramError(param,
                 "Number of entries must be equal to 4, but ",
                 getParam<std::vector<Real>>(param).size(),
                 " were provided.");
}

InputParameters
ShaftConnectedJunction1PhaseAC::shaftConnectedUserObjectParameters(
    const std::string & class_name) const
{
  InputParameters params = junctionUserObjectParameters(class_name);
  params.set<Point>("di_out") = outletDirection();
  params.set<std::vector<VariableName>>("omega") = {omegaVariableName()};
  params.set<bool>("apply_velocity_scaling") = _apply_velocity_scaling;
  for (const auto & param :
       {"omega_rated", "speed_cr_fr", "tau_fr_const", "speed_cr_I", "inertia_const"})
    params.set<Real>(param) = getParam<Real>(param);
  for (const auto & param : {"tau_fr_coeff", "inertia_coeff"})
    params.set<std::vector<Real>>(param) = getParam<std::vector<Real>>(param);
  return params;
}

void
ShaftConnectedJunction1PhaseAC::addMeshGenerators()
{
  // A component that no shaft has connected is reported on the first task on which every component
  // exists, so that the shaft has already told this one that it is connected, if it is
  if (!connectedToShaft())
    paramError("inlet", "This component must be connected to a shaft.");

  VolumeJunction1PhaseAC::addMeshGenerators();
}

VariableName
ShaftConnectedJunction1PhaseAC::omegaVariableName() const
{
  return _awh.getAction<ShaftAC>(_shaft_name).omegaVariableName();
}

void
ShaftConnectedJunction1PhaseAC::addConstantJunctionVariable(const VariableName & var)
{
  _coordinator.requestVariable(false, var, "MONOMIAL", "CONSTANT", 1.0, _blocks);
}

void
ShaftConnectedJunction1PhaseAC::addAuxiliaryVariables()
{
  VolumeJunction1PhaseAC::addAuxiliaryVariables();

  for (const auto & [var, value] : shaftConnectedVariables())
    addConstantJunctionVariable(var);
  addConstantJunctionVariable(_friction_torque_var_name);
  addConstantJunctionVariable(_moment_of_inertia_var_name);
}

void
ShaftConnectedJunction1PhaseAC::addInitialConditions()
{
  VolumeJunction1PhaseAC::addInitialConditions();

  std::vector<std::pair<VariableName, Real>> variables = shaftConnectedVariables();
  variables.emplace_back(_friction_torque_var_name, 0.0);
  variables.emplace_back(_moment_of_inertia_var_name, _inertia_const);
  for (const auto & [var, value] : variables)
  {
    const std::string class_name = "ConstantIC";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<VariableName>("variable") = var;
    params.set<std::vector<SubdomainName>>("block") = _blocks;
    params.set<Real>("value") = value;
    getProblem().addInitialCondition(class_name, genName(name(), var, "ic"), params);
  }
}

void
ShaftConnectedJunction1PhaseAC::addAuxiliaryKernels()
{
  VolumeJunction1PhaseAC::addAuxiliaryKernels();

  auto quantities = auxKernelQuantities();
  addQuantityAuxKernels(auxKernelClass(), auxKernelUserObjectParam(), quantities);
}

void
ShaftConnectedJunction1PhaseAC::addQuantityAuxKernels(
    const std::string & class_name,
    const std::string & uo_param,
    const std::vector<std::pair<std::string, VariableName>> & quantities)
{
  for (const auto & [quantity, var] : quantities)
  {
    InputParameters params = _factory.getValidParams(class_name);
    params.set<AuxVariableName>("variable") = var;
    params.set<MooseEnum>("quantity") = quantity;
    params.set<UserObjectName>(uo_param) = shaftConnectedUserObjectName();
    params.set<std::vector<SubdomainName>>("block") = _blocks;
    getProblem().addAuxKernel(class_name, genName(name(), quantity + "_aux"), params);
  }
}
