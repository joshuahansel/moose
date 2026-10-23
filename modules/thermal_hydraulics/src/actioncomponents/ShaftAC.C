//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ShaftAC.h"
#include "ShaftConnectableAC.h"
#include "THMUtils.h"
#include "MooseVariableScalar.h"

// Registered under the classic Component's name ('Shaft'); see FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp", ShaftAC, "Shaft", "add_variable");
registerMooseActionAliased("ThermalHydraulicsApp", ShaftAC, "Shaft", "add_ic");
registerMooseActionAliased("ThermalHydraulicsApp", ShaftAC, "Shaft", "add_scalar_kernel");
registerActionComponentAliased("ThermalHydraulicsApp", ShaftAC, "Shaft");

InputParameters
ShaftAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addParam<Real>("scaling_factor_omega", 1.0, "Scaling factor for omega [-]");
  params.addParam<Real>("initial_speed", "Initial shaft speed");
  params.addParam<bool>("ad", true, "Use AD version or not");
  params.addRequiredParam<std::vector<std::string>>("connected_components",
                                                    "Names of the connected components");
  params.addClassDescription("Component that connects torque of turbomachinery components");
  return params;
}

ShaftAC::ShaftAC(const InputParameters & params)
  : ActionComponent(params),
    _scaling_factor_omega(getParam<Real>("scaling_factor_omega")),
    _omega_var_name(genName(name(), "omega")),
    _connected_components(getParam<std::vector<std::string>>("connected_components"))
{
  THM::informExperimental(*this);
  addRequiredTask("add_variable");
  addRequiredTask("add_ic");
  addRequiredTask("add_scalar_kernel");

  // The connected components need to know their shaft before they build their own objects, so this
  // is done as soon as all components exist
  connectComponents();
}

void
ShaftAC::connectComponents()
{
  if (_connected_components.empty())
    return;

  for (const auto & comp_name : _connected_components)
  {
    const ActionComponent * comp = nullptr;
    for (const auto & action : _awh.getActionListByName("list_component"))
      if (action->name() == comp_name)
        comp = dynamic_cast<const ActionComponent *>(action);
    if (!comp)
      continue;
    if (const auto * connectable = dynamic_cast<const ShaftConnectableAC *>(comp))
      connectable->setShaftName(name());
  }
}

std::vector<UserObjectName>
ShaftAC::connectedUserObjectNames() const
{
  std::vector<UserObjectName> uo_names;
  for (const auto & comp_name : _connected_components)
  {
    bool found = false;
    for (const auto & action : _awh.getActionListByName("list_component"))
      if (action->name() == comp_name)
      {
        const auto * connectable = dynamic_cast<const ShaftConnectableAC *>(action);
        if (!connectable)
          paramError("connected_components",
                     "The component '",
                     comp_name,
                     "' is not a component that can be connected to a shaft.");
        uo_names.push_back(connectable->shaftConnectedUserObjectName());
        found = true;
      }
    if (!found)
      paramError("connected_components", "The component '", comp_name, "' does not exist.");
  }
  return uo_names;
}

void
ShaftAC::addSolverVariables()
{
  if (_connected_components.empty())
    paramError("connected_components", "No components are connected to the shaft.");
  if (!isParamValid("initial_speed") && !_app.isRestarting())
    paramError("initial_speed", "The `initial_speed` parameter is missing.");

  const std::string class_name = "MooseVariableScalar";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<MooseEnum>("family") = "SCALAR";
  params.set<MooseEnum>("order") = "FIRST";
  params.set<std::vector<Real>>("scaling") = {_scaling_factor_omega};
  getProblem().addVariable(class_name, _omega_var_name, params);
}

void
ShaftAC::actOnAdditionalTasks()
{
  if (_current_task == "add_ic")
  {
    if (!isParamValid("initial_speed"))
      return;
    const std::string class_name = "ScalarConstantIC";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<VariableName>("variable") = _omega_var_name;
    params.set<Real>("value") = getParam<Real>("initial_speed");
    getProblem().addInitialCondition(class_name, genName(_omega_var_name, "ic"), params);
  }
  else if (_current_task == "add_scalar_kernel")
  {
    const auto uo_names = connectedUserObjectNames();
    const bool ad = getParam<bool>("ad");
    const std::string prefix = ad ? "AD" : "";

    {
      const std::string class_name = prefix + "ShaftTimeDerivativeScalarKernel";
      InputParameters params = getFactory().getValidParams(class_name);
      params.set<NonlinearVariableName>("variable") = _omega_var_name;
      params.set<std::vector<UserObjectName>>("uo_names") = uo_names;
      getProblem().addScalarKernel(class_name, genName(name(), "td"), params);
    }
    for (const auto i : index_range(uo_names))
    {
      const std::string class_name = prefix + "ShaftComponentTorqueScalarKernel";
      InputParameters params = getFactory().getValidParams(class_name);
      params.set<NonlinearVariableName>("variable") = _omega_var_name;
      params.set<UserObjectName>("shaft_connected_component_uo") = uo_names[i];
      getProblem().addScalarKernel(class_name, genName(name(), i, "shaft_speed"), params);
    }
  }
}
