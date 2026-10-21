//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TotalPowerAC.h"
#include "THMUtils.h"
#include "InputParameterWarehouse.h"
#include "MooseObjectParameterName.h"

// Registered under the classic Component's name ('TotalPower') rather than this class's own name;
// see FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp", TotalPowerAC, "TotalPower", "add_aux_variable");
registerMooseActionAliased("ThermalHydraulicsApp", TotalPowerAC, "TotalPower", "add_ic");
registerMooseActionAliased("ThermalHydraulicsApp",
                           TotalPowerAC,
                           "TotalPower",
                           "add_aux_scalar_kernel");
registerActionComponentAliased("ThermalHydraulicsApp", TotalPowerAC, "TotalPower");

InputParameters
TotalPowerAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addRequiredParam<Real>("power", "Total power [W]");
  params.declareControllable("power");
  params.addClassDescription("Prescribes total power via a user supplied value");
  return params;
}

TotalPowerAC::TotalPowerAC(const InputParameters & params)
  : ActionComponent(params),
    _power(getParam<Real>("power")),
    _power_var_name(genName(name(), "power"))
{
  THM::informExperimental(*this);

  // Postprocessors that read the power are built before the aux kernel that uses it, so the classic
  // 'component/<name>/power' address is pointed at this component's own parameter now
  getMooseApp().getInputParameterWarehouse().addControllableParameterAlias(
      MooseObjectParameterName("component", name(), "power", "::"),
      MooseObjectParameterName(uniqueName(), "power"));

  addRequiredTask("add_aux_variable");
  addRequiredTask("add_ic");
  addRequiredTask("add_aux_scalar_kernel");
}

void
TotalPowerAC::actOnAdditionalTasks()
{
  if (_current_task == "add_aux_variable")
  {
    const std::string class_name = "MooseVariableScalar";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<MooseEnum>("family") = "SCALAR";
    params.set<MooseEnum>("order") = "FIRST";
    getProblem().addAuxVariable(class_name, _power_var_name, params);
  }
  else if (_current_task == "add_ic")
  {
    const std::string class_name = "ScalarConstantIC";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<VariableName>("variable") = _power_var_name;
    params.set<Real>("value") = _power;
    getProblem().addInitialCondition(class_name, genName(_power_var_name, "ic"), params);
  }
  else if (_current_task == "add_aux_scalar_kernel")
  {
    const std::string class_name = "ConstantScalarAux";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<AuxVariableName>("variable") = _power_var_name;
    params.set<Real>("value") = _power;
    const std::string kernel_name = genName(name(), "power_aux");
    getProblem().addAuxScalarKernel(class_name, kernel_name, params);
    // The classic alias was registered in the constructor (see there); connect the component's
    // own parameter to the object that uses it
    connectControllableParams("power", class_name, kernel_name, "value");
  }
}
