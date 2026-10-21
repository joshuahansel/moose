//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatSourceFromPowerDensityAC.h"
#include "HeatStructureInterfaceAC.h"
#include "THMUtils.h"

// Registered under the classic Component's name ('HeatSourceFromPowerDensity'); see
// FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatSourceFromPowerDensityAC,
                           "HeatSourceFromPowerDensity",
                           "add_kernel");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HeatSourceFromPowerDensityAC,
                               "HeatSourceFromPowerDensity");

InputParameters
HeatSourceFromPowerDensityAC::validParams()
{
  InputParameters params = HeatSourceBaseAC::validParams();
  params.addRequiredParam<VariableName>("power_density", "Power density variable");
  params.addClassDescription("Heat generation from a power density variable");
  return params;
}

HeatSourceFromPowerDensityAC::HeatSourceFromPowerDensityAC(const InputParameters & params)
  : ActionComponent(params), HeatSourceBaseAC(params)
{
  THM::informExperimental(*this);
  addRequiredTask("add_kernel");
}

void
HeatSourceFromPowerDensityAC::actOnAdditionalTasks()
{
  if (_current_task != "add_kernel")
    return;

  // This is the first task that uses the heat structure, so check the inputs that refer to it here
  checkRegions();

  const std::string class_name = "CoupledForce";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") =
      hs().getPhysics().getParam<VariableName>("temperature_name");
  params.set<std::vector<SubdomainName>>("block") = subdomains();
  params.set<std::vector<VariableName>>("v") = {getParam<VariableName>("power_density")};
  getProblem().addKernel(class_name, genName(name(), "heat_src"), params);
}
