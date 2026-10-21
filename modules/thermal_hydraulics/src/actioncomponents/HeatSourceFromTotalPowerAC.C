//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatSourceFromTotalPowerAC.h"
#include "HeatStructureInterfaceAC.h"
#include "HeatStructurePlateAC.h"
#include "TotalPowerAC.h"
#include "THMUtils.h"
#include "THMNames.h"

// Registered under the classic Component's name ('HeatSourceFromTotalPower'); see
// FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatSourceFromTotalPowerAC,
                           "HeatSourceFromTotalPower",
                           "add_function");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatSourceFromTotalPowerAC,
                           "HeatSourceFromTotalPower",
                           "add_postprocessor");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatSourceFromTotalPowerAC,
                           "HeatSourceFromTotalPower",
                           "add_kernel");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HeatSourceFromTotalPowerAC,
                               "HeatSourceFromTotalPower");

InputParameters
HeatSourceFromTotalPowerAC::validParams()
{
  InputParameters params = HeatSourceBaseAC::validParams();
  params.addRequiredParam<std::string>("power", "Component that provides total power");
  params.addParam<Real>(
      "power_fraction", 1., "Fraction of the total power that goes into the heat structure [-]");
  params.addParam<FunctionName>("power_shape_function", "Axial power shape [-]");
  params.declareControllable("power_fraction");
  params.addClassDescription("Heat generation from total power");
  return params;
}

HeatSourceFromTotalPowerAC::HeatSourceFromTotalPowerAC(const InputParameters & params)
  : ActionComponent(params),
    HeatSourceBaseAC(params),
    _has_psf(isParamValid("power_shape_function")),
    _power_shape_fn(_has_psf ? getParam<FunctionName>("power_shape_function")
                             : FunctionName(genName(name(), "power_shape_fn"))),
    _power_shape_integral_name(_has_psf ? genName(name(), _power_shape_fn, "integral")
                                        : genName(_power_shape_fn, "integral"))
{
  THM::informExperimental(*this);
  addRequiredTask("add_function");
  addRequiredTask("add_postprocessor");
  addRequiredTask("add_kernel");
}

const TotalPowerAC &
HeatSourceFromTotalPowerAC::totalPower() const
{
  return component<TotalPowerAC>("power", "TotalPowerBase");
}

void
HeatSourceFromTotalPowerAC::actOnAdditionalTasks()
{
  if (_current_task == "add_function")
  {
    // This is the first task that uses the heat structure, so check the inputs that refer to it
    // here, before a later task fails with a less specific message
    totalPower();
    checkRegions();

    // Without a power shape function the source is uniform along the heat structure length
    if (!_has_psf)
    {
      const std::string class_name = "ConstantFunction";
      InputParameters params = getFactory().getValidParams(class_name);
      params.set<Real>("value") = 1. / hs().length();
      getProblem().addFunction(class_name, _power_shape_fn, params);
    }
  }
  else if (_current_task == "add_postprocessor")
  {
    const std::string class_name = "FunctionElementIntegral";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::vector<SubdomainName>>("block") = subdomains();
    params.set<FunctionName>("function") = _power_shape_fn;
    params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL};
    params.set<std::vector<OutputName>>("outputs") = {"none"};
    getProblem().addPostprocessor(class_name, _power_shape_integral_name, params);
  }
  else if (_current_task == "add_kernel")
  {
    const std::string class_name = "ADHeatStructureHeatSource";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<NonlinearVariableName>("variable") =
        hs().getPhysics().getParam<VariableName>("temperature_name");
    params.set<std::vector<SubdomainName>>("block") = subdomains();
    params.set<Real>("num_units") = hs().numberOfUnits();
    params.set<Real>("power_fraction") = getParam<Real>("power_fraction");
    params.set<FunctionName>("power_shape_function") = _power_shape_fn;
    params.set<std::vector<VariableName>>("total_power") = {totalPower().getPowerVariableName()};
    params.set<PostprocessorName>("power_shape_integral_pp") = _power_shape_integral_name;
    // The element integral of the power shape covers only the plane of the mesh, not the plate's
    // depth, so the depth still has to be applied (see HeatSourceFromTotalPower)
    if (const auto plate = dynamic_cast<const HeatStructurePlateAC *>(&hs()))
      params.set<Real>("scale") = 1.0 / plate->depth();
    const std::string kernel_name = genName(name(), "heat_src");
    getProblem().addKernel(class_name, kernel_name, params);
    THM::connectComponentControllableParam(*this, "power_fraction", class_name, kernel_name);
  }
}
