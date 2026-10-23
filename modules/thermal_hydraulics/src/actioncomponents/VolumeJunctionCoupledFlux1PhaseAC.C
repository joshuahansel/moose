//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VolumeJunctionCoupledFlux1PhaseAC.h"
#include "VolumeJunction1PhaseAC.h"
#include "THMUtils.h"
#include "MooseApp.h"
#include "Factory.h"

// Registered under the classic Component's name ('VolumeJunctionCoupledFlux1Phase'); see
// FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           VolumeJunctionCoupledFlux1PhaseAC,
                           "VolumeJunctionCoupledFlux1Phase",
                           "add_kernel");
registerMooseActionAliased("ThermalHydraulicsApp",
                           VolumeJunctionCoupledFlux1PhaseAC,
                           "VolumeJunctionCoupledFlux1Phase",
                           "add_postprocessor");
registerMooseActionAliased("ThermalHydraulicsApp",
                           VolumeJunctionCoupledFlux1PhaseAC,
                           "VolumeJunctionCoupledFlux1Phase",
                           "add_transfer");
registerActionComponentAliased("ThermalHydraulicsApp",
                               VolumeJunctionCoupledFlux1PhaseAC,
                               "VolumeJunctionCoupledFlux1Phase");

InputParameters
VolumeJunctionCoupledFlux1PhaseAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addRequiredParam<std::string>("volume_junction",
                                       "Name of volume junction component on which to apply flux");
  params.addRequiredParam<RealVectorValue>("normal_from_junction",
                                           "Direction vector pointing from the junction to the "
                                           "coupled volume. This vector will be normalized.");
  params.addRequiredParam<Real>("A_coupled", "Area of the coupled surface [m^2]");
  params.addParam<MultiAppName>(
      "multi_app",
      "If provided, transfers occur with this MultiApp. The following would be transferred: mass "
      "flow rate, energy flow rate, temperature, and pressure at the junction.");
  params.addRequiredParam<std::string>("pp_suffix", "Suffix to append to post-processor names");
  params.addClassDescription(
      "Applies a flux between a VolumeJunction1Phase component and an external application.");
  return params;
}

VolumeJunctionCoupledFlux1PhaseAC::VolumeJunctionCoupledFlux1PhaseAC(const InputParameters & params)
  : ActionComponent(params),
    _volume_junction_name(getParam<std::string>("volume_junction")),
    _normal_from_junction(getParam<RealVectorValue>("normal_from_junction").unit()),
    _pp_suffix(getParam<std::string>("pp_suffix"))
{
  THM::informExperimental(*this);
  addRequiredTask("add_kernel");
  addRequiredTask("add_postprocessor");
  addRequiredTask("add_transfer");
}

const VolumeJunction1PhaseAC &
VolumeJunctionCoupledFlux1PhaseAC::volumeJunction() const
{
  for (const auto & action : _awh.getActionListByName("list_component"))
    if (action->name() == _volume_junction_name)
    {
      if (const auto * junction = dynamic_cast<const VolumeJunction1PhaseAC *>(action))
        return *junction;
      paramError("volume_junction",
                 "The component '",
                 _volume_junction_name,
                 "' is not of type 'VolumeJunction1Phase'");
    }
  paramError("volume_junction", "The component '", _volume_junction_name, "' does not exist");
}

void
VolumeJunctionCoupledFlux1PhaseAC::actOnAdditionalTasks()
{
  const std::vector<std::string> equations = {"mass", "energy"};
  const std::vector<std::string> properties = {"p", "T"};

  if (_current_task == "add_kernel")
  {
    const std::vector<std::string> vars = {"rhoV", "rhouV", "rhovV", "rhowV", "rhoEV"};
    for (const auto i : index_range(vars))
      addVolumeJunctionKernel(vars[i], i);
  }
  else if (_current_task == "add_postprocessor")
  {
    for (const auto & equation : equations)
      addFluxPostprocessor(equation);
    for (const auto & property : properties)
      addPropertyPostprocessor(property);
  }
  else if (_current_task == "add_transfer" && isParamValid("multi_app"))
  {
    for (const auto & equation : equations)
      addFluxTransfer(equation);
    for (const auto & property : properties)
      addPropertyTransfer(property);
  }
}

void
VolumeJunctionCoupledFlux1PhaseAC::addVolumeJunctionKernel(const std::string & var, unsigned int i)
{
  const auto & junction = volumeJunction();

  const std::string class_name = "VolumeJunctionCoupledFlux1PhaseKernel";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<NonlinearVariableName>("variable") = var;
  params.set<std::vector<SubdomainName>>("block") = {_volume_junction_name};
  params.set<unsigned int>("equation_index") = i;
  params.set<PostprocessorName>("pressure") = addPostprocessorSuffix("p");
  params.set<PostprocessorName>("temperature") = addPostprocessorSuffix("T");
  params.set<Real>("A_coupled") = getParam<Real>("A_coupled");
  params.set<RealVectorValue>("normal_from_junction") = _normal_from_junction;
  params.set<UserObjectName>("volume_junction_uo") = junction.junctionUserObjectName();
  params.set<UserObjectName>("numerical_flux_uo") = junction.numericalFluxName(0);
  params.set<UserObjectName>("fluid_properties") = junction.junctionFluidPropertiesName();
  getProblem().addKernel(class_name, genName(name(), var + "_kernel"), params);
}

void
VolumeJunctionCoupledFlux1PhaseAC::addFluxPostprocessor(const std::string & equation)
{
  const auto & junction = volumeJunction();

  const std::string quantity = equation + "_rate";
  const std::string class_name = "VolumeJunctionCoupledFlux1PhasePostprocessor";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<MooseEnum>("equation") = equation;
  params.set<PostprocessorName>("pressure") = addPostprocessorSuffix("p");
  params.set<PostprocessorName>("temperature") = addPostprocessorSuffix("T");
  params.set<Real>("A_coupled") = getParam<Real>("A_coupled");
  params.set<RealVectorValue>("normal_from_junction") = _normal_from_junction;
  params.set<UserObjectName>("volume_junction_uo") = junction.junctionUserObjectName();
  params.set<UserObjectName>("numerical_flux_uo") = junction.numericalFluxName(0);
  params.set<UserObjectName>("fluid_properties") = junction.junctionFluidPropertiesName();
  getProblem().addPostprocessor(class_name, addPostprocessorSuffix(quantity), params);
}

void
VolumeJunctionCoupledFlux1PhaseAC::addPropertyPostprocessor(const std::string & property)
{
  const std::string class_name = "Receiver";
  InputParameters params = getFactory().getValidParams(class_name);
  getProblem().addPostprocessor(class_name, addPostprocessorSuffix(property), params);
}

void
VolumeJunctionCoupledFlux1PhaseAC::addFluxTransfer(const std::string & equation)
{
  const std::string quantity = equation + "_rate";
  const std::string class_name = "MultiAppPostprocessorTransfer";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<MultiAppName>("to_multi_app") = getParam<MultiAppName>("multi_app");
  params.set<PostprocessorName>("from_postprocessor") = addPostprocessorSuffix(quantity);
  params.set<PostprocessorName>("to_postprocessor") = addPostprocessorSuffix(quantity);
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
  getProblem().addTransfer(class_name, genName(name(), quantity + "_transfer"), params);
}

void
VolumeJunctionCoupledFlux1PhaseAC::addPropertyTransfer(const std::string & property)
{
  const std::string class_name = "MultiAppPostprocessorTransfer";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<MultiAppName>("from_multi_app") = getParam<MultiAppName>("multi_app");
  params.set<PostprocessorName>("from_postprocessor") = addPostprocessorSuffix(property);
  params.set<PostprocessorName>("to_postprocessor") = addPostprocessorSuffix(property);
  params.set<MooseEnum>("reduction_type") = "average";
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
  getProblem().addTransfer(class_name, genName(name(), property + "_transfer"), params);
}

PostprocessorName
VolumeJunctionCoupledFlux1PhaseAC::addPostprocessorSuffix(const std::string & base_name) const
{
  return base_name + "_" + _pp_suffix;
}
