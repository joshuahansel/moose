//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiracJunction1Phase.h"
#include "FlowChannel1Phase.h"
#include "THMNames.h"

registerMooseObject("ThermalHydraulicsApp", DiracJunction1Phase);

InputParameters
DiracJunction1Phase::validParams()
{
  InputParameters params = Component::validParams();

  params.addRequiredParam<std::string>("flow_channel",
                                       "Name of flow channel component on which to apply source");
  params.addRequiredParam<Point>("point", "Point at which to apply source");
  params.addRequiredParam<Real>("A_junction", "Junction surface area [m^2]");
  params.addRequiredParam<MooseFunctorName>("pressure", "Pressure functor [Pa]");
  params.addRequiredParam<MooseFunctorName>("temperature", "Temperature functor [K]");

  params.addClassDescription("Dirac flux source for a FlowChannel1Phase.");

  return params;
}

DiracJunction1Phase::DiracJunction1Phase(const InputParameters & params)
  : Component(params),
    _flow_channel_name(getParam<std::string>("flow_channel")),
    _point(getParam<Point>("point")),
    _uo_name(genName(name(), "uo"))
{
}

void
DiracJunction1Phase::check() const
{
  checkComponentOfTypeExistsByName<FlowChannel1Phase>(_flow_channel_name);
}

void
DiracJunction1Phase::addMooseObjects()
{
  const FlowChannel1Phase & flow_channel =
      getTHMProblem().getComponentByName<FlowChannel1Phase>(_flow_channel_name);

  // user object
  {
    const std::string class_name = "DiracJunction1PhaseUserObject";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<UserObjectName>("numerical_flux") = flow_channel.getNumericalFluxUserObjectName();
    getTHMProblem().addUserObject(class_name, _uo_name, params);
  }

  // kernels
  addDiracJunction1PhaseKernel(THM::RHOA);
  addDiracJunction1PhaseKernel(THM::RHOEA);
}

void
DiracJunction1Phase::addDiracJunction1PhaseKernel(const std::string & var)
{
    const FlowChannel1Phase & flow_channel =
      getTHMProblem().getComponentByName<FlowChannel1Phase>(_flow_channel_name);

    const std::string class_name = "DiracJunction1PhaseKernel";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<NonlinearVariableName>("variable") = var;
    params.set<std::vector<SubdomainName>>("block") = flow_channel.getSubdomainNames();
    params.set<std::vector<VariableName>>("rhoA") = {THM::RHOA};
    params.set<std::vector<VariableName>>("rhoEA") = {THM::RHOEA};
    params.set<std::vector<VariableName>>("A") = {THM::AREA};
    params.set<Real>("A_junction") = getParam<Real>("A_junction");
    params.set<Point>("point") = _point;
    params.set<MooseFunctorName>("pressure") = getParam<MooseFunctorName>("pressure");
    params.set<MooseFunctorName>("temperature") = getParam<MooseFunctorName>("temperature");
    params.set<UserObjectName>("dirac_junction_1phase_uo") = _uo_name;
    params.set<UserObjectName>("fluid_properties") = flow_channel.getFluidPropertiesName();
    const std::string obj_name = genName(name(), var + "_kernel");
    getTHMProblem().addDiracKernel(class_name, obj_name, params);
}

