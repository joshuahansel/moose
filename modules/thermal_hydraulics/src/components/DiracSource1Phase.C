//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiracSource1Phase.h"
#include "FlowChannel1Phase.h"
#include "THMNames.h"

registerMooseObject("ThermalHydraulicsApp", DiracSource1Phase);

InputParameters
DiracSource1Phase::validParams()
{
  InputParameters params = Component::validParams();

  params.addRequiredParam<std::string>("flow_channel",
                                       "Name of flow channel component on which to apply source");
  params.addRequiredParam<Point>("point", "Point at which to apply source");
  params.addRequiredParam<MooseFunctorName>("mass_source_rate", "Mass source rate functor");
  params.addRequiredParam<MooseFunctorName>("energy_source_rate", "Energy source rate functor");

  params.addClassDescription("Dirac mass and energy sources for a FlowChannel1Phase.");

  return params;
}

DiracSource1Phase::DiracSource1Phase(const InputParameters & params)
  : Component(params),
    _flow_channel_name(getParam<std::string>("flow_channel")),
    _point(getParam<Point>("point"))
{
}

void
DiracSource1Phase::check() const
{
  checkComponentOfTypeExistsByName<FlowChannel1Phase>(_flow_channel_name);
}

void
DiracSource1Phase::addMooseObjects()
{
  addFunctorDiracKernel(THM::RHOA, getParam<MooseFunctorName>("mass_source_rate"));
  addFunctorDiracKernel(THM::RHOEA, getParam<MooseFunctorName>("energy_source_rate"));
}

void
DiracSource1Phase::addFunctorDiracKernel(const std::string & var, const MooseFunctorName & source_functor)
{
    const FlowChannel1Phase & flow_channel =
      getTHMProblem().getComponentByName<FlowChannel1Phase>(_flow_channel_name);

    const std::string class_name = "FunctorDiracKernel";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<NonlinearVariableName>("variable") = var;
    params.set<std::vector<SubdomainName>>("block") = flow_channel.getSubdomainNames();
    params.set<Point>("point") = _point;
    params.set<MooseFunctorName>("functor") = source_functor;
    const std::string obj_name = genName(name(), var + "_kernel");
    getTHMProblem().addDiracKernel(class_name, obj_name, params);
}
