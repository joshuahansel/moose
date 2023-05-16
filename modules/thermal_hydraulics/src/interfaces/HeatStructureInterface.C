//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatStructureInterface.h"
#include "HeatConductionModel.h"
#include "MeshComponent.h"

InputParameters
HeatStructureInterface::validParams()
{
  InputParameters params = emptyInputParameters();

  params.addParam<FunctionName>("initial_T", "Initial temperature [K]");
  params.addParam<Real>(
      "scaling_factor_temperature", 1.0, "Scaling factor for solid temperature variable.");

  params.addPrivateParam<std::string>("component_type", "heat_struct");

  return params;
}

HeatStructureInterface::HeatStructureInterface(MeshComponent * mesh_component)
  : _mesh_component_hsi(*mesh_component)
{
}

std::shared_ptr<HeatConductionModel>
HeatStructureInterface::buildModel()
{
  auto & factory = _mesh_component_hsi.getMooseApp().getFactory();

  const std::string class_name = "HeatConductionModel";
  InputParameters params = factory.getValidParams(class_name);
  params.set<THMProblem *>("_thm_problem") = &_mesh_component_hsi.getTHMProblem();
  params.set<HeatStructureInterface *>("_hs") = this;
  params.applyParameters(_mesh_component_hsi.parameters());
  return factory.create<HeatConductionModel>(class_name, _mesh_component_hsi.name(), params, 0);
}

void
HeatStructureInterface::init()
{
  _hc_model = buildModel();
}

void
HeatStructureInterface::check() const
{
  auto & moose_app = _mesh_component_hsi.getMooseApp();
  bool ics_set = _mesh_component_hsi.getTHMProblem().hasInitialConditionsFromFile() ||
                 _mesh_component_hsi.isParamValid("initial_T");
  if (!ics_set && !moose_app.isRestarting())
    _mesh_component_hsi.logError("Missing initial condition for temperature.");
}

void
HeatStructureInterface::addVariables()
{
  _hc_model->addVariables();
  if (_mesh_component_hsi.isParamValid("initial_T"))
    _hc_model->addInitialConditions();
}

void
HeatStructureInterface::addMooseObjects()
{
  _hc_model->addHeatEquation();
}

FunctionName
HeatStructureInterface::getInitialT() const
{
  if (_mesh_component_hsi.isParamValid("initial_T"))
    return _mesh_component_hsi.getParam<FunctionName>("initial_T");
  else
    _mesh_component_hsi.mooseError(_mesh_component_hsi.name(),
                                   ": The parameter 'initial_T' was requested but not supplied");
}
