//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TestRegisterTaskComponent.h"

registerMooseAction("MooseApp", TestRegisterTaskComponent, "add_mesh_generator");
// TestRegisterTaskComponent is an example of ComponentPhysicsInterface
registerMooseAction("MooseApp", TestRegisterTaskComponent, "init_component_physics");
// TestRegisterTaskComponent is an example of ComponentMaterialPropertyInterface
registerMooseAction("MooseApp", TestRegisterTaskComponent, "add_material");
// TestRegisterTaskComponent is an example of ComponentInitialConditionInterface
registerMooseAction("MooseApp", TestRegisterTaskComponent, "check_integrity");
registerActionComponent("MooseApp", TestRegisterTaskComponent);

InputParameters
TestRegisterTaskComponent::validParams()
{
  InputParameters params = ActionComponent::validParams();
  return params;
}

TestRegisterTaskComponent::TestRegisterTaskComponent(const InputParameters & params)
  : ActionComponent(params)
{
  _dimension = getParam<MooseEnum>("dimension");
  addRequiredTask("add_mesh_generator");

  registerMooseAction("MooseApp", TestRegisterTaskComponent, name() + ":setup");
}

void
TestRegisterTaskComponent::addMeshGenerators()
{
  // Create the base mesh for the component using a mesh generator
  if (_dimension == 0)
    paramError("dimension", "0D cylinder not implemented");
  else if (_dimension == 1 || _dimension == 2)
  {
    InputParameters params = _factory.getValidParams("GeneratedMeshGenerator");
    params.set<MooseEnum>("dim") = _dimension;
    params.set<Real>("xmax") = {getParam<Real>("length")};
    params.set<unsigned int>("nx") = {getParam<unsigned int>("n_axial")};
    params.set<std::string>("boundary_name_prefix") = name();
    if (_dimension == 2)
    {
      params.set<Real>("ymax") = {getParam<Real>("radius")};
      if (!isParamValid("n_radial"))
        paramError("n_radial", "Should be provided for a 2D cylinder");
      params.set<unsigned int>("ny") = {getParam<unsigned int>("n_radial")};
    }
    else if (isParamValid("n_radial"))
      paramError("n_radial", "Should not be provided for a 1D cylinder");
    if (isParamValid("block"))
    {
      const auto block_name = getParam<SubdomainName>("block");
      params.set<SubdomainName>("subdomain_name") = block_name;
      _blocks.push_back(block_name);
    }
    _app.getMeshGeneratorSystem().addMeshGenerator(
        "GeneratedMeshGenerator", name() + "_base", params);
    _mg_names.push_back(name() + "_base");
  }
  else
  {
    paramError("dimension", "3D cylinder is not implemented");
    if (!isParamValid("n_radial"))
      paramError("n_radial", "Should be provided for a 3D cylinder");
    if (!isParamValid("n_azimuthal"))
      paramError("n_azimuthal", "Should be provided in 3D");
  }

  ComponentMeshTransformHelper::addMeshGenerators();
}

void
TestRegisterTaskComponent::setupComponent()
{
  if (_dimension == 2)
    _awh.getMesh()->setCoordSystem(_blocks, MultiMooseEnum("COORD_RZ"));
}

void
TestRegisterTaskComponent::checkIntegrity()
{
  ComponentInitialConditionInterface::checkIntegrity();
  ComponentBoundaryConditionInterface::checkIntegrity();
}
