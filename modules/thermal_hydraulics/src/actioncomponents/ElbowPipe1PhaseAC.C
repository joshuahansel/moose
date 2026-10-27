//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElbowPipe1PhaseAC.h"
#include "MooseUtils.h"
#include "MooseApp.h"
#include "Factory.h"

// Registered under the classic Component's name ('ElbowPipe1Phase'); see FlowChannel1PhaseAC.C for
// why. The tasks are the same as for FlowChannel1PhaseAC.
registerMooseActionAliased("ThermalHydraulicsApp",
                           ElbowPipe1PhaseAC,
                           "ElbowPipe1Phase",
                           "add_mesh_generator");
registerMooseActionAliased("ThermalHydraulicsApp",
                           ElbowPipe1PhaseAC,
                           "ElbowPipe1Phase",
                           "init_component_physics");
registerMooseActionAliased("ThermalHydraulicsApp",
                           ElbowPipe1PhaseAC,
                           "ElbowPipe1Phase",
                           "add_functor_material");
registerMooseActionAliased("ThermalHydraulicsApp",
                           ElbowPipe1PhaseAC,
                           "ElbowPipe1Phase",
                           "add_postprocessor");
registerMooseActionAliased("ThermalHydraulicsApp",
                           ElbowPipe1PhaseAC,
                           "ElbowPipe1Phase",
                           "add_convergence");
registerMooseActionAliased("ThermalHydraulicsApp",
                           ElbowPipe1PhaseAC,
                           "ElbowPipe1Phase",
                           "add_vector_postprocessor");
registerActionComponentAliased("ThermalHydraulicsApp", ElbowPipe1PhaseAC, "ElbowPipe1Phase");

InputParameters
ElbowPipe1PhaseAC::validParams()
{
  InputParameters params = FlowChannel1PhaseAC::validParams();
  params.addRequiredParam<Real>("radius", "Radius of the pipe [m]");
  params.addRequiredParam<Real>("start_angle", "Angle at which the pipe starts [degrees]");
  params.addRequiredParam<Real>("end_angle", "Angle at which the pipe ends [degrees]");

  // The length is determined by the radius and the angles. It is required by the parent class, so
  // it is given a value and hidden.
  params.set<std::vector<Real>>("length") = {1.0};
  params.suppressParameter<std::vector<Real>>("length");

  params.addClassDescription("Bent pipe for 1-phase flow");
  return params;
}

ElbowPipe1PhaseAC::ElbowPipe1PhaseAC(const InputParameters & params)
  : ActionComponent(params),
    FlowChannel1PhaseAC(params),
    _radius(getParam<Real>("radius")),
    _start_angle(getParam<Real>("start_angle")),
    _end_angle(getParam<Real>("end_angle"))
{
  if (_start_angle > _end_angle)
    _end_angle += 360.;
  _central_angle = _end_angle - _start_angle;
  if (_central_angle > 360)
    paramError("end_angle",
               "The difference between the angle parameters 'end_angle' and 'start_angle' (",
               _central_angle,
               ") is greater than 360 degrees");

  _length = 2 * libMesh::pi * _radius * (_central_angle / 360.);
}

void
ElbowPipe1PhaseAC::addShapeMeshGenerators()
{
  // Maps the position x along the straight mesh onto the arc, as classic ElbowPipe1Phase does: the
  // angle varies linearly with the distance along the arc, and the pipe lies in the xy plane.
  const std::string alpha = "(" + Moose::stringifyExact(_start_angle) + " + " +
                            Moose::stringifyExact(_central_angle) + " * x / " +
                            Moose::stringifyExact(_length) + ") * pi / 180";

  InputParameters params = _factory.getValidParams("ParsedNodeTransformGenerator");
  params.set<MeshGeneratorName>("input") = _mg_names.back();
  params.set<ParsedFunctionExpression>("x_function") =
      Moose::stringifyExact(_radius) + " * cos(" + alpha + ")";
  params.set<ParsedFunctionExpression>("y_function") =
      Moose::stringifyExact(_radius) + " * sin(" + alpha + ")";
  params.set<ParsedFunctionExpression>("z_function") = "0";
  params.set<std::vector<std::string>>("constant_names") = {"pi"};
  params.set<std::vector<std::string>>("constant_expressions") = {
      Moose::stringifyExact(libMesh::pi)};
  _app.getMeshGeneratorSystem().addMeshGenerator(
      "ParsedNodeTransformGenerator", name() + "_elbow", params);
  _mg_names.push_back(name() + "_elbow");
}
