//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FlowComponentNS2D.h"

registerMooseObject("ThermalHydraulicsApp", FlowComponentNS2D);

InputParameters
FlowComponentNS2D::validParams()
{
  InputParameters params = NSFVBase<Component2D>::validParams();

  params.addRequiredParam<bool>("is_cylindrical", "Set to true if the component is cylindrical");
  params.addParam<Real>("depth", 1.0, "Depth of the component if not cylindrical");
  params.addParam<Real>("inner_radius", 0.0, "Inner radius of the component if cylindrical");


  params.addClassDescription("Navier-Stokes flow component with a 2D generated mesh.");

  return params;
}

FlowComponentNS2D::FlowComponentNS2D(const InputParameters & parameters)
  : NSFVBase<Component2D>(parameters),

    _is_cylindrical(getParam<bool>("is_cylindrical")),
    _depth(getParam<Real>("depth")),
    _inner_radius(getParam<Real>("inner_radius"))
{
  checkCopyNSNodalVariables();
}

void FlowComponentNS2D::addRelationshipManagers(Moose::RelationshipManagerType /*input_rm_type*/)
{
  addRelationshipManagersFromParameters(getGhostParametersForRM());
}

void
FlowComponentNS2D::addVariables()
{
  addNSVariables();
  addNSInitialConditions();
}

void
FlowComponentNS2D::addMooseObjects()
{
  addNSUserObjects();
  addNSKernels();
  addNSBoundaryConditions();
  addNSMaterials();
  addNSPostprocessors();
}

void
FlowComponentNS2D::init()
{
  Component2D::init();

  processMesh();
  copyNSNodalVariables();
}

void
FlowComponentNS2D::check() const
{
  Component2D::check();

  if (isCylindrical() && !isParamSetByUser("inner_radius"))
    logError("If the component is cylindrical, then 'inner_radius' must be provided.");
  if (!isCylindrical() && !isParamSetByUser("depth"))
    logError("If the component is not cylindrical, then 'depth' must be provided.");
}
