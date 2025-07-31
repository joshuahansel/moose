//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NodalExtremeVariableChange.h"

registerMooseObject("MooseApp", NodalExtremeVariableChange);

InputParameters
NodalExtremeVariableChange::validParams()
{
  InputParameters params = ExtremeVariableChangeBase<NodalVariablePostprocessor>::validParams();
  params.addClassDescription("Computes the extreme over nodes of the change of a variable over a "
                             "time step, nonlinear iteration, or MultiApp fixed point iteration.");
  return params;
}

NodalExtremeVariableChange::NodalExtremeVariableChange(const InputParameters & parameters)
  : ExtremeVariableChangeBase<NodalVariablePostprocessor>(parameters)
{
}
