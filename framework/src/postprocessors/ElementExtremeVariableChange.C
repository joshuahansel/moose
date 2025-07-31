//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementExtremeVariableChange.h"

registerMooseObject("MooseApp", ElementExtremeVariableChange);

InputParameters
ElementExtremeVariableChange::validParams()
{
  InputParameters params = ExtremeVariableChangeBase<ElementVariablePostprocessor>::validParams();
  params.addClassDescription(
      "Computes the extreme over elements of the change of a variable over a "
      "time step, nonlinear iteration, or MultiApp fixed point iteration.");
  return params;
}

ElementExtremeVariableChange::ElementExtremeVariableChange(const InputParameters & parameters)
  : ExtremeVariableChangeBase<ElementVariablePostprocessor>(parameters)
{
}
