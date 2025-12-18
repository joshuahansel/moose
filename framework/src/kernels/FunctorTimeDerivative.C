//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FunctorTimeDerivative.h"

registerMooseObject("MooseApp", FunctorTimeDerivative);

InputParameters
FunctorTimeDerivative::validParams()
{
  InputParameters params = FunctorKernel::validParams();

  params.addClassDescription("Adds a time derivative from a functor.");

  params.addRequiredParam<bool>(
      "evaluate_dot",
      "If set to 'true', the dot() operator will be evaluated on the functor instead of the "
      "evaluate() operator. Thus 'true' corresponds to the provided functor being the quantity y, "
      "and 'false' corresponds to the provided functor being dy/dt.");

  params.set<MultiMooseEnum>("vector_tags") = "time";
  params.set<MultiMooseEnum>("matrix_tags") = "system time";

  // Only the 'add' mode is supported in the time derivative
  params.suppressParameter<MooseEnum>("mode");

  return params;
}

FunctorTimeDerivative::FunctorTimeDerivative(const InputParameters & parameters)
  : FunctorKernel(parameters), _evaluate_dot(getParam<bool>("evaluate_dot"))
{
}

ADReal
FunctorTimeDerivative::precomputeQpResidual()
{
  const Moose::ElemQpArg space_arg = {_current_elem, _qp, _qrule, _q_point[_qp]};
  if (_evaluate_dot)
    return _sign * _functor.dot(space_arg, Moose::currentState());
  else
    return _sign * _functor(space_arg, Moose::currentState());
}
