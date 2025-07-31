//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ExtremeFunctorChange.h"

registerMooseObject("MooseApp", ExtremeFunctorChange);

InputParameters
ExtremeFunctorChange::validParams()
{
  InputParameters params = ExtremeValueBase<ElementPostprocessor>::validParams();

  MooseEnum change_over("time_step nonlinear_iteration multiapp_fp_iteration");
  params.addRequiredParam<MooseEnum>(
      "change_over", change_over, "Interval over which to compute the change");

  params.addRequiredParam<MooseFunctorName>(
      "functor", "The name of the functor for which to find the extrema");
  params.addParam<MooseFunctorName>(
      "proxy_functor",
      "The name of the functor to use to identify the location at which "
      "the functor value should be taken; if not provided, this defaults "
      "to the 'functor' parameter.");

  params.addClassDescription("Computes the maximum absolute difference of a functor over a "
                             "time step, nonlinear iteration, or MultiApp fixed point iteration.");
  return params;
}

ExtremeFunctorChange::ExtremeFunctorChange(const InputParameters & parameters)
  : ExtremeValueBase<ElementPostprocessor>(parameters),
    _ref_state(referenceState(getParam<MooseEnum>("change_over"))),
    _functor(getFunctor<Real>("functor")),
    _proxy_functor(isParamValid("proxy_functor") ? getFunctor<Real>("proxy_functor")
                                                 : getFunctor<Real>("functor")),
    _qp(0)
{
}

void
ExtremeFunctorChange::execute()
{
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    computeExtremeValue();
}

std::pair<Real, Real>
ExtremeFunctorChange::getProxyValuePair()
{
  Moose::ElemQpArg elem_qp = {_current_elem, _qp, _qrule, _q_point[_qp]};
  const Real proxy = MetaPhysicL::raw_value(_proxy_functor(elem_qp, Moose::currentState()));
  const Real change = MetaPhysicL::raw_value(
      std::abs(_functor(elem_qp, Moose::currentState()) - _functor(elem_qp, _ref_state)));
  return std::make_pair(proxy, change);
}

Moose::StateArg
ExtremeFunctorChange::referenceState(const MooseEnum & change_over) const
{
  if (change_over == "time_step")
    return Moose::oldState();
  else if (change_over == "nonlinear_iteration")
    return Moose::previousNonlinearState();
  else if (change_over == "multiapp_fp_iteration")
    return Moose::previousFixedPointState();
  else
    mooseError("Invalid value");
}
