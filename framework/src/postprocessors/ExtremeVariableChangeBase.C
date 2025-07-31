//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ExtremeVariableChangeBase.h"

template <class T>
InputParameters
ExtremeVariableChangeBase<T>::validParams()
{
  InputParameters params = ExtremeValueBase<T>::validParams();

  MooseEnum change_over("time_step nonlinear_iteration multiapp_fp_iteration");
  params.addRequiredParam<MooseEnum>(
      "change_over", change_over, "Interval over which to compute the change");

  params.addCoupledVar("proxy_variable",
                       "The name of the variable to use to identify the location at which "
                       "the variable value should be taken; if not provided, this defaults "
                       "to the 'variable'.");

  params.addClassDescription("Computes the maximum absolute difference of a nodal variable over a "
                             "time step, nonlinear iteration, or MultiApp fixed point iteration.");
  return params;
}

template <class T>
ExtremeVariableChangeBase<T>::ExtremeVariableChangeBase(const InputParameters & parameters)
  : ExtremeValueBase<T>(parameters),
    _change_over(this->template getParam<MooseEnum>("change_over")),
    _u_old(previousValue(_change_over)),
    _proxy_variable(this->isParamValid("proxy_variable") ? this->coupledValue("proxy_variable")
                                                         : this->_u)
{
}

template <class T>
std::pair<Real, Real>
ExtremeVariableChangeBase<T>::getProxyValuePair()
{
  return std::make_pair(_proxy_variable[this->_qp], this->_u[this->_qp] - _u_old[this->_qp]);
}

template <class T>
const VariableValue &
ExtremeVariableChangeBase<T>::previousValue(const MooseEnum & /*change_over*/) const
{
  // if (change_over == "time_step")
  //   return this->coupledValueOld("variable");
  // else if (change_over == "nonlinear_iteration")
  //   return this->coupledValuePreviousNL("variable");
  // else if (change_over == "multiapp_fp_iteration")
  //   return this->coupledValuePreviousFP("variable");
  // else
  mooseError("Invalid value");
}

template class ExtremeVariableChangeBase<ElementVariablePostprocessor>;
