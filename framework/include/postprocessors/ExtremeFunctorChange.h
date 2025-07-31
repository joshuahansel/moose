//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ExtremeValueBase.h"
#include "ElementPostprocessor.h"

/**
 * Computes the extreme of the absolute value of the change of a functor over a time step,
 * nonlinear iteration, or MultiApp fixed point iteration.
 */
class ExtremeFunctorChange : public ExtremeValueBase<ElementPostprocessor>
{
public:
  static InputParameters validParams();

  ExtremeFunctorChange(const InputParameters & parameters);

  virtual void execute() override;

protected:
  virtual std::pair<Real, Real> getProxyValuePair() override;

  /// Returns the state argument to use for the reference value
  Moose::StateArg referenceState(const MooseEnum & change_over) const;

  /// Reference value state argument
  const Moose::StateArg _ref_state;

  /// Functor to search the extrema for
  const Moose::Functor<Real> & _functor;
  /// Proxy functor used to find the extreme location
  const Moose::Functor<Real> & _proxy_functor;

  /// Quadrature point index
  unsigned int _qp;
};
