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

/**
 * Base class for computing the extreme of the change of a variable over a time step,
 * nonlinear iteration, or MultiApp fixed point iteration.
 */
template <class T>
class ExtremeVariableChangeBase : public ExtremeValueBase<T>
{
public:
  static InputParameters validParams();

  ExtremeVariableChangeBase(const InputParameters & parameters);

protected:
  virtual std::pair<Real, Real> getProxyValuePair() override;

  /// Gets the previous (old, previous NL, or previous FP) variable value
  const VariableValue & previousValue(const MooseEnum & change_over) const;

  /// Interval over which to compute change
  const MooseEnum & _change_over;

  /// The previous (old, previous NL, or previous FP) variable value
  const VariableValue & _u_old;

  /// Proxy variable used to find the extreme location
  const VariableValue & _proxy_variable;
};
