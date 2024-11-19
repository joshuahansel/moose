//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SmoothTransition.h"

/**
 * Weighted transition between two functions of one variable
 */
template <bool is_ad>
class WeightedTransitionTempl : public SmoothTransition<is_ad>
{
public:
  /**
   * Constructor.
   *
   * @param[in] x_center   Center point of transition
   * @param[in] transition_width   Width of transition
   */
  WeightedTransitionTempl(const GenericReal<is_ad> & x_center,
                          const GenericReal<is_ad> & transition_width);

  virtual GenericReal<is_ad> value(const GenericReal<is_ad> & x,
                                   const GenericReal<is_ad> & f1,
                                   const GenericReal<is_ad> & f2) const override;

  /**
   * Computes the derivative of the transition value
   *
   * @param[in] x   Point at which to evaluate transition
   * @param[in] f1   First function value
   * @param[in] f2   Second function value
   * @param[in] df1dx   First function derivative
   * @param[in] df2dx   Second function derivative
   */
  GenericReal<is_ad> derivative(const GenericReal<is_ad> & x,
                                const GenericReal<is_ad> & f1,
                                const GenericReal<is_ad> & f2,
                                const GenericReal<is_ad> & df1dx,
                                const GenericReal<is_ad> & df2dx) const;

  /**
   * Computes the weight of the first function
   *
   * @param[in] x   Point at which to evaluate weight
   */
  GenericReal<is_ad> weight(const GenericReal<is_ad> & x) const;

protected:
  using SmoothTransition<is_ad>::_x1;
  using SmoothTransition<is_ad>::_x2;
};

typedef WeightedTransitionTempl<false> WeightedTransition;
typedef WeightedTransitionTempl<true> ADWeightedTransition;
