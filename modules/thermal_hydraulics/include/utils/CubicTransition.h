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
 *  Cubic polynomial transition between two functions of one variable
 */
template <bool is_ad>
class CubicTransitionTempl : public SmoothTransition<is_ad>
{
public:
  /**
   * Constructor.
   *
   * @param[in] x_center   Center point of transition
   * @param[in] transition_width   Width of transition
   */
  CubicTransitionTempl(const GenericReal<is_ad> & x_center,
                       const GenericReal<is_ad> & transition_width);

  virtual GenericReal<is_ad> value(const GenericReal<is_ad> & x,
                                   const GenericReal<is_ad> & f1,
                                   const GenericReal<is_ad> & f2) const override;

  /**
   * Computes the derivative of the transition value
   *
   * @param[in] x   Point at which to evaluate transition
   * @param[in] df1dx   First function derivative
   * @param[in] df2dx   Second function derivative
   */
  GenericReal<is_ad> derivative(const GenericReal<is_ad> & x,
                                const GenericReal<is_ad> & df1dx,
                                const GenericReal<is_ad> & df2dx) const;

  /**
   * Initializes the polynomial coefficients
   *
   * @param[in] f1_end_value   Value of left function at left transition end point
   * @param[in] f2_end_value   Value of right function at right transition end point
   * @param[in] df1dx_end_value   Value of left function derivative at left transition end point
   * @param[in] df2dx_end_value   Value of right function derivative at right transition end point
   */
  void initialize(const GenericReal<is_ad> & f1_end_value,
                  const GenericReal<is_ad> & f2_end_value,
                  const GenericReal<is_ad> & df1dx_end_value,
                  const GenericReal<is_ad> & df2dx_end_value);

protected:
  using SmoothTransition<is_ad>::_x1;
  using SmoothTransition<is_ad>::_x2;

  // Polynomial coefficients
  GenericReal<is_ad> _A;
  GenericReal<is_ad> _B;
  GenericReal<is_ad> _C;
  GenericReal<is_ad> _D;

  /// Flag that transition has been initialized
  bool _initialized;
};

typedef CubicTransitionTempl<false> CubicTransition;
typedef CubicTransitionTempl<true> ADCubicTransition;
