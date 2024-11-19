//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"

/**
 * Base class for smooth transitions between two functions of one variable
 */
template <bool is_ad>
class SmoothTransition
{
public:
  /**
   * Constructor.
   *
   * @param[in] x_center   Center point of transition
   * @param[in] transition_width   Width of transition
   */
  SmoothTransition(const GenericReal<is_ad> & x_center,
                   const GenericReal<is_ad> & transition_width);

  /**
   * Updates transition center, width, and dependent quantities
   *
   * @param[in] x_center   Center point of transition
   * @param[in] transition_width   Width of transition
   */
  virtual void updateCenterAndWidth(const GenericReal<is_ad> & x_center,
                                    const GenericReal<is_ad> & transition_width);

  /**
   * Computes the transition value
   *
   * @param[in] x    Point at which to evaluate function
   * @param[in] f1   Left function
   * @param[in] f2   Right function
   */
  virtual GenericReal<is_ad> value(const GenericReal<is_ad> & x,
                                   const GenericReal<is_ad> & f1,
                                   const GenericReal<is_ad> & f2) const = 0;

  /**
   * Returns the coordinate of the left end of the transition
   */
  const GenericReal<is_ad> & leftEnd() const { return _x1; }

  /**
   * Returns the coordinate of the right end of the transition
   */
  const GenericReal<is_ad> & rightEnd() const { return _x2; }

protected:
  /// Center point of transition
  GenericReal<is_ad> _x_center;
  /// Width of transition
  GenericReal<is_ad> _transition_width;

  /// Left end point of transition
  GenericReal<is_ad> _x1;
  /// Right end point of transition
  GenericReal<is_ad> _x2;
};
