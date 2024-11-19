//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CubicTransition.h"
#include "MooseError.h"
#include "ADReal.h"

#include "libmesh/dense_vector.h"
#include "DenseMatrix.h"

template <bool is_ad>
CubicTransitionTempl<is_ad>::CubicTransitionTempl(const GenericReal<is_ad> & x_center,
                                                  const GenericReal<is_ad> & transition_width)
  : SmoothTransition<is_ad>(x_center, transition_width),

    _A(0.0),
    _B(0.0),
    _C(0.0),
    _D(0.0),

    _needs_init(true)
{
}

template <bool is_ad>
void
CubicTransitionTempl<is_ad>::updateCenterAndWidth(const GenericReal<is_ad> & x_center,
                                                  const GenericReal<is_ad> & transition_width)
{
  SmoothTransition<is_ad>::updateCenterAndWidth(x_center, transition_width);
  _needs_init = true;
}

template <bool is_ad>
void
CubicTransitionTempl<is_ad>::initialize(const GenericReal<is_ad> & f1_end_value,
                                        const GenericReal<is_ad> & f2_end_value,
                                        const GenericReal<is_ad> & df1dx_end_value,
                                        const GenericReal<is_ad> & df2dx_end_value)
{
  // compute cubic polynomial coefficients

  DenseMatrix<GenericReal<is_ad>> mat(4, 4);

  mat(0, 0) = std::pow(_x1, 3);
  mat(0, 1) = std::pow(_x1, 2);
  mat(0, 2) = _x1;
  mat(0, 3) = 1.0;

  mat(1, 0) = std::pow(_x2, 3);
  mat(1, 1) = std::pow(_x2, 2);
  mat(1, 2) = _x2;
  mat(1, 3) = 1.0;

  mat(2, 0) = 3.0 * std::pow(_x1, 2);
  mat(2, 1) = 2.0 * _x1;
  mat(2, 2) = 1.0;
  mat(2, 3) = 0.0;

  mat(3, 0) = 3.0 * std::pow(_x2, 2);
  mat(3, 1) = 2.0 * _x2;
  mat(3, 2) = 1.0;
  mat(3, 3) = 0.0;

  DenseVector<GenericReal<is_ad>> rhs(4);
  rhs(0) = f1_end_value;
  rhs(1) = f2_end_value;
  rhs(2) = df1dx_end_value;
  rhs(3) = df2dx_end_value;

  DenseVector<GenericReal<is_ad>> coefs(4);
  mat.lu_solve(rhs, coefs);

  _A = coefs(0);
  _B = coefs(1);
  _C = coefs(2);
  _D = coefs(3);

  _needs_init = false;
}

template <bool is_ad>
GenericReal<is_ad>
CubicTransitionTempl<is_ad>::value(const GenericReal<is_ad> & x,
                                   const GenericReal<is_ad> & f1,
                                   const GenericReal<is_ad> & f2) const
{
  mooseAssert(!_needs_init, "initialize() must be called.");

  if (x <= _x1)
    return f1;
  else if (x >= _x2)
    return f2;
  else
    return _A * std::pow(x, 3) + _B * std::pow(x, 2) + _C * x + _D;
}

template <bool is_ad>
GenericReal<is_ad>
CubicTransitionTempl<is_ad>::derivative(const GenericReal<is_ad> & x,
                                        const GenericReal<is_ad> & df1dx,
                                        const GenericReal<is_ad> & df2dx) const
{
  mooseAssert(!_needs_init, "initialize() must be called.");

  if (x <= _x1)
    return df1dx;
  else if (x >= _x2)
    return df2dx;
  else
    return 3.0 * _A * std::pow(x, 2) + 2.0 * _B * x + _C;
}

template class CubicTransitionTempl<false>;
template class CubicTransitionTempl<true>;
