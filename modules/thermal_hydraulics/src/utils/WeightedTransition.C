//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WeightedTransition.h"

template <bool is_ad>
WeightedTransitionTempl<is_ad>::WeightedTransitionTempl(const GenericReal<is_ad> & x_center,
                                                        const GenericReal<is_ad> & transition_width)
  : SmoothTransition<is_ad>(x_center, transition_width)
{
}

template <bool is_ad>
GenericReal<is_ad>
WeightedTransitionTempl<is_ad>::value(const GenericReal<is_ad> & x,
                                      const GenericReal<is_ad> & f1,
                                      const GenericReal<is_ad> & f2) const
{
  if (x <= _x1)
    return f1;
  else if (x >= _x2)
    return f2;
  else
  {
    const auto w = weight(x);
    return w * f1 + (1.0 - w) * f2;
  }
}

template <bool is_ad>
GenericReal<is_ad>
WeightedTransitionTempl<is_ad>::derivative(const GenericReal<is_ad> & x,
                                           const GenericReal<is_ad> & f1,
                                           const GenericReal<is_ad> & f2,
                                           const GenericReal<is_ad> & df1dx,
                                           const GenericReal<is_ad> & df2dx) const
{
  if (x <= _x1)
    return df1dx;
  else if (x >= _x2)
    return df2dx;
  else
  {
    const auto w = weight(x);
    const auto dwdx = -0.5 * std::sin(M_PI / (_x2 - _x1) * (x - _x1)) * M_PI / (_x2 - _x1);
    return w * df1dx + (1.0 - w) * df2dx + dwdx * (f1 - f2);
  }
}

template <bool is_ad>
GenericReal<is_ad>
WeightedTransitionTempl<is_ad>::weight(const GenericReal<is_ad> & x) const
{
  return 0.5 * (std::cos(M_PI / (_x2 - _x1) * (x - _x1)) + 1.0);
}

template class WeightedTransitionTempl<false>;
template class WeightedTransitionTempl<true>;
