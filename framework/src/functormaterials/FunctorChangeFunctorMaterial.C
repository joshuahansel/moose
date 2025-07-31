//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FunctorChangeFunctorMaterial.h"

registerMooseObject("MooseApp", FunctorChangeFunctorMaterial);
registerMooseObject("MooseApp", ADFunctorChangeFunctorMaterial);

template <bool is_ad>
InputParameters
FunctorChangeFunctorMaterialTempl<is_ad>::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.set<ExecFlagEnum>("execute_on") = {EXEC_ALWAYS};
  // params.addClassDescription(
  //     "FunctorMaterial object for declaring properties that are populated by evaluation of a "
  //     "Functor (a constant, variable, function or functor material property) objects.");
  params.addParam<std::string>("prop_name", "The name to give the functor material property");
  return params;
}

template <bool is_ad>
FunctorChangeFunctorMaterialTempl<is_ad>::FunctorChangeFunctorMaterialTempl(
    const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _functor(getFunctor<GenericReal<is_ad>>("functor")),

    _take_absolute_value(getParam<bool>("take_absolute_value")),
    _prop_name(getParam<std::string>("prop_name"))
{
  const std::set<ExecFlagType> clearance_schedule(_execute_enum.begin(), _execute_enum.end());
  addFunctorProperty<GenericReal<is_ad>>(
      _prop_name,
      [this](const auto & r, const auto & t) -> GenericReal<is_ad>
      {
        mooseAssert(t == Moose::currentState(),
                    "The functor properties defined by (AD)FunctorChangeFunctorMaterial objects "
                    "may only be evaluated at the current state.");

        const auto change = _functor(r, t) - _functor(r, _ref_state);
        if (_take_absolute_value)
          return std::abs(change);
        else
          return change;
      },
      clearance_schedule);
}

template class FunctorChangeFunctorMaterialTempl<false>;
template class FunctorChangeFunctorMaterialTempl<true>;
