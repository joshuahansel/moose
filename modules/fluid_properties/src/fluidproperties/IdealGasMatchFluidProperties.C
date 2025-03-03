//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IdealGasMatchFluidProperties.h"
#include "SinglePhaseFluidProperties.h"

registerMooseObject("FluidPropertiesApp", IdealGasMatchFluidProperties);

InputParameters
IdealGasMatchFluidProperties::validParams()
{
  InputParameters params = IdealGasFluidPropertiesBase::validParams();

  params.addRequiredParam<UserObjectName>("fluid_properties",
                                          "SinglePhaseFluidProperties object to match");
  params.addRequiredParam<Real>("p", "Pressure at which to match properties [Pa]");
  params.addRequiredParam<Real>("T", "Temperature at which to match properties [K]");
  params.addParam<bool>("require_critical_properties",
                        false,
                        "If true, require that the reference SinglePhaseFluidProperties has "
                        "implemented the critical property methods.");

  params.addClassDescription("Ideal gas fluid properties that fit another fluid at a (p,T) state.");

  return params;
}

IdealGasMatchFluidProperties::IdealGasMatchFluidProperties(const InputParameters & parameters)
  : IdealGasFluidPropertiesBase(parameters),
    _fp(getUserObject<SinglePhaseFluidProperties>("fluid_properties")),
    _p(getParam<Real>("p")),
    _T(getParam<Real>("T"))
{
}

void
IdealGasMatchFluidProperties::initialSetupInner()
{
  _molar_mass = _fp.molarMass();
  _R_specific = _R / _molar_mass;

  _cv = _fp.cv_from_p_T(_p, _T);
  _gamma = _R_specific / _cv + 1.0;
  _cp = _gamma * _cv;

  // from h(T) = gamma * cv * T + e_ref
  _e_ref = _fp.h_from_p_T(_p, _T) - _gamma * _cv * _T;

  _mu = _fp.mu_from_p_T(_p, _T);
  _k = _fp.k_from_p_T(_p, _T);

  if (getParam<bool>("require_critical_properties"))
  {
    _T_c = _fp.criticalTemperature();
    _rho_c = _fp.criticalDensity();
    _e_c = _fp.criticalInternalEnergy();
  }
  else
  {
    _T_c = getNaN();
    _rho_c = getNaN();
    _e_c = getNaN();
  }
}
