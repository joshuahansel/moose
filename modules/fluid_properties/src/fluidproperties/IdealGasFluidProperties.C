//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IdealGasFluidProperties.h"

registerMooseObject("FluidPropertiesApp", IdealGasFluidProperties);

InputParameters
IdealGasFluidProperties::validParams()
{
  InputParameters params = IdealGasFluidPropertiesBase::validParams();

  params.addRangeCheckedParam<Real>("gamma", 1.4, "gamma > 1", "gamma value (cp/cv)");
  params.addParam<Real>("molar_mass", 29.0e-3, "Constant molar mass of the fluid (kg/mol)");
  params.addParam<Real>("e_ref", 0, "Reference specific internal energy [J/kg]");
  params.addParam<Real>("mu", 18.23e-6, "Dynamic viscosity, Pa.s");
  params.addParam<Real>("k", 25.68e-3, "Thermal conductivity, W/(m-K)");
  params.addParam<Real>("T_c", 0, "Critical temperature, K");
  params.addParam<Real>("rho_c", 0, "Critical density, kg/m3");
  params.addParam<Real>("e_c", 0, "Internal energy at the critical point, J/kg");

  params.addClassDescription("Ideal gas fluid properties from user-specified parameters.");

  return params;
}

IdealGasFluidProperties::IdealGasFluidProperties(const InputParameters & parameters)
  : IdealGasFluidPropertiesBase(parameters)
{
  _gamma = getParam<Real>("gamma");
  _molar_mass = getParam<Real>("molar_mass");
  _e_ref = getParam<Real>("e_ref");

  _R_specific = _R / _molar_mass;
  _cp = _gamma * _R_specific / (_gamma - 1.0);
  _cv = _cp / _gamma;

  _mu = getParam<Real>("mu");
  _k = getParam<Real>("k");

  _T_c = getParam<Real>("T_c");
  _rho_c = getParam<Real>("rho_c");
  _e_c = getParam<Real>("e_c");
}

void
IdealGasFluidProperties::initialSetupInner()
{
  // do nothing since parameters have already been set in the constructor
}
