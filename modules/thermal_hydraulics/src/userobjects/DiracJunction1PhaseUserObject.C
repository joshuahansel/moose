//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiracJunction1PhaseUserObject.h"
#include "ADNumericalFlux3EqnBase.h"

registerMooseObject("ThermalHydraulicsApp", DiracJunction1PhaseUserObject);

InputParameters
DiracJunction1PhaseUserObject::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  params.addRequiredParam<UserObjectName>(
      "numerical_flux",
      "Numerical flux user object");

  params.addClassDescription("Computes and caches the flux vector for DiracJunction1Phase.");

  return params;
}

DiracJunction1PhaseUserObject::DiracJunction1PhaseUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters),
  _numerical_flux_uo(getUserObject<ADNumericalFlux3EqnBase>("numerical_flux")), _nLR_dot_d(1.0), _flux_is_cached(false)
{
}

void
DiracJunction1PhaseUserObject::residualSetup()
{
  _flux_is_cached = false;
}

ADReal
DiracJunction1PhaseUserObject::getPrimaryFlux(unsigned int equation_index, const std::vector<ADReal> & UL_1d, const std::vector<ADReal> & UR_1d) const
{
  if (!_flux_is_cached)
    computeFluxVectors(UL_1d, UR_1d);

  return _FL_1d[equation_index];
}

ADReal
DiracJunction1PhaseUserObject::retrievePrimaryFlux(unsigned int equation_index) const
{
  return _FR_1d[equation_index];
}

ADReal
DiracJunction1PhaseUserObject::retrieveSecondaryFlux(unsigned int equation_index) const
{
  return _FR_1d[equation_index];
}

void
DiracJunction1PhaseUserObject::computeFluxVectors(const std::vector<ADReal> & UL_1d, const std::vector<ADReal> & UR_1d) const
{
  _numerical_flux_uo.calcFlux1D(UL_1d, UR_1d, _nLR_dot_d, _FL_1d, _FR_1d, _FL_3d, _FR_3d);
  _flux_is_cached = true;
}
