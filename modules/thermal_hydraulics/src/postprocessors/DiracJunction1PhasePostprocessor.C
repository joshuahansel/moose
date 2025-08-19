//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiracJunction1PhasePostprocessor.h"
#include "DiracJunction1PhaseUserObject.h"

registerMooseObject("ThermalHydraulicsApp", DiracJunction1PhasePostprocessor);

InputParameters
DiracJunction1PhasePostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();

  params.addRequiredParam<bool>("get_primary_side", "Whether to get primary side or secondary side flux");
  MooseEnum equation("mass energy");
  params.addRequiredParam<MooseEnum>(
      "equation", equation, "Equation for which to query flux vector");
  params.addRequiredParam<UserObjectName>("dirac_junction_1phase_uo", "DiracJunction1PhaseUserObject object");

  params.addClassDescription(
      "Retrieves a mass or energy flux from a DiracJunction1PhaseUserObject.");

  return params;
}

DiracJunction1PhasePostprocessor::DiracJunction1PhasePostprocessor(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _get_primary_side(getParam<bool>("get_primary_side")),
    _equation_index(getEquationIndex(getParam<MooseEnum>("equation"))),
    _dirac_junction_1phase_uo(getUserObject<DiracJunction1PhaseUserObject>("dirac_junction_1phase_uo"))
{
}

void
DiracJunction1PhasePostprocessor::finalize()
{
  if (!_dirac_junction_1phase_uo.fluxIsCached())
    mooseError("DiracJunction1PhasePostprocessor should only be executed when the flux is cached.");

  if (_get_primary_side)
    _value = raw_value(_dirac_junction_1phase_uo.retrievePrimaryFlux(_equation_index));
  else
    _value = raw_value(_dirac_junction_1phase_uo.retrieveSecondaryFlux(_equation_index));
}

PostprocessorValue
DiracJunction1PhasePostprocessor::getValue() const
{
  return _value;
}

unsigned int
DiracJunction1PhasePostprocessor::getEquationIndex(const MooseEnum & equation) const
{
  if (equation == "mass")
    return THMVACE1D::MASS;
  else if (equation == "energy")
    return THMVACE1D::ENERGY;
  else
    mooseError("Invalid MooseEnum value.");
}
