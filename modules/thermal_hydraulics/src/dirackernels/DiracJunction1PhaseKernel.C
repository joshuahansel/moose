//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiracJunction1PhaseKernel.h"
#include "DiracJunction1PhaseUserObject.h"
#include "SinglePhaseFluidProperties.h"
#include "THMNames.h"
#include "THMIndicesVACE.h"

registerMooseObject("ThermalHydraulicsApp", DiracJunction1PhaseKernel);

InputParameters
DiracJunction1PhaseKernel::validParams()
{
  InputParameters params = ADDiracKernel::validParams();

  params.addRequiredCoupledVar("rhoA", "rho*A");
  params.addRequiredCoupledVar("rhoEA", "rho*E*A");
  params.addRequiredCoupledVar("A", "Cross-sectional area");

  params.addRequiredParam<MooseFunctorName>("pressure", "Pressure functor");
  params.addRequiredParam<MooseFunctorName>("temperature", "Temperature functor");
  params.addRequiredParam<Real>("A_junction", "Junction surface area");
  params.addRequiredParam<Point>("point", "Junction point");
  params.addRequiredParam<UserObjectName>("dirac_junction_1phase_uo", "DiracJunction1PhaseUserObject object");
  params.addRequiredParam<UserObjectName>("fluid_properties", "SinglePhaseFluidProperties object");

  params.addClassDescription("Adds the source for DiracJunction1Phase.");

  return params;
}

DiracJunction1PhaseKernel::DiracJunction1PhaseKernel(const InputParameters & parameters)
  : ADDiracKernel(parameters),
    _rhoA(adCoupledValue("rhoA")),
    _rhoEA(adCoupledValue("rhoEA")),
    _A(coupledValue("A")),
    _p_functor(getFunctor<ADReal>("pressure")),
    _T_functor(getFunctor<ADReal>("temperature")),
    _A_junction(getParam<Real>("A_junction")),
    _point(getParam<Point>("point")),
    _dirac_junction_1phase_uo(getUserObject<DiracJunction1PhaseUserObject>("dirac_junction_1phase_uo")),
    _fp(getUserObject<SinglePhaseFluidProperties>("fluid_properties")),
    _equation_index(getEquationIndex())
{
}

void
DiracJunction1PhaseKernel::addPoints()
{
  addPoint(_point);
}

ADReal
DiracJunction1PhaseKernel::computeQpResidual()
{
  std::vector<ADReal> UL_1d(THMVACE1D::N_FLUX_INPUTS);
  UL_1d[THMVACE1D::RHOA] = _rhoA[_qp] / _A[_qp] * _A_junction;
  UL_1d[THMVACE1D::RHOUA] = 0;
  UL_1d[THMVACE1D::RHOEA] = _rhoEA[_qp] / _A[_qp] * _A_junction;
  UL_1d[THMVACE1D::AREA] = _A_junction;

  const Moose::ElemQpArg space_arg = {_current_elem, _qp, _qrule, _current_point};
  const auto pR = _p_functor(space_arg, Moose::currentState());
  const auto TR = _T_functor(space_arg, Moose::currentState());
  const auto rhoR = _fp.rho_from_p_T(pR, TR);
  const auto ER = _fp.e_from_p_T(pR, TR);

  std::vector<ADReal> UR_1d(THMVACE1D::N_FLUX_INPUTS);
  UR_1d[THMVACE1D::RHOA] = rhoR * _A_junction;
  UR_1d[THMVACE1D::RHOUA] = 0;
  UR_1d[THMVACE1D::RHOEA] = rhoR * ER * _A_junction;
  UR_1d[THMVACE1D::AREA] = _A_junction;

  const auto flux = _dirac_junction_1phase_uo.getPrimaryFlux(_equation_index, UL_1d, UR_1d);

  return -_test[_i][_qp] * flux;
}

unsigned int
DiracJunction1PhaseKernel::getEquationIndex() const
{
  if (_var.name() == THM::RHOA)
    return THMVACE1D::MASS;
  else if (_var.name() == THM::RHOUA)
    return THMVACE1D::MOMENTUM;
  else if (_var.name() == THM::RHOEA)
    return THMVACE1D::ENERGY;
  else
    mooseError("The kernel variable must be 'rhoA', 'rhouA', or 'rhoEA'.");
}
