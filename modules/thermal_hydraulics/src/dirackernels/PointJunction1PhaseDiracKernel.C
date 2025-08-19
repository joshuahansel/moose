//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PointJunction1PhaseDiracKernel.h"
#include "THMNames.h"
#include "THMIndicesVACE.h"

registerMooseObject("ThermalHydraulicsApp", PointJunction1PhaseDiracKernel);

InputParameters
PointJunction1PhaseDiracKernel::validParams()
{
  InputParameters params = ADDiracKernel::validParams();

  params.addRequiredParam<MooseFunctorName>("pressure", "Pressure functor");
  params.addRequiredParam<MooseFunctorName>("temperature", "Temperature functor");
  params.addRequiredParam<Point>("point", "Junction point");

  // params.addClassDescription("Computes a dirac source using a functor.");

  return params;
}

PointJunction1PhaseDiracKernel::PointJunction1PhaseDiracKernel(const InputParameters & parameters)
  : ADDiracKernel(parameters),
    _p(getFunctor<ADReal>("pressure")),
    _T(getFunctor<ADReal>("temperature")),
    _point(getParam<Point>("point")),
    _equation_index(getEquationIndex()),
    _flux_cached(false)
{
}

void
PointJunction1PhaseDiracKernel::residualSetup()
{
  std::cout << "residualSetup" << std::endl;
  _flux_cached = false;
}

void
PointJunction1PhaseDiracKernel::jacobianSetup()
{
  std::cout << "jacobianSetup" << std::endl;
  _flux_cached = false;
}

void
PointJunction1PhaseDiracKernel::addPoints()
{
  addPoint(_point);
}

ADReal
PointJunction1PhaseDiracKernel::computeQpResidual()
{
  std::cout << "elem=" << _current_elem->id() << ", qp=" << _qp << ", var=" << _var.name()
            << std::endl;
  if (!_flux_cached)
  {
    std::cout << "  recomputing flux" << std::endl;
    // compute flux
    _flux_cached = true;
  }
  // const Moose::ElemQpArg space_arg = {_current_elem, _qp, _qrule, _current_point};
  // return -_test[_i][_qp] * _functor(space_arg, Moose::currentState());
  return 0;
}

unsigned int
PointJunction1PhaseDiracKernel::getEquationIndex() const
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
