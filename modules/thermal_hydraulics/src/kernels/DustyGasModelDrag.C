//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DustyGasModelDrag.h"
#include "THMNames.h"
#include "MooseUtils.h"

registerMooseObject("ThermalHydraulicsApp", DustyGasModelDrag);

InputParameters
DustyGasModelDrag::validParams()
{
  InputParameters params = ADKernel::validParams();

  params.addRequiredCoupledVar("A", "Cross-sectional area");
  params.addParam<MaterialPropertyName>("density", THM::DENSITY, "Density property");
  params.addParam<MaterialPropertyName>("velocity", THM::VELOCITY, "Velocity property");
  params.addParam<MaterialPropertyName>("pressure", THM::PRESSURE, "Pressure property");
  params.addParam<MaterialPropertyName>(
      "dynamic_viscosity", THM::DYNAMIC_VISCOSITY, "Dynamic viscosity property");

  params.addClassDescription("Adds drags according to the Dusty Gas Model.");

  return params;
}

DustyGasModelDrag::DustyGasModelDrag(const InputParameters & parameters)
  : ADKernel(parameters),
    _A(adCoupledValue("A")),
    _rho(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("density"))),
    _vel(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("velocity"))),
    _p(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("pressure"))),
    _mu(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("dynamic_viscosity")))
{
}

ADReal
DustyGasModelDrag::computeQpResidual()
{
  const auto & A = _A[_qp];
  const auto & rho = _rho[_qp];
  const auto v = 1.0 / rho;
  const auto & vel = _vel[_qp];
  const auto & p = _p[_qp];
  const auto & mu = _mu[_qp];

  // Use permeability corresponding to Hagen-Poiseuille flow in cylinder
  const auto R = std::sqrt(A / libMesh::pi);
  const auto K = std::pow(R, 2) / 8.0;
  const auto D_visc = K * p / mu;

  const auto vel_avg = std::sqrt(8.0 * p * v / libMesh::pi);
  const auto bK = 2.0 / 3.0 * vel_avg * R;
  const Real C2 = 4.0;
  const Real C1 = 0.81 * C2;
  const auto cK1 = C1 * R / mu * std::sqrt(1.0 / (p * v));
  const auto cK2 = C2 / C1 * cK1;
  const auto D_K = bK * (1 + cK1 * p) / (1 + cK2 * p);

  const auto dpdx = -p * vel / (D_visc + D_K);
  if (MooseUtils::absoluteFuzzyEqual(vel, 0.0))
    return 0;
  else
    return -dpdx * A * _test[_i][_qp];
}
