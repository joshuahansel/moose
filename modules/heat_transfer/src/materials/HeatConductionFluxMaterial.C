//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatConductionFluxMaterial.h"

registerMooseObject("HeatTransferApp", HeatConductionFluxMaterial);
registerMooseObject("HeatTransferApp", ADHeatConductionFluxMaterial);

template <bool is_ad>
InputParameters
HeatConductionFluxMaterialTempl<is_ad>::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredCoupledVar("temperature", "Temperature");
  params.addRequiredParam<MaterialPropertyName>("thermal_conductivity",
                                                "Thermal conductivity material property");
  params.addRequiredParam<MaterialPropertyName>("property_name",
                                                "Name to give heat flux material property");

  // params.addClassDescription("Blah");

  return params;
}

template <bool is_ad>
HeatConductionFluxMaterialTempl<is_ad>::HeatConductionFluxMaterialTempl(
    const InputParameters & parameters)
  : Material(parameters),
    _grad_T(coupledGenericGradient<is_ad>("temperature")),
    _k(getGenericMaterialProperty<Real, is_ad>("thermal_conductivity")),
    _heat_flux(declareGenericProperty<RealVectorValue, is_ad>(
        getParam<MaterialPropertyName>("property_name")))
{
}

template <bool is_ad>
void
HeatConductionFluxMaterialTempl<is_ad>::computeQpProperties()
{
  _heat_flux[_qp] = -_k[_qp] * _grad_T[_qp];
}

template class HeatConductionFluxMaterialTempl<false>;
template class HeatConductionFluxMaterialTempl<true>;
