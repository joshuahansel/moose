//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatRateExternalAppConvection.h"
#include "Function.h"

registerMooseObject("ThermalHydraulicsApp", HeatRateExternalAppConvection);

InputParameters
HeatRateExternalAppConvection::validParams()
{
  InputParameters params = SideIntegralPostprocessor::validParams();

  params.addRequiredCoupledVar("T", "Temperature");
  params.addRequiredCoupledVar("T_ext", "Temperature from external application");
  params.addRequiredCoupledVar("htc_ext", "Heat transfer coefficient from external application");
  params.addParam<FunctionName>("scale", 1.0, "Function by which to scale the heat flux");

  params.addClassDescription("Integrates a heat structure boundary convective heat "
                             "flux from an external application");

  return params;
}

HeatRateExternalAppConvection::HeatRateExternalAppConvection(const InputParameters & parameters)
  : SideIntegralPostprocessor(parameters),

    _T(coupledValue("T")),
    _T_ext(coupledValue("T_ext")),
    _htc_ext(coupledValue("htc_ext")),
    _scale_fn(getFunction("scale"))
{
}

Real
HeatRateExternalAppConvection::computeQpIntegral()
{
  return _scale_fn.value(_t, _q_point[_qp]) * _htc_ext[_qp] * (_T_ext[_qp] - _T[_qp]);
}
