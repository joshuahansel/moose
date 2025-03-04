//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StiffenedGasMatchTwoPhaseFluidProperties.h"
#include "TwoPhaseFluidProperties.h"
#include "StiffenedGasMatchFluidProperties.h"
#include "IdealGasMatchFluidProperties.h"

registerMooseObject("FluidPropertiesApp", StiffenedGasMatchTwoPhaseFluidProperties);

InputParameters
StiffenedGasMatchTwoPhaseFluidProperties::validParams()
{
  InputParameters params = TwoPhaseFluidProperties::validParams();
  params += NaNInterface::validParams();

  params.addRequiredParam<UserObjectName>("fluid_properties",
                                          "TwoPhaseFluidProperties object to match");
  params.addRequiredParam<Real>("T", "Temperature at which to match properties [K]");

  params.addClassDescription(
      "Two-phase fluid properties with stiffened gas for liquid, ideal gas for vapor that matches "
      "another two-phase fluid properties at the saturation state at a given temperature.");

  return params;
}

StiffenedGasMatchTwoPhaseFluidProperties::StiffenedGasMatchTwoPhaseFluidProperties(
    const InputParameters & parameters)
  : TwoPhaseFluidProperties(parameters),
    NaNInterface(this),
    _fp_2phase(getUserObject<TwoPhaseFluidProperties>("fluid_properties")),
    _T(getParam<Real>("T")),
    _p(_fp_2phase.p_sat(_T))
{
  if (_tid == 0)
  {
    std::string class_name = "StiffenedGasMatchFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    params.set<UserObjectName>("fluid_properties") = _fp_2phase.getLiquidName();
    params.set<Real>("p") = _p;
    params.set<Real>("T") = _T;
    params.set<MooseEnum>("emit_on_nan") = getParam<MooseEnum>("emit_on_nan");
    params.set<bool>("allow_nonphysical_states") = false;
    _fe_problem.addUserObject(class_name, _liquid_name, params);
  }
  _fp_liquid = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(_liquid_name, _tid);

  if (_tid == 0)
  {
    std::string class_name = "IdealGasMatchFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    params.set<UserObjectName>("fluid_properties") = _fp_2phase.getVaporName();
    params.set<Real>("p") = _p;
    params.set<Real>("T") = _T;
    params.set<MooseEnum>("emit_on_nan") = getParam<MooseEnum>("emit_on_nan");
    _fe_problem.addUserObject(class_name, _vapor_name, params);
  }
  _fp_vapor = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(_vapor_name, _tid);
}

Real
StiffenedGasMatchTwoPhaseFluidProperties::p_critical() const
{
  return _fp_2phase.p_critical();
}

Real
StiffenedGasMatchTwoPhaseFluidProperties::T_triple() const
{
  return _fp_2phase.T_triple();
}

Real
StiffenedGasMatchTwoPhaseFluidProperties::L_fusion() const
{
  return _fp_2phase.L_fusion();
}

Real
StiffenedGasMatchTwoPhaseFluidProperties::T_sat(Real p) const
{
  return _fp_2phase.T_sat(p);
}

Real
StiffenedGasMatchTwoPhaseFluidProperties::p_sat(Real T) const
{
  return _fp_2phase.p_sat(T);
}

Real
StiffenedGasMatchTwoPhaseFluidProperties::dT_sat_dp(Real p) const
{
  return _fp_2phase.dT_sat_dp(p);
}

Real
StiffenedGasMatchTwoPhaseFluidProperties::sigma_from_T(Real T) const
{
  return _fp_2phase.sigma_from_T(T);
}

Real
StiffenedGasMatchTwoPhaseFluidProperties::dsigma_dT_from_T(Real T) const
{
  return _fp_2phase.dsigma_dT_from_T(T);
}
