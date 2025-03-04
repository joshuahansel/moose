//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TwoPhaseFluidProperties.h"
#include "NaNInterface.h"

class TwoPhaseFluidProperties;

/**
 * Two-phase fluid properties with stiffened gas for liquid, ideal gas for vapor that matches
 * another two-phase fluid properties at the saturation state at a given temperature.
 */
class StiffenedGasMatchTwoPhaseFluidProperties : public TwoPhaseFluidProperties, public NaNInterface
{
public:
  static InputParameters validParams();

  StiffenedGasMatchTwoPhaseFluidProperties(const InputParameters & parameters);

  virtual Real p_critical() const override;
  virtual Real T_triple() const override;
  virtual Real L_fusion() const override;
  virtual Real T_sat(Real p) const override;
  virtual Real p_sat(Real T) const override;
  virtual Real dT_sat_dp(Real p) const override;
  virtual Real sigma_from_T(Real T) const override;
  virtual Real dsigma_dT_from_T(Real T) const override;

  virtual bool supportsPhaseChange() const override { return true; }

protected:
  /// Fluid properties to match
  const TwoPhaseFluidProperties & _fp_2phase;

  /// Temperature of state to match
  const Real _T;
  /// Pressure of state to match
  const Real _p;
};
