//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADDiracKernel.h"

class DiracJunction1PhaseUserObject;
class SinglePhaseFluidProperties;

/**
 * Adds the source for DiracJunction1Phase.
 */
class DiracJunction1PhaseKernel : public ADDiracKernel
{
public:
  static InputParameters validParams();

  DiracJunction1PhaseKernel(const InputParameters & parameters);

  virtual void addPoints() override;

protected:
  virtual ADReal computeQpResidual() override;

  /// Gets the PDE equation index corresponding to the solution variable
  unsigned int getEquationIndex() const;

  /// rho*A
  const ADVariableValue & _rhoA;
  /// rho*E*A
  const ADVariableValue & _rhoEA;
  /// Cross-sectional area
  const VariableValue & _A;

  /// Pressure
  const Moose::Functor<ADReal> & _p_functor;
  /// Temperature
  const Moose::Functor<ADReal> & _T_functor;
  /// Junction surface area
  const Real _A_junction;
  /// Junction point
  const Point & _point;

  /// User object to perform flux calculation
  const DiracJunction1PhaseUserObject & _dirac_junction_1phase_uo;
  /// Fluid properties
  const SinglePhaseFluidProperties & _fp;

  /// Equation index in PDE system (mass, momentum, energy)
  const unsigned int _equation_index;
};
