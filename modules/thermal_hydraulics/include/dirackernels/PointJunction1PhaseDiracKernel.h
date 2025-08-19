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

/**
 */
class PointJunction1PhaseDiracKernel : public ADDiracKernel
{
public:
  static InputParameters validParams();

  PointJunction1PhaseDiracKernel(const InputParameters & parameters);

  virtual void residualSetup() override;
  virtual void jacobianSetup() override;
  virtual void addPoints() override;

protected:
  virtual ADReal computeQpResidual() override;

  /**
   * Gets the PDE equation index corresponding to the solution variable
   */
  unsigned int getEquationIndex() const;

  /// Pressure
  const Moose::Functor<ADReal> & _p;
  /// Temperature
  const Moose::Functor<ADReal> & _T;
  /// Junction point
  const Point & _point;

  /// Equation index in PDE system (mass, momentum, energy)
  const unsigned int _equation_index;

  /// Flag that the flux vector has already been computed and cached this residual/Jacobian
  bool _flux_cached;
};
