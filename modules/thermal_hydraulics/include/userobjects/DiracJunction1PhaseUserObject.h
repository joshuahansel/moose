//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"

class ADNumericalFlux3EqnBase;

/**
 * Computes and caches the flux vector for DiracJunction1Phase.
 */
class DiracJunction1PhaseUserObject : public GeneralUserObject
{
public:
  static InputParameters validParams();

  DiracJunction1PhaseUserObject(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void residualSetup() override;
  virtual void execute() override {}
  virtual void finalize() override {}

  ADReal getPrimaryFlux(unsigned int equation_index, const std::vector<ADReal> & UL_1d, const std::vector<ADReal> & UR_1d) const;
  ADReal retrievePrimaryFlux(unsigned int equation_index) const;
  ADReal retrieveSecondaryFlux(unsigned int equation_index) const;
  bool fluxIsCached() const { return _flux_is_cached; }

protected:
  void computeFluxVectors(const std::vector<ADReal> & UL_1d, const std::vector<ADReal> & UR_1d) const;

  /// Numerical flux
  const ADNumericalFlux3EqnBase & _numerical_flux_uo;

  /// Dot product of direction from "left" to "right" with the flow channel direction
  const Real _nLR_dot_d;

  /// flux vector for the "left" cell for 1D
  mutable std::vector<ADReal> _FL_1d;
  /// flux vector for the "right" cell for 1D
  mutable std::vector<ADReal> _FR_1d;
  /// flux vector for the "left" cell for 3D
  mutable std::vector<ADReal> _FL_3d;
  /// flux vector for the "right" cell for 3D
  mutable std::vector<ADReal> _FR_3d;

  /// Flag that the flux vector has already been computed this residual evaluation
  mutable bool _flux_is_cached;
};
