//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"

class DiracJunction1PhaseUserObject;

/**
 * Retrieves a mass or energy flux from a DiracJunction1PhaseUserObject.
 */
class DiracJunction1PhasePostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  DiracJunction1PhasePostprocessor(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual void finalize() override;
  virtual PostprocessorValue getValue() const override;

protected:
  /// Gets the requested PDE equation index
  unsigned int getEquationIndex(const MooseEnum & equation) const;

  /// Whether to get primary side or secondary side flux
  const bool _get_primary_side;

  /// Index within flux vector to query
  const unsigned int _equation_index;

  /// User object to perform flux calculation
  const DiracJunction1PhaseUserObject & _dirac_junction_1phase_uo;

  /// Value of this PP
  Real _value;
};
