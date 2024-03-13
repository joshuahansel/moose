//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ReferenceResidualProblem.h"
#include "FEProblem.h"
#include "NodeFaceConstraint.h"
#include "MechanicalContactConstraint.h"
#include "Convergence.h"

class AugmentedLagrangianContactProblemInterface
{
public:
  static InputParameters validParams();
  AugmentedLagrangianContactProblemInterface(const InputParameters & params);
  virtual const unsigned int & getLagrangianIterationNumber() const
  {
    return _lagrangian_iteration_number;
  }

protected:
  /// maximum mumber of augmented lagrange iterations
  const unsigned int _maximum_number_lagrangian_iterations;

  /// current augmented lagrange iteration number
  unsigned int _lagrangian_iteration_number;
};
