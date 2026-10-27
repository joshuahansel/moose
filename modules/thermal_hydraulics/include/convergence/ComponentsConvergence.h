//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IterationCountConvergence.h"

class THMProblem;
class FEProblemBase;

/**
 * Assesses convergence of all Component objects in a simulation, and of all the ActionComponents
 * if the simulation has any. Every ActionComponent must declare whether it has a convergence object
 * (see ActionComponentConvergenceInterface).
 */
class ComponentsConvergence : public IterationCountConvergence
{
public:
  static InputParameters validParams();

  ComponentsConvergence(const InputParameters & parameters);

  virtual void initialSetup() override;

  virtual void preLoop() override;
  virtual void preIteration() override;

protected:
  virtual MooseConvergenceStatus checkConvergenceInner(unsigned int iter) override;

  /// The problem
  FEProblemBase & _fe_problem;
  /// THM problem, if the problem is one (it is not if the simulation only has ActionComponents)
  const THMProblem * const _thm_problem;

  /// Convergence objects for all Components that provide one
  std::vector<Convergence *> _convergence_objects;
};
