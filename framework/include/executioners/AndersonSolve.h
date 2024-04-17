//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FixedPointSolve.h"

class AndersonSolve : public FixedPointSolve
{
public:
  AndersonSolve(Executioner & ex);

  static InputParameters validParams();

  virtual void allocateStorage(const bool primary) override final;

private:
  virtual void saveVariableValues(const bool primary) override final;
  virtual void savePostprocessorValues(const bool /*primary*/) override final {}
  virtual bool useFixedPointAlgorithmUpdateInsteadOfPicard(const bool /*primary*/) override final
  {
    return true;
  };
  virtual void transformPostprocessors(const bool /*primary*/) override final{};
  virtual void transformVariables(const std::set<dof_id_type> & transformed_dofs,
                                  const bool primary) override final;
  virtual void printFixedPointConvergenceHistory() override final;

  TagID _Uold_tagid;
  TagID _f_tagid;
  TagID _fold_tagid;
  TagID _R_tagid;
  TagID _Rold_tagid;
};
