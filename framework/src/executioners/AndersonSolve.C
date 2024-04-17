//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AndersonSolve.h"

#include "Executioner.h"
#include "FEProblemBase.h"
#include "NonlinearSystem.h"
#include "Console.h"

InputParameters
AndersonSolve::validParams()
{
  InputParameters params = FixedPointSolve::validParams();

  return params;
}

AndersonSolve::AndersonSolve(Executioner & ex) : FixedPointSolve(ex) { allocateStorage(true); }

void
AndersonSolve::allocateStorage(const bool /*primary*/)
{
  _Uold_tagid = _problem.addVectorTag("anderson_Uold", Moose::VECTOR_TAG_SOLUTION);
  _f_tagid = _problem.addVectorTag("anderson_f", Moose::VECTOR_TAG_SOLUTION);
  _fold_tagid = _problem.addVectorTag("anderson_fold", Moose::VECTOR_TAG_SOLUTION);
  _R_tagid = _problem.addVectorTag("anderson_R", Moose::VECTOR_TAG_SOLUTION);
  _Rold_tagid = _problem.addVectorTag("anderson_Rold", Moose::VECTOR_TAG_SOLUTION);

  _solver_sys.addVector(_Uold_tagid, false, PARALLEL);
  _solver_sys.addVector(_f_tagid, false, PARALLEL);
  _solver_sys.addVector(_fold_tagid, false, PARALLEL);
  _solver_sys.addVector(_R_tagid, false, PARALLEL);
  _solver_sys.addVector(_Rold_tagid, false, PARALLEL);
}

void
AndersonSolve::saveVariableValues(const bool /*primary*/)
{
  auto & U = _solver_sys.solution();
  auto & Uold = _solver_sys.getVector(_Uold_tagid);
  auto & f = _solver_sys.getVector(_f_tagid);
  auto & fold = _solver_sys.getVector(_fold_tagid);
  auto & R = _solver_sys.getVector(_R_tagid);
  auto & Rold = _solver_sys.getVector(_Rold_tagid);

  Uold = U;
  fold = f;
  Rold = R;
}

void
AndersonSolve::transformVariables(const std::set<dof_id_type> & /*target_dofs*/,
                                  const bool /*primary*/)
{
  auto & U = _solver_sys.solution();
  auto & Uold = _solver_sys.getVector(_Uold_tagid);
  auto & R = _solver_sys.getVector(_R_tagid);
  auto & Rold = _solver_sys.getVector(_Rold_tagid);
  auto & f = _solver_sys.getVector(_f_tagid);
  auto & fold = _solver_sys.getVector(_fold_tagid);

  f = U;

  R = Uold;
  R -= f;

  if (_fixed_point_it > 1)
  {
    // For rest of iteration, have Uold store U-Uold
    auto & dU = Uold;
    dU *= -1.0;
    dU.add(U);

    // For rest of iteration, have Rold store R-Rold
    auto & dR = Rold;
    dR *= -1.0;
    dR.add(R);

    const auto gamma = dU.dot(R) / dU.dot(dR);
    const auto alphaold = gamma;
    const auto alpha = 1.0 - alphaold;

    U = fold;
    U *= alphaold;
    U.add(alpha, f);
  }

  U.close();
  _solver_sys.update();
}

void
AndersonSolve::printFixedPointConvergenceHistory()
{
  _console << "\n 0 Anderson |R| = "
           << Console::outputNorm(std::numeric_limits<Real>::max(), _fixed_point_initial_norm)
           << '\n';

  Real max_norm_old = _fixed_point_initial_norm;
  for (unsigned int i = 0; i <= _fixed_point_it; ++i)
  {
    Real max_norm =
        std::max(_fixed_point_timestep_begin_norm[i], _fixed_point_timestep_end_norm[i]);
    _console << std::setw(2) << i + 1
             << " Anderson |R| = " << Console::outputNorm(max_norm_old, max_norm) << '\n';
    max_norm_old = max_norm;
  }
}
