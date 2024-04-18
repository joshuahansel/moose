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

AndersonSolve::AndersonSolve(Executioner & ex)
  : FixedPointSolve(ex), _m(1), _U_tagid(_m + 1), _f_tagid(_m + 1), _g_tagid(_m + 1)
{
  allocateStorage(true);
}

void
AndersonSolve::allocateStorage(const bool /*primary*/)
{
  for (unsigned int k = 0; k < _m + 1; k++)
  {
    const auto kstr = std::to_string(k);
    _U_tagid[k] = _problem.addVectorTag("anderson_U" + kstr, Moose::VECTOR_TAG_SOLUTION);
    _f_tagid[k] = _problem.addVectorTag("anderson_f" + kstr, Moose::VECTOR_TAG_SOLUTION);
    _g_tagid[k] = _problem.addVectorTag("anderson_g" + kstr, Moose::VECTOR_TAG_SOLUTION);

    _solver_sys.addVector(_U_tagid[k], false, PARALLEL);
    _solver_sys.addVector(_f_tagid[k], false, PARALLEL);
    _solver_sys.addVector(_g_tagid[k], false, PARALLEL);
  }
}

void
AndersonSolve::saveVariableValues(const bool /*primary*/)
{
  for (unsigned int k = 0; k < _m; k++)
  {
    auto & U_old = _solver_sys.getVector(_U_tagid[k]);
    auto & U_new = (k == _m - 1) ? _solver_sys.solution() : _solver_sys.getVector(_U_tagid[k + 1]);
    U_old = U_new;

    auto & f_old = _solver_sys.getVector(_f_tagid[k]);
    auto & f_new = _solver_sys.getVector(_f_tagid[k + 1]);
    f_old = f_new;

    auto & g_old = _solver_sys.getVector(_g_tagid[k]);
    auto & g_new = _solver_sys.getVector(_g_tagid[k + 1]);
    g_old = g_new;
  }
}

void
AndersonSolve::transformVariables(const std::set<dof_id_type> & /*target_dofs*/,
                                  const bool /*primary*/)
{
  auto & U = _solver_sys.solution();
  auto & Uold = _solver_sys.getVector(_U_tagid[0]);
  auto & g = _solver_sys.getVector(_g_tagid[1]);
  auto & gold = _solver_sys.getVector(_g_tagid[0]);
  auto & f = _solver_sys.getVector(_f_tagid[1]);
  auto & fold = _solver_sys.getVector(_f_tagid[0]);

  f = U;

  g = Uold;
  g -= f;

  if (_fixed_point_it > 1)
  {
    // For rest of iteration, have Uold store U-Uold
    auto & dU = Uold;
    dU *= -1.0;
    dU.add(U);

    // For rest of iteration, have Rold store R-Rold
    auto & dg = gold;
    dg *= -1.0;
    dg.add(g);

    const auto gamma = dU.dot(g) / dU.dot(dg);
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
