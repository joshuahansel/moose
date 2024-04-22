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
  : FixedPointSolve(ex), _m(5), _U_tagid(_m + 1), _f_tagid(_m + 1), _g_tagid(_m + 1)
{
  allocateStorage(true);
}

void
AndersonSolve::allocateStorage(const bool /*primary*/)
{
  for (unsigned int k = 0; k <= _m; k++)
  {
    const auto kstr = std::to_string(k);
    _U_tagid[k] = _problem.addVectorTag("anderson_U" + kstr, Moose::VECTOR_TAG_SOLUTION);
    _f_tagid[k] = _problem.addVectorTag("anderson_f" + kstr, Moose::VECTOR_TAG_SOLUTION);
    _g_tagid[k] = _problem.addVectorTag("anderson_g" + kstr, Moose::VECTOR_TAG_SOLUTION);

    _solver_sys.addVector(_U_tagid[k], false, PARALLEL);
    _solver_sys.addVector(_f_tagid[k], false, PARALLEL);
    _solver_sys.addVector(_g_tagid[k], false, PARALLEL);
  }

  _s_tagid = _problem.addVectorTag("anderson_s", Moose::VECTOR_TAG_SOLUTION);
  _y_tagid = _problem.addVectorTag("anderson_y", Moose::VECTOR_TAG_SOLUTION);

  _solver_sys.addVector(_s_tagid, false, PARALLEL);
  _solver_sys.addVector(_y_tagid, false, PARALLEL);
}

void
AndersonSolve::saveVariableValues(const bool /*primary*/)
{
  for (unsigned int k = 0; k <= _m - 1; k++)
  {
    auto & U_old = _solver_sys.getVector(_U_tagid[k]);
    auto & U_new = _solver_sys.getVector(_U_tagid[k + 1]);
    U_old = U_new;

    auto & f_old = _solver_sys.getVector(_f_tagid[k]);
    auto & f_new = _solver_sys.getVector(_f_tagid[k + 1]);
    f_old = f_new;

    auto & g_old = _solver_sys.getVector(_g_tagid[k]);
    auto & g_new = _solver_sys.getVector(_g_tagid[k + 1]);
    g_old = g_new;
  }

  auto & U = _solver_sys.getVector(_U_tagid[_m]);
  auto & solution = _solver_sys.solution();
  U = solution;
}

void
AndersonSolve::transformVariables(const std::set<dof_id_type> & /*target_dofs*/,
                                  const bool /*primary*/)
{
  auto & U = _solver_sys.solution();
  auto & f = _solver_sys.getVector(_f_tagid[_m]);
  auto & g = _solver_sys.getVector(_g_tagid[_m]);
  const auto & Uold = _solver_sys.getVector(_U_tagid[_m]);
  auto & s = _solver_sys.getVector(_s_tagid);
  auto & y = _solver_sys.getVector(_y_tagid);

  // f_k
  f = U;

  // g_k
  g = Uold;
  g -= f;

  if (_fixed_point_it > 0)
  {
    const unsigned int mk = std::min(_m, _fixed_point_it);
    DenseMatrix<Number> s_times_y(mk, mk);
    DenseVector<Number> s_times_g(mk);

    for (unsigned int i = 0; i <= mk - 1; i++)
    {
      // s_{k-mk+i} = U_{k-mk+i+1} - U_{k-mk+i}
      const auto & Unewk = _solver_sys.getVector(_U_tagid[_m - mk + i + 1]);
      const auto & Uoldk = _solver_sys.getVector(_U_tagid[_m - mk + i]);
      s = Unewk;
      s.add(-1.0, Uoldk);

      s_times_g(i) = s.dot(g);

      for (unsigned int j = 0; j <= mk - 1; j++)
      {
        // y_{k-mk+j} = g_{k-mk+j+1} - g_{k-mk+j}
        const auto & gnewk = _solver_sys.getVector(_g_tagid[_m - mk + j + 1]);
        const auto & goldk = _solver_sys.getVector(_g_tagid[_m - mk + j]);
        y = gnewk;
        y.add(-1.0, goldk);

        s_times_y(i, j) = s.dot(y);
      }
    }

    DenseVector<Number> gamma(mk);
    std::cout << "A = " << s_times_y(0, 0) << ", b = " << s_times_g(0) << std::endl;
    s_times_y.lu_solve(s_times_g, gamma);

    DenseVector<Number> alpha(mk + 1);
    alpha(0) = gamma(0);
    for (unsigned int i = 1; i < mk; i++)
      alpha(i) = gamma(i) - gamma(i - 1);
    alpha(mk) = 1.0 - gamma(mk - 1);

    U.zero();
    for (unsigned int i = 0; i <= mk; i++)
    {
      // f_{k-mk+i}
      auto & fi = _solver_sys.getVector(_f_tagid[_m - mk + i]);
      U.add(alpha(i), fi);
    }

    // s = Uold;
    // const auto & Uolder = _solver_sys.getVector(_U_tagid[_m - 1]);
    // s.add(-1.0, Uolder);

    // y = g;
    // const auto & gold = _solver_sys.getVector(_g_tagid[_m - 1]);
    // y.add(-1.0, gold);

    // const Real A = s.dot(y);
    // const Real b = s.dot(g);
    // std::cout<<"A = "<<A<<", b = "<<b<<std::endl;
    // // if (std::abs(A) < 1e-12)
    // //   mooseError("A is zero");
    // const auto gamma = b / A;
    // const auto alphaold = gamma;
    // const auto alpha = 1.0 - alphaold;

    // const auto & fold = _solver_sys.getVector(_f_tagid[_m - 1]);
    // U = fold;
    // U *= alphaold;
    // U.add(alpha, f);
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
