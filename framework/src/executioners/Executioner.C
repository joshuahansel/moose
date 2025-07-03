//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// Moose includes
#include "Executioner.h"

#include "MooseApp.h"
#include "MooseMesh.h"
#include "FEProblem.h"
#include "NonlinearSystem.h"
#include "SlepcSupport.h"
#include "SecantSolve.h"
#include "SteffensenSolve.h"

// C++ includes
#include <vector>
#include <limits>

InputParameters
Executioner::validParams()
{
  InputParameters params = MooseObject::validParams();
  params += FixedPointSolve::validParams();
  params += Reporter::validParams();
  params += ReporterInterface::validParams();

  params.addParam<MooseEnum>("fixed_point_algorithm",
                             iterationMethods(),
                             "The fixed point algorithm to converge the sequence of problems.");

  params.addParam<bool>("verbose", false, "Set to true to print additional information");

  params.addDeprecatedParam<FileNameNoExtension>(
      "restart_file_base",
      "",
      "File base name used for restart",
      "Please use \"Problem/restart_file_base\" instead");

  params.registerBase("Executioner");

  params.addParamNamesToGroup("fixed_point_algorithm", "Fixed point iterations");
  params.addParamNamesToGroup("restart_file_base", "Restart");

  return params;
}

Executioner::Executioner(const InputParameters & parameters)
  : MooseObject(parameters),
    Reporter(this),
    ReporterInterface(this),
    UserObjectInterface(this),
    PostprocessorInterface(this),
    Restartable(this, "Executioners"),
    PerfGraphInterface(this),
    _fe_problem(*getCheckedPointerParam<FEProblemBase *>(
        "_fe_problem_base", "This might happen if you don't have a mesh")),
    _iteration_method(getParam<MooseEnum>("fixed_point_algorithm")),
    _restart_file_base(getParam<FileNameNoExtension>("restart_file_base")),
    _verbose(getParam<bool>("verbose"))
{
  for (const auto i : make_range(_fe_problem.numNonlinearSystems()))
    _fe_problem.getNonlinearSystemBase(i).setVerboseFlag(_verbose);

  if (!_restart_file_base.empty())
    _fe_problem.setRestartFile(_restart_file_base);

  // Instantiate the SolveObject for the MultiApp fixed point iteration algorithm
  if (_iteration_method == "picard")
    _fixed_point_solve = std::make_unique<PicardSolve>(*this);
  else if (_iteration_method == "secant")
    _fixed_point_solve = std::make_unique<SecantSolve>(*this);
  else if (_iteration_method == "steffensen")
    _fixed_point_solve = std::make_unique<SteffensenSolve>(*this);

  // Propagate the verbosity down to the problem
  if (_verbose)
    _fe_problem.setVerboseProblem(_verbose);
}

Executioner::Executioner(const InputParameters & parameters, bool)
  : MooseObject(parameters),
    Reporter(this),
    ReporterInterface(this),
    UserObjectInterface(this),
    PostprocessorInterface(this),
    Restartable(this, "Executioners"),
    PerfGraphInterface(this),
    _fe_problem(*getCheckedPointerParam<FEProblemBase *>(
        "_fe_problem_base", "This might happen if you don't have a mesh")),
    _iteration_method(getParam<MooseEnum>("fixed_point_algorithm")),
    _restart_file_base(getParam<FileNameNoExtension>("restart_file_base")),
    _verbose(getParam<bool>("verbose"))
{
  if (!_restart_file_base.empty())
    _fe_problem.setRestartFile(_restart_file_base);
}

Problem &
Executioner::problem()
{
  mooseDoOnce(mooseWarning("This method is deprecated, use feProblem() instead"));
  return _fe_problem;
}

FEProblemBase &
Executioner::feProblem()
{
  return _fe_problem;
}

PostprocessorValue &
Executioner::addAttributeReporter(const std::string & name, Real initial_value)
{
  // Get a reference to the value
  PostprocessorValue & value = declareValueByName<PostprocessorValue>(name, initial_value);

  // Create storage for the old/older values
  ReporterName r_name(this->name(), name);
  getReporterValueByName<PostprocessorValue>(r_name, 1);
  getReporterValueByName<PostprocessorValue>(r_name, 2);

  return value;
}

void
Executioner::init()
{
  checkMultiAppFixedPointParameters();
}

void
Executioner::checkMultiAppFixedPointParameters() const
{
  // const bool executing_fp_loop = _fixed_point_solve->maxFixedPointIts() > 1;
  // if (executing_fp_loop && !_fe_problem.hasMultiApps())
  //   mooseError("'fixed_point_max_its' was set to be > 1, but there are no MultiApps.");

  const std::vector<std::string> fp_params{"accept_on_max_fixed_point_iteration",
                                           "disable_fixed_point_residual_norm_check",
                                           "fixed_point_min_its",
                                           "fixed_point_max_its",
                                           "fixed_point_rel_tol",
                                           "fixed_point_abs_tol",
                                           "fixed_point_force_norms",
                                           "custom_pp",
                                           "custom_rel_tol",
                                           "custom_abs_tol",
                                           "direct_pp_value"};
  std::vector<std::string> provided_fp_params;
  for (const auto & param : fp_params)
    if (isParamSetByUser(param))
      provided_fp_params.push_back(param);

  if (provided_fp_params.size() > 0 && !executing_fp_loop)
  {
    std::ostringstream oss;
    oss << "The following MultiApp fixed point parameters were provided:\n";
    for (const auto & param : provided_fp_params)
      oss << "  " << param << "\n";
    oss << "These parameters are not allowed when 'fixed_point_max_its' <= 1.";

    mooseError(oss.str());
  }
}
