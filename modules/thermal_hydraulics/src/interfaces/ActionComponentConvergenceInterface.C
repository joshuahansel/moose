//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ActionComponentConvergenceInterface.h"
#include "FEProblemBase.h"
#include "Factory.h"

#include <limits>

namespace THM
{
void
addNonlinearStepFunctorMaterial(FEProblemBase & problem,
                                Factory & factory,
                                const std::string & name,
                                const std::vector<SubdomainName> & blocks,
                                const std::string & functor_name,
                                const std::string & property,
                                bool functor_is_ad)
{
  const std::string class_name =
      functor_is_ad ? "ADFunctorChangeFunctorMaterial" : "FunctorChangeFunctorMaterial";
  InputParameters params = factory.getValidParams(class_name);
  params.set<std::vector<SubdomainName>>("block") = blocks;
  params.set<MooseFunctorName>("functor") = functor_name;
  params.set<MooseEnum>("change_over") = "nonlinear";
  params.set<std::string>("prop_name") = property;
  params.set<bool>("take_absolute_value") = true;
  problem.addFunctorMaterial(class_name, name, params);
}

void
addMaximumFunctorPostprocessor(FEProblemBase & problem,
                               Factory & factory,
                               const std::string & functor_name,
                               const std::string & pp_name,
                               Real normalization,
                               const std::vector<SubdomainName> & subdomains)
{
  const std::string class_name = "ElementExtremeFunctorValue";
  InputParameters params = factory.getValidParams(class_name);
  params.set<std::vector<SubdomainName>>("block") = subdomains;
  params.set<MooseEnum>("value_type") = "max";
  params.set<MooseFunctorName>("functor") = functor_name;
  params.set<Real>("scale") = 1.0 / normalization;
  params.set<ExecFlagEnum>("execute_on") = EXEC_NONLINEAR_CONVERGENCE;
  params.set<std::vector<OutputName>>("outputs") = {"none"};
  problem.addPostprocessor(class_name, pp_name, params);
}

void
addMultiPostprocessorConvergence(FEProblemBase & problem,
                                 Factory & factory,
                                 const ConvergenceName & name,
                                 const std::vector<PostprocessorName> & postprocessors,
                                 const std::vector<std::string> & descriptions,
                                 const std::vector<Real> & tolerances)
{
  const std::string class_name = "MultiPostprocessorConvergence";
  InputParameters params = factory.getValidParams(class_name);
  params.set<std::vector<PostprocessorName>>("postprocessors") = postprocessors;
  params.set<std::vector<std::string>>("descriptions") = descriptions;
  params.set<std::vector<Real>>("tolerances") = tolerances;
  params.set<unsigned int>("min_iterations") = 1;
  params.set<unsigned int>("max_iterations") = std::numeric_limits<unsigned int>::max();
  problem.addConvergence(class_name, name, params);
}
}
