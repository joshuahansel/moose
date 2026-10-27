//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"

#include <optional>
#include <string>
#include <vector>

class FEProblemBase;
class Factory;

namespace THM
{
/**
 * Adds a functor material that computes the absolute change of a functor over a nonlinear
 * iteration, for ActionComponents; the counterpart of Component::addNonlinearStepFunctorMaterial()
 *
 * @param[in] name   Name of the functor material
 * @param[in] blocks   Blocks of the functor material
 * @param[in] functor_name   Name of the functor
 * @param[in] property   Name of the functor material property holding the change
 * @param[in] functor_is_ad   Whether the functor is AD
 */
void addNonlinearStepFunctorMaterial(FEProblemBase & problem,
                                     Factory & factory,
                                     const std::string & name,
                                     const std::vector<SubdomainName> & blocks,
                                     const std::string & functor_name,
                                     const std::string & property,
                                     bool functor_is_ad);

/**
 * Adds a postprocessor computing the maximum of a functor, normalized, on each nonlinear
 * convergence check; the counterpart of Component::addMaximumFunctorPostprocessor()
 */
void addMaximumFunctorPostprocessor(FEProblemBase & problem,
                                    Factory & factory,
                                    const std::string & functor_name,
                                    const std::string & pp_name,
                                    Real normalization,
                                    const std::vector<SubdomainName> & subdomains);

/**
 * Adds a convergence object that combines postprocessors, each with a tolerance; the counterpart of
 * Component::addMultiPostprocessorConvergence()
 */
void addMultiPostprocessorConvergence(FEProblemBase & problem,
                                      Factory & factory,
                                      const ConvergenceName & name,
                                      const std::vector<PostprocessorName> & postprocessors,
                                      const std::vector<std::string> & descriptions,
                                      const std::vector<Real> & tolerances);
}

/**
 * Interface for an ActionComponent to declare whether it has a nonlinear Convergence object, which
 * ComponentsConvergence combines with those of the other components. This is the counterpart of
 * Component::getNonlinearConvergence() for classic components.
 *
 * An ActionComponent that does not implement this interface is an error for ComponentsConvergence,
 * since whether its convergence matters is unknown. A component that has no convergence criterion
 * of its own derives from NoNonlinearConvergenceAC.
 */
class ActionComponentConvergenceInterface
{
public:
  virtual ~ActionComponentConvergenceInterface() = default;

  /// Name of the nonlinear Convergence object of the component, or nothing if it has none
  virtual std::optional<ConvergenceName> nonlinearConvergenceName() const = 0;
};

/**
 * For an ActionComponent that does not have a nonlinear Convergence object, such as a boundary
 * condition, which is then ignored by ComponentsConvergence.
 */
class NoNonlinearConvergenceAC : public ActionComponentConvergenceInterface
{
public:
  virtual std::optional<ConvergenceName> nonlinearConvergenceName() const override
  {
    return std::nullopt;
  }
};
