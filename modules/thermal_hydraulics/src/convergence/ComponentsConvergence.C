//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComponentsConvergence.h"
#include "THMProblem.h"
#include "Component.h"
#include "ActionComponent.h"
#include "ActionComponentConvergenceInterface.h"
#include "ActionWarehouse.h"
#include "MooseApp.h"
#include "FEProblemBase.h"

registerMooseObject("ThermalHydraulicsApp", ComponentsConvergence);

InputParameters
ComponentsConvergence::validParams()
{
  InputParameters params = IterationCountConvergence::validParams();

  params.addClassDescription(
      "Assesses convergence of all Component and ActionComponent objects in a simulation.");

  return params;
}

ComponentsConvergence::ComponentsConvergence(const InputParameters & parameters)
  : IterationCountConvergence(parameters),
    _fe_problem(*getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")),
    _thm_problem(dynamic_cast<THMProblem *>(&_fe_problem))
{
}

void
ComponentsConvergence::initialSetup()
{
  IterationCountConvergence::initialSetup();

  if (_thm_problem)
    for (const auto & comp : _thm_problem->getComponents())
    {
      const auto nl_conv = comp->getNonlinearConvergence();
      if (nl_conv)
        _convergence_objects.push_back(nl_conv);
    }

  for (const auto * component : _app.actionWarehouse().getActions<ActionComponent>())
  {
    const auto * interface = dynamic_cast<const ActionComponentConvergenceInterface *>(component);
    if (!interface)
      mooseError("The component '",
                 component->name(),
                 "' does not provide a nonlinear convergence, which ComponentsConvergence "
                 "requires of all ActionComponents.");

    if (const auto name = interface->nonlinearConvergenceName())
      _convergence_objects.push_back(&_fe_problem.getConvergence(*name));
  }
}

void
ComponentsConvergence::preLoop()
{
  for (auto & conv : _convergence_objects)
    conv->preLoop();
}

void
ComponentsConvergence::preIteration()
{
  for (auto & conv : _convergence_objects)
    conv->preIteration();
}

Convergence::MooseConvergenceStatus
ComponentsConvergence::checkConvergenceInner(unsigned int iter)
{
  bool all_converged = true;
  for (auto & conv : _convergence_objects)
  {
    const auto status = conv->checkConvergence(iter);
    if (status == Convergence::MooseConvergenceStatus::DIVERGED)
      return Convergence::MooseConvergenceStatus::DIVERGED;
    else if (status == Convergence::MooseConvergenceStatus::ITERATING)
      all_converged = false;
    else
    {
      // checking in case additional status added in future
      mooseAssert(status == Convergence::MooseConvergenceStatus::CONVERGED,
                  "Status not implemented");
    }
  }

  if (all_converged)
    return Convergence::MooseConvergenceStatus::CONVERGED;
  else
    return Convergence::MooseConvergenceStatus::ITERATING;
}
