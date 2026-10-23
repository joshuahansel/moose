//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FormLossFromExternalApp1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"
#include "THMVariableCoordinator.h"

// Registered under the classic Component's name ('FormLossFromExternalApp1Phase'); see
// FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           FormLossFromExternalApp1PhaseAC,
                           "FormLossFromExternalApp1Phase",
                           "add_aux_variable");
registerMooseActionAliased("ThermalHydraulicsApp",
                           FormLossFromExternalApp1PhaseAC,
                           "FormLossFromExternalApp1Phase",
                           "add_material");
registerMooseActionAliased("ThermalHydraulicsApp",
                           FormLossFromExternalApp1PhaseAC,
                           "FormLossFromExternalApp1Phase",
                           "add_kernel");
registerActionComponentAliased("ThermalHydraulicsApp",
                               FormLossFromExternalApp1PhaseAC,
                               "FormLossFromExternalApp1Phase");

InputParameters
FormLossFromExternalApp1PhaseAC::validParams()
{
  InputParameters params = FormLoss1PhaseBaseAC::validParams();
  params.addClassDescription("Apply a distributed form loss over a 1-phase flow channel computed "
                             "by an external application.");
  return params;
}

FormLossFromExternalApp1PhaseAC::FormLossFromExternalApp1PhaseAC(const InputParameters & params)
  : ActionComponent(params), FormLoss1PhaseBaseAC(params), _K_prime_var_name("K_prime")
{
}

void
FormLossFromExternalApp1PhaseAC::addFormLossVariables()
{
  THMVariableCoordinator::findOrCreate(_awh, _action_factory)
      .requestVariable(false,
                       _K_prime_var_name,
                       "LAGRANGE",
                       "FIRST",
                       1.0,
                       connectedPhysics().getClosuresBlocks());
}

void
FormLossFromExternalApp1PhaseAC::addFormLossMaterials()
{
  const std::string class_name = "ADCoupledVariableValueMaterial";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
  params.set<MaterialPropertyName>("prop_name") = "K_prime";
  params.set<std::vector<VariableName>>("coupled_variable") = {_K_prime_var_name};
  getProblem().addMaterial(class_name, genName(name(), "k_prime_material"), params);
}
