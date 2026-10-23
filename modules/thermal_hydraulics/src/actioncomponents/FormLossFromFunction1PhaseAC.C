//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FormLossFromFunction1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"

// Registered under the classic Component's name ('FormLossFromFunction1Phase'); see
// FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           FormLossFromFunction1PhaseAC,
                           "FormLossFromFunction1Phase",
                           "add_aux_variable");
registerMooseActionAliased("ThermalHydraulicsApp",
                           FormLossFromFunction1PhaseAC,
                           "FormLossFromFunction1Phase",
                           "add_material");
registerMooseActionAliased("ThermalHydraulicsApp",
                           FormLossFromFunction1PhaseAC,
                           "FormLossFromFunction1Phase",
                           "add_kernel");
registerActionComponentAliased("ThermalHydraulicsApp",
                               FormLossFromFunction1PhaseAC,
                               "FormLossFromFunction1Phase");

InputParameters
FormLossFromFunction1PhaseAC::validParams()
{
  InputParameters params = FormLoss1PhaseBaseAC::validParams();
  params.addRequiredParam<FunctionName>("K_prime",
                                        "Form loss coefficient per unit length function [1/m]");
  params.addClassDescription(
      "Prescribe a form loss over a 1-phase flow channel given by a function");
  return params;
}

FormLossFromFunction1PhaseAC::FormLossFromFunction1PhaseAC(const InputParameters & params)
  : ActionComponent(params), FormLoss1PhaseBaseAC(params)
{
}

void
FormLossFromFunction1PhaseAC::addFormLossMaterials()
{
  const std::string class_name = "ADGenericFunctionMaterial";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<SubdomainName>>("block") = connectedPhysics().getClosuresBlocks();
  params.set<std::vector<std::string>>("prop_names") = {"K_prime"};
  params.set<std::vector<FunctionName>>("prop_values") = {getParam<FunctionName>("K_prime")};
  getProblem().addMaterial(class_name, genName(name(), "k_prime_material"), params);
}
