//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatStructure2DCouplerAC.h"
#include "MeshAlignment.h"

// Registered under the classic Component's name ('HeatStructure2DCoupler'); see
// FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructure2DCouplerAC,
                           "HeatStructure2DCoupler",
                           "add_bc");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HeatStructure2DCouplerAC,
                               "HeatStructure2DCoupler");

InputParameters
HeatStructure2DCouplerAC::validParams()
{
  InputParameters params = HeatStructure2DCouplerBaseAC::validParams();
  params.addRequiredParam<FunctionName>("heat_transfer_coefficient",
                                        "Heat transfer coefficient function [W/(m^2-K)]");
  params.addClassDescription(
      "Couples boundaries of two 2D heat structures via a heat transfer coefficient");
  return params;
}

HeatStructure2DCouplerAC::HeatStructure2DCouplerAC(const InputParameters & params)
  : ActionComponent(params), HeatStructure2DCouplerBaseAC(params)
{
}

void
HeatStructure2DCouplerAC::checkCoupling()
{
  if (!meshAlignment().meshesAreAligned())
    paramError("primary_boundary", "The primary and secondary boundaries must be aligned.");
}

void
HeatStructure2DCouplerAC::addBCs()
{
  for (const auto i : make_range(2u))
  {
    const std::string class_name = "HeatStructure2DCouplerBC";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<NonlinearVariableName>("variable") =
        hs(i).getPhysics().getParam<VariableName>("temperature_name");
    params.set<std::string>("coupled_variable") =
        hs(i == 0 ? 1 : 0).getPhysics().getParam<VariableName>("temperature_name");
    params.set<std::vector<BoundaryName>>("boundary") = {_hs_boundaries[i]};
    params.set<MeshAlignment *>("_mesh_alignment") = &meshAlignment();
    params.set<FunctionName>("heat_transfer_coefficient") =
        getParam<FunctionName>("heat_transfer_coefficient");
    params.set<Real>("coupling_area_fraction") = _coupling_area_fractions[i];
    getProblem().addBoundaryCondition(class_name, genName(name(), class_name, i), params);
  }
}
