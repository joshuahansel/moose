//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatStructure2DRadiationCouplerRZAC.h"
#include "HeatConductionNames.h"
#include "MeshAlignment.h"

// Registered under the classic Component's name ('HeatStructure2DRadiationCouplerRZ'); see
// FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructure2DRadiationCouplerRZAC,
                           "HeatStructure2DRadiationCouplerRZ",
                           "add_bc");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HeatStructure2DRadiationCouplerRZAC,
                               "HeatStructure2DRadiationCouplerRZ");

InputParameters
HeatStructure2DRadiationCouplerRZAC::validParams()
{
  InputParameters params = HeatStructure2DCouplerBaseAC::validParams();
  params.addRequiredParam<Real>("primary_emissivity", "Emissivity for the primary side");
  params.addRequiredParam<Real>("secondary_emissivity", "Emissivity for the secondary side");
  params.addParam<Real>("stefan_boltzmann_constant",
                        HeatConduction::Constants::sigma,
                        "Stefan Boltzmann constant [W/(m^2-K^4)]. This constant is provided as a "
                        "parameter to allow different precisions.");
  params.addClassDescription(
      "Couples boundaries of two 2D cylindrical heat structures via radiation");
  return params;
}

HeatStructure2DRadiationCouplerRZAC::HeatStructure2DRadiationCouplerRZAC(
    const InputParameters & params)
  : ActionComponent(params),
    HeatStructure2DCouplerBaseAC(params),
    _emissivities({getParam<Real>("primary_emissivity"), getParam<Real>("secondary_emissivity")})
{
}

void
HeatStructure2DRadiationCouplerRZAC::checkCoupling()
{
  if (!isCylindrical(0) || !isCylindrical(1))
    paramError("primary_heat_structure",
               "The primary and secondary heat structures must be of a type inherited from "
               "'HeatStructureCylindricalBase'.");

  if (!meshAlignment().meshesAreAligned())
    paramError("primary_boundary", "The primary and secondary boundaries must be aligned.");

  for (const auto i : make_range(2u))
  {
    const auto side = hs(i).getExternalBoundaryType(_hs_boundaries[i]);
    if (side == HeatStructureBaseAC::ExternalBoundaryType::START ||
        side == HeatStructureBaseAC::ExternalBoundaryType::END)
      paramError("primary_boundary",
                 "The primary and secondary boundaries must be radial boundaries.");
  }
}

void
HeatStructure2DRadiationCouplerRZAC::addBCs()
{
  for (const auto i : make_range(2u))
  {
    const unsigned int j = i == 0 ? 1 : 0;

    // The view factor of the smaller of the two boundaries is the ratio of the areas
    const Real view_factor = _coupling_area_fractions[i];

    const std::string class_name = "HeatStructure2DRadiationCouplerBC";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<NonlinearVariableName>("variable") =
        hs(i).getPhysics().getParam<VariableName>("temperature_name");
    params.set<std::string>("coupled_variable") =
        hs(j).getPhysics().getParam<VariableName>("temperature_name");
    params.set<std::vector<BoundaryName>>("boundary") = {_hs_boundaries[i]};
    params.set<MeshAlignment *>("_mesh_alignment") = &meshAlignment();
    params.set<Real>("emissivity") = _emissivities[i];
    params.set<Real>("coupled_emissivity") = _emissivities[j];
    params.set<Real>("view_factor") = view_factor;
    params.set<Real>("area") = _areas[i];
    params.set<Real>("coupled_area") = _areas[j];
    params.set<Real>("stefan_boltzmann_constant") = getParam<Real>("stefan_boltzmann_constant");
    getProblem().addBoundaryCondition(class_name, genName(name(), class_name, i), params);
  }
}
