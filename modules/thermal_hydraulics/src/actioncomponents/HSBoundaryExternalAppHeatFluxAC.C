//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HSBoundaryExternalAppHeatFluxAC.h"
#include "HeatStructureBaseAC.h"
#include "Conversion.h"

// Registered under the classic Component's name; see FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryExternalAppHeatFluxAC,
                           "HSBoundaryExternalAppHeatFlux",
                           "add_aux_variable");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryExternalAppHeatFluxAC,
                           "HSBoundaryExternalAppHeatFlux",
                           "add_function");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryExternalAppHeatFluxAC,
                           "HSBoundaryExternalAppHeatFlux",
                           "add_bc");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryExternalAppHeatFluxAC,
                           "HSBoundaryExternalAppHeatFlux",
                           "add_postprocessor");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HSBoundaryExternalAppHeatFluxAC,
                           "HSBoundaryExternalAppHeatFlux",
                           "check_integrity");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HSBoundaryExternalAppHeatFluxAC,
                               "HSBoundaryExternalAppHeatFlux");

InputParameters
HSBoundaryExternalAppHeatFluxAC::validParams()
{
  InputParameters params = HSBoundaryAC::validParams();
  params.addRequiredParam<VariableName>(
      "heat_flux_name",
      "Name to give the heat flux variable transferred from the external application");
  params.addParam<bool>(
      "heat_flux_is_monomial",
      true,
      "If true, makes the heat flux variable transferred from the external application to have the "
      "FE type 'CONSTANT MONOMIAL'. Else, the FE type is 'FIRST LAGRANGE'.");
  params.addRequiredParam<PostprocessorName>(
      "perimeter_ext", "Name to give the external application perimeter post-processor");
  params.addRequiredParam<bool>(
      "heat_flux_is_inward",
      "Set to true if the transferred heat flux corresponds to the inward direction on the heat "
      "structure boundary; else the outward direction");
  params.addClassDescription("Heat structure boundary condition to apply a heat flux transferred "
                             "from another application.");
  return params;
}

HSBoundaryExternalAppHeatFluxAC::HSBoundaryExternalAppHeatFluxAC(const InputParameters & params)
  : ActionComponent(params),
    HSBoundaryAC(params),
    _heat_flux_name(getParam<VariableName>("heat_flux_name")),
    _perimeter_ext_pp_name(getParam<PostprocessorName>("perimeter_ext")),
    _bc_coefficient_fn_name(genName(name(), "bc_coefficient_fn")),
    _heat_rate_prefactor_fn_name(genName(name(), "scale_fn"))
{
  addRequiredTask("add_aux_variable");
  addRequiredTask("add_function");
  addRequiredTask("add_postprocessor");
}

Real
HSBoundaryExternalAppHeatFluxAC::boundaryUnitPerimeter() const
{
  // HeatStructurePlate is not yet supported, matching classic HSBoundaryExternalAppHeatFlux
  if (!cylindricalHS())
    paramError("hs", "The heat structure '", _hs_name, "' must be a cylindrical heat structure.");

  const auto & heat_structure = *cylindricalHS();
  for (const auto & boundary : _boundary)
    if (!heat_structure.hasExternalBoundary(boundary))
      paramError("boundary", "The boundaries given in 'boundary' must all be external.");

  const auto side = heat_structure.getExternalBoundaryType(_boundary[0]);
  for (const auto & boundary : _boundary)
    if (heat_structure.getExternalBoundaryType(boundary) != side)
      paramError("boundary",
                 "The boundaries supplied in 'boundary' do not have a common external boundary "
                 "type.");
  if (side != HeatStructureBaseAC::ExternalBoundaryType::INNER &&
      side != HeatStructureBaseAC::ExternalBoundaryType::OUTER)
    paramError("boundary",
               "The boundaries in 'boundary' must be of an inner/outer type, not of a start/end "
               "type.");

  return heat_structure.getUnitPerimeter(side);
}

void
HSBoundaryExternalAppHeatFluxAC::addAuxVariables()
{
  if (getParam<bool>("heat_flux_is_monomial"))
    requestAuxVariable(_heat_flux_name, "MONOMIAL", "CONSTANT");
  else
    requestAuxVariable(_heat_flux_name, "LAGRANGE", "FIRST");
}

void
HSBoundaryExternalAppHeatFluxAC::addFunctions()
{
  // Scale by P_ext / P_hs rather than P_ext; see class documentation
  const std::string P_hs = Moose::stringifyExact(boundaryUnitPerimeter());
  const std::string sign = getParam<bool>("heat_flux_is_inward") ? "1" : "-1";

  const std::string class_name = "ParsedFunction";
  {
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::string>("expression") = "P / P_hs";
    params.set<std::vector<std::string>>("symbol_names") = {"P", "P_hs"};
    params.set<std::vector<std::string>>("symbol_values") = {_perimeter_ext_pp_name, P_hs};
    getProblem().addFunction(class_name, _bc_coefficient_fn_name, params);
  }
  {
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::string>("expression") = "sign * P / P_hs";
    params.set<std::vector<std::string>>("symbol_names") = {"sign", "P", "P_hs"};
    params.set<std::vector<std::string>>("symbol_values") = {sign, _perimeter_ext_pp_name, P_hs};
    getProblem().addFunction(class_name, _heat_rate_prefactor_fn_name, params);
  }
}

void
HSBoundaryExternalAppHeatFluxAC::addBCs()
{
  const std::string class_name = "FunctorNeumannBC";
  InputParameters params = getFactory().getValidParams(class_name);
  params.set<std::vector<BoundaryName>>("boundary") = _boundary;
  params.set<NonlinearVariableName>("variable") = temperatureName();
  params.set<MooseFunctorName>("functor") = _heat_flux_name;
  params.set<MooseFunctorName>("coefficient") = _bc_coefficient_fn_name;
  params.set<bool>("flux_is_inward") = getParam<bool>("heat_flux_is_inward");
  getProblem().addBoundaryCondition(class_name, genName(name(), "bc"), params);
}

void
HSBoundaryExternalAppHeatFluxAC::addPostprocessors()
{
  {
    const std::string class_name = "Receiver";
    InputParameters params = getFactory().getValidParams(class_name);
    getProblem().addPostprocessor(class_name, _perimeter_ext_pp_name, params);
  }
  {
    const std::string class_name = "SideIntegralFunctorPostprocessor";
    InputParameters params = getFactory().getValidParams(class_name);
    params.set<std::vector<BoundaryName>>("boundary") = _boundary;
    params.set<MooseFunctorName>("functor") = _heat_flux_name;
    params.set<MooseFunctorName>("prefactor") = _heat_rate_prefactor_fn_name;
    params.set<MooseEnum>("functor_argument") = "qp";
    params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_END};
    getProblem().addPostprocessor(class_name, genSafeName(name(), "integral"), params);
  }
}
