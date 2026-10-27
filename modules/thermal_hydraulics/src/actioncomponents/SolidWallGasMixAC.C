//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolidWallGasMixAC.h"
#include "THMNames.h"
#include "THMUtils.h"
#include "MooseUtils.h"

// Registered under the classic Component's name ('SolidWallGasMix'); see FlowChannel1PhaseAC.C for
// why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           SolidWallGasMixAC,
                           "SolidWallGasMix",
                           "add_user_object");
registerMooseActionAliased("ThermalHydraulicsApp", SolidWallGasMixAC, "SolidWallGasMix", "add_bc");
registerActionComponentAliased("ThermalHydraulicsApp", SolidWallGasMixAC, "SolidWallGasMix");

InputParameters
SolidWallGasMixAC::validParams()
{
  InputParameters params = FlowBoundary1PhaseAC::validParams();
  params.addClassDescription("Wall boundary condition component for FlowChannelGasMixAC.");
  return params;
}

SolidWallGasMixAC::SolidWallGasMixAC(const InputParameters & params)
  : ActionComponent(params), FlowBoundary1PhaseAC(params)
{
}

void
SolidWallGasMixAC::addUserObjects()
{
  const std::string class_name = "BoundaryFluxGasMixGhostWall";
  InputParameters params = _factory.getValidParams(class_name);
  params.set<UserObjectName>("numerical_flux") = connectedPhysics().numericalFluxUserObjectName();
  params.set<Real>("normal") = _normal;
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};
  getProblem().addUserObject(class_name, _boundary_uo_name, params);
}

void
SolidWallGasMixAC::addWeakBCs()
{
  const std::string class_name = "BoundaryFluxGasMixBC";
  InputParameters params = _factory.getValidParams(class_name);
  params.set<std::vector<BoundaryName>>("boundary") = {_input};
  params.set<Real>("normal") = _normal;
  params.set<UserObjectName>("boundary_flux") = _boundary_uo_name;
  params.set<std::vector<VariableName>>("A_linear") = {THM::AREA_LINEAR};
  params.set<std::vector<VariableName>>("xirhoA") = {THM::XIRHOA};
  params.set<std::vector<VariableName>>("rhoA") = {THM::RHOA};
  params.set<std::vector<VariableName>>("rhouA") = {THM::RHOUA};
  params.set<std::vector<VariableName>>("rhoEA") = {THM::RHOEA};
  params.set<bool>("implicit") = THM::implicitTimeIntegrationFlag(_app);

  for (const auto & var : connectedPhysics().solutionVariableNames())
  {
    params.set<NonlinearVariableName>("variable") = var;
    getProblem().addBoundaryCondition(class_name, genName(name(), var, "bc"), params);
  }
}
