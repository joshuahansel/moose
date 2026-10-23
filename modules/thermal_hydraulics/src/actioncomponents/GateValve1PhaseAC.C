//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GateValve1PhaseAC.h"
#include "FlowChannel1PhaseAC.h"
#include "FlowChannel1PhasePhysics.h"
#include "THMNames.h"
#include "THMUtils.h"
#include "MooseApp.h"
#include "Factory.h"

// Registered under the classic Component's name ('GateValve1Phase'); see FlowChannel1PhaseAC.C for
// why. The tasks are the same as for JunctionOneToOne1PhaseAC.
registerMooseActionAliased("ThermalHydraulicsApp",
                           GateValve1PhaseAC,
                           "GateValve1Phase",
                           "add_user_object");
registerMooseActionAliased("ThermalHydraulicsApp", GateValve1PhaseAC, "GateValve1Phase", "add_bc");
registerActionComponentAliased("ThermalHydraulicsApp", GateValve1PhaseAC, "GateValve1Phase");

InputParameters
GateValve1PhaseAC::validParams()
{
  InputParameters params = JunctionOneToOne1PhaseAC::validParams();
  params.addRequiredParam<Real>("open_area_fraction", "Fraction of flow area that is open [-]");
  params.declareControllable("open_area_fraction");
  params.addClassDescription("Gate valve component for 1-phase flow");
  return params;
}

GateValve1PhaseAC::GateValve1PhaseAC(const InputParameters & params)
  : ActionComponent(params), JunctionOneToOne1PhaseAC(params)
{
  // Log a warning if slope reconstruction is used on one or more of the adjacent flow channels
  bool slope_reconstruction_used = false;
  for (const auto & fch_name : _connected_component_names)
    for (const auto & action : _awh.getActionListByName("list_component"))
      if (action->name() == fch_name)
        if (const auto * fch = dynamic_cast<const FlowChannel1PhaseAC *>(action))
          slope_reconstruction_used =
              slope_reconstruction_used ||
              (fch->getParam<MooseEnum>("rdg_slope_reconstruction") != "none");
  if (slope_reconstruction_used)
    mooseWarning(name(),
                 ": Currently GateValve1Phase cannot perform slope reconstruction across the "
                 "junction, so the slopes on the adjacent elements will be zero.");
}

void
GateValve1PhaseAC::buildJunctionUserObject()
{
  ExecFlagEnum execute_on(MooseUtils::getDefaultExecFlagEnum());
  execute_on = {EXEC_INITIAL, EXEC_LINEAR, EXEC_NONLINEAR};

  // It is assumed that each channel has the same numerical flux, so the first one's is used
  const std::string class_name = "ADGateValve1PhaseUserObject";
  InputParameters params = _factory.getValidParams(class_name);
  params.set<std::vector<BoundaryName>>("boundary") = _connections;
  params.set<std::vector<Real>>("normals") = _normals;
  params.set<UserObjectName>("numerical_flux") = connectedPhysics(0).numericalFluxUserObjectName();
  params.set<Real>("open_area_fraction") = getParam<Real>("open_area_fraction");
  params.set<std::vector<VariableName>>("A") = {THM::AREA};
  params.set<std::vector<VariableName>>("rhoA") = {THM::RHOA};
  params.set<std::vector<VariableName>>("rhouA") = {THM::RHOUA};
  params.set<std::vector<VariableName>>("rhoEA") = {THM::RHOEA};
  params.set<std::string>("component_name") = name();
  params.set<ExecFlagEnum>("execute_on") = execute_on;
  getProblem().addUserObject(class_name, _junction_uo_name, params);
  THM::connectComponentControllableParam(
      *this, "open_area_fraction", class_name, _junction_uo_name);
}

void
GateValve1PhaseAC::addWeakBC(unsigned int i, const NonlinearVariableName & var)
{
  const std::string class_name = "ADGateValve1PhaseBC";
  InputParameters params = _factory.getValidParams(class_name);
  params.set<std::vector<BoundaryName>>("boundary") = {_connections[i]};
  params.set<Real>("normal") = _normals[i];
  params.set<NonlinearVariableName>("variable") = var;
  params.set<UserObjectName>("gate_valve_uo") = _junction_uo_name;
  params.set<unsigned int>("connection_index") = i;
  params.set<std::vector<VariableName>>("A_elem") = {THM::AREA};
  params.set<std::vector<VariableName>>("A_linear") = {THM::AREA_LINEAR};
  params.set<std::vector<VariableName>>("rhoA") = {THM::RHOA};
  params.set<std::vector<VariableName>>("rhouA") = {THM::RHOUA};
  params.set<std::vector<VariableName>>("rhoEA") = {THM::RHOEA};
  params.set<bool>("implicit") = THM::implicitTimeIntegrationFlag(_app);
  getProblem().addBoundaryCondition(class_name, genName(name(), i, var + ":" + class_name), params);
}
