//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HSBoundaryAC.h"
#include "THMUtils.h"
#include "HeatStructureBaseAC.h"
#include "HeatStructureInterfaceAC.h"
#include "HeatStructureCylindricalAC.h"
#include "THMVariableCoordinator.h"
#include "MooseUtils.h"

InputParameters
HSBoundaryAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addRequiredParam<std::vector<BoundaryName>>(
      "boundary", "List of boundary names for which this component applies");
  params.addRequiredParam<std::string>("hs", "Heat structure name");
  return params;
}

HSBoundaryAC::HSBoundaryAC(const InputParameters & params)
  : ActionComponent(params),
    _boundary(getParam<std::vector<BoundaryName>>("boundary")),
    _hs_name(getParam<std::string>("hs"))
{
  THM::informExperimental(*this);
  addRequiredTask("add_bc");
  addRequiredTask("check_integrity");
}

const HeatStructureInterfaceAC &
HSBoundaryAC::hs() const
{
  return _awh.getAction<HeatStructureInterfaceAC>(_hs_name);
}

const HeatStructureBaseAC *
HSBoundaryAC::hs2D() const
{
  return dynamic_cast<const HeatStructureBaseAC *>(&hs());
}

const HeatStructureCylindricalAC *
HSBoundaryAC::cylindricalHS() const
{
  return dynamic_cast<const HeatStructureCylindricalAC *>(&hs());
}

const VariableName &
HSBoundaryAC::temperatureName() const
{
  return hs().getPhysics().getParam<VariableName>("temperature_name");
}

void
HSBoundaryAC::requestAuxVariable(const VariableName & var_name,
                                 const std::string & family,
                                 const std::string & order)
{
  THMVariableCoordinator::findOrCreate(_awh, _action_factory)
      .requestVariable(false, var_name, family, order, 1.0, hs().blocks());
}

void
HSBoundaryAC::actOnAdditionalTasks()
{
  if (_current_task == "add_aux_variable")
    addAuxVariables();
  else if (_current_task == "add_function")
    addFunctions();
  else if (_current_task == "add_bc")
    addBCs();
  else if (_current_task == "add_postprocessor")
    addPostprocessors();
}

void
HSBoundaryAC::checkIntegrity()
{
  // Check that no perimeter is zero; if so, there is not physically a boundary. This only applies
  // to 2D heat structures, which have a perimeter.
  const auto * heat_structure_ptr = hs2D();
  if (!heat_structure_ptr)
    return;
  const auto & heat_structure = *heat_structure_ptr;
  for (const auto & boundary : _boundary)
    if (heat_structure.hasExternalBoundary(boundary))
    {
      const auto side = heat_structure.getExternalBoundaryType(boundary);
      if ((side == HeatStructureBaseAC::ExternalBoundaryType::INNER ||
           side == HeatStructureBaseAC::ExternalBoundaryType::OUTER) &&
          MooseUtils::absoluteFuzzyEqual(heat_structure.getUnitPerimeter(side), 0))
        paramError("boundary",
                   "The heat structure side of the heat structure '",
                   _hs_name,
                   "' corresponding to the boundary name '",
                   boundary,
                   "' has a zero perimeter. This can be caused by applying the boundary on the "
                   "axis of symmetry of a cylindrical heat structure.");
    }
}
