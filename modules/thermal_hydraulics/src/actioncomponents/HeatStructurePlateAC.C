//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatStructurePlateAC.h"

// Registered under the classic Component's name ('HeatStructurePlate') rather than this class's
// own name; see FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructurePlateAC,
                           "HeatStructurePlate",
                           "add_mesh_generator");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructurePlateAC,
                           "HeatStructurePlate",
                           "init_component_physics");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructurePlateAC,
                           "HeatStructurePlate",
                           "add_material");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructurePlateAC,
                           "HeatStructurePlate",
                           "check_integrity");
registerActionComponentAliased("ThermalHydraulicsApp", HeatStructurePlateAC, "HeatStructurePlate");

InputParameters
HeatStructurePlateAC::validParams()
{
  InputParameters params = HeatStructureBaseAC::validParams();
  params.addRequiredParam<Real>("depth", "Dimension of the plate in the third direction [m]");
  params.addClassDescription("Plate heat structure.");
  return params;
}

HeatStructurePlateAC::HeatStructurePlateAC(const InputParameters & params)
  : ActionComponent(params), HeatStructureBaseAC(params), _depth(getParam<Real>("depth"))
{
}

Real
HeatStructurePlateAC::getUnitPerimeter(ExternalBoundaryType side) const
{
  switch (side)
  {
    case ExternalBoundaryType::INNER:
    case ExternalBoundaryType::OUTER:
      return _depth;
    case ExternalBoundaryType::START:
    case ExternalBoundaryType::END:
      return std::numeric_limits<Real>::quiet_NaN();
  }
  mooseError(name(), ": Unknown external boundary type.");
}

Real
HeatStructurePlateAC::computeRadialBoundaryArea(Real length, Real /*y*/) const
{
  return length * _depth;
}

Real
HeatStructurePlateAC::computeAxialBoundaryArea(Real y_min, Real y_max) const
{
  return (y_max - y_min) * _depth;
}
