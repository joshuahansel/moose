//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatStructureCylindrical.h"

registerMooseObject("ThermalHydraulicsApp", HeatStructureCylindrical);

InputParameters
HeatStructureCylindrical::validParams()
{
  InputParameters params = HeatStructureCylindricalBase::validParams();

  params.addRequiredParam<std::vector<std::string>>("names", "Name of each radial region");
  params.addRequiredParam<std::vector<Real>>("widths", "Width of each radial region [m]");
  params.addRequiredParam<std::vector<unsigned int>>("n_part_elems",
                                                     "Number of elements of each radial region");
  params.addParam<std::vector<std::string>>("materials", "Material name for each radial region");
  params.addParam<Real>("num_rods", 1.0, "Number of rods represented by this heat structure");
  params.addParam<Real>("inner_radius", 0., "Inner radius of the heat structure [m]");

  params.addClassDescription("Cylindrical heat structure");

  return params;
}

HeatStructureCylindrical::HeatStructureCylindrical(const InputParameters & params)
  : HeatStructureCylindricalBase(params)
{
  _material_names = getParam<std::vector<std::string>>("materials");
}

void
HeatStructureCylindrical::check() const
{
  HeatStructureCylindricalBase::check();

  checkEqualSize<std::string, unsigned int>("names", "n_part_elems");
  checkEqualSize<std::string, Real>("names", "widths");
  if (isParamValid("materials"))
    checkEqualSize<std::string, std::string>("names", "materials");
}

std::vector<std::string>
HeatStructureCylindrical::getRegionNames() const
{
  return getParam<std::vector<std::string>>("names");
}

std::vector<Real>
HeatStructureCylindrical::getRegionWidths() const
{
  return getParam<std::vector<Real>>("widths");
}

std::vector<unsigned int>
HeatStructureCylindrical::getRegionNumbersOfElements() const
{
  return getParam<std::vector<unsigned int>>("n_part_elems");
}

Real
HeatStructureCylindrical::getNumberOfUnits() const
{
  return getParam<Real>("num_rods");
}

Real
HeatStructureCylindrical::getInnerRadius() const
{
  return getParam<Real>("inner_radius");
}
