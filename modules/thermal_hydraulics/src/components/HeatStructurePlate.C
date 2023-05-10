//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatStructurePlate.h"

registerMooseObject("ThermalHydraulicsApp", HeatStructurePlate);

InputParameters
HeatStructurePlate::validParams()
{
  InputParameters params = HeatStructureBase::validParams();

  params.addRequiredParam<std::vector<std::string>>("names", "Name of each transverse region");
  params.addRequiredParam<std::vector<Real>>("widths", "Width of each transverse region [m]");
  params.addRequiredParam<std::vector<unsigned int>>(
      "n_part_elems", "Number of elements of each transverse region");
  params.addParam<std::vector<std::string>>("materials",
                                            "Material name for each transverse region");
  params.addParam<Real>("num_rods", 1.0, "Number of rods represented by this heat structure");
  params.addRequiredParam<Real>("depth", "Dimension of plate fuel in the third direction [m]");

  params.addClassDescription("Plate heat structure");

  return params;
}

HeatStructurePlate::HeatStructurePlate(const InputParameters & params)
  : HeatStructureBase(params), _depth(getParam<Real>("depth"))
{
  _material_names = getParam<std::vector<std::string>>("materials");
}

void
HeatStructurePlate::check() const
{
  HeatStructureBase::check();

  checkEqualSize<std::string, unsigned int>("names", "n_part_elems");
  checkEqualSize<std::string, Real>("names", "widths");
  if (isParamValid("materials"))
    checkEqualSize<std::string, std::string>("names", "materials");
}

std::vector<std::string>
HeatStructurePlate::getRegionNames() const
{
  return getParam<std::vector<std::string>>("names");
}

std::vector<Real>
HeatStructurePlate::getRegionWidths() const
{
  return getParam<std::vector<Real>>("widths");
}

std::vector<unsigned int>
HeatStructurePlate::getRegionNumbersOfElements() const
{
  return getParam<std::vector<unsigned int>>("n_part_elems");
}

Real
HeatStructurePlate::getNumberOfUnits() const
{
  return getParam<Real>("num_rods");
}
