//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatSourceBaseAC.h"
#include "HeatStructureInterfaceAC.h"

InputParameters
HeatSourceBaseAC::validParams()
{
  InputParameters params = ActionComponent::validParams();
  params.addRequiredParam<std::string>("hs", "Heat structure in which to apply heat source");
  params.addRequiredParam<std::vector<std::string>>(
      "regions", "Heat structure regions where heat generation is to be applied");
  return params;
}

HeatSourceBaseAC::HeatSourceBaseAC(const InputParameters & params) : ActionComponent(params)
{
  if (getParam<std::vector<std::string>>("regions").empty())
    paramError("regions", "At least one region must be given.");
}

const HeatStructureInterfaceAC &
HeatSourceBaseAC::hs() const
{
  return component<HeatStructureInterfaceAC>("hs", "HeatStructureInterface");
}

void
HeatSourceBaseAC::checkRegions() const
{
  for (const auto & region : getParam<std::vector<std::string>>("regions"))
    if (!hs().hasRegion(region))
      paramError("regions",
                 "Region '",
                 region,
                 "' does not exist in heat structure '",
                 getParam<std::string>("hs"),
                 "'.");
}

std::vector<SubdomainName>
HeatSourceBaseAC::subdomains() const
{
  std::vector<SubdomainName> blocks;
  for (const auto & region : getParam<std::vector<std::string>>("regions"))
    blocks.push_back(genName(getParam<std::string>("hs"), region));
  return blocks;
}
