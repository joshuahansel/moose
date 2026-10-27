//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ActionComponentConvergenceInterface.h"
#include "ActionComponent.h"
#include "NamingInterface.h"

class HeatStructureInterfaceAC;

/**
 * Base class for ActionComponents that apply a heat source in regions of a heat structure
 * ActionComponent, mirroring the classic HeatSourceBase.
 */
class HeatSourceBaseAC : public virtual ActionComponent,
                         public NamingInterface,
                         public NoNonlinearConvergenceAC
{
public:
  static InputParameters validParams();
  HeatSourceBaseAC(const InputParameters & params);

protected:
  /// The heat structure, looked up by name lazily so 'hs' may be declared after this component
  const HeatStructureInterfaceAC & hs() const;

  /// The component named by parameter \p param, which must be of type \p T (named \p type_name
  /// in the error message)
  template <typename T>
  const T & component(const std::string & param, const std::string & type_name) const;

  /// Errors if one of the regions does not exist in the heat structure
  void checkRegions() const;

  /// Names of the heat structure subdomains in the regions this source is applied to
  std::vector<SubdomainName> subdomains() const;
};

template <typename T>
const T &
HeatSourceBaseAC::component(const std::string & param, const std::string & type_name) const
{
  const auto & comp_name = getParam<std::string>(param);
  for (const auto & action : _awh.getActionListByName("list_component"))
    if (action->name() == comp_name)
    {
      if (const auto * comp = dynamic_cast<const T *>(action))
        return *comp;
      paramError(param, "The component '", comp_name, "' is not of type '", type_name, "'");
    }
  paramError(param, "The component '", comp_name, "' does not exist");
}
