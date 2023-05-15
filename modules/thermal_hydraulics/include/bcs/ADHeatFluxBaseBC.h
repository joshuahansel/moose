//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADIntegratedBC.h"

class ADHeatFluxFromHeatStructureBaseUserObject;

/**
 * Base class for handling heat flux between flow channels and heat structures
 */
class ADHeatFluxBaseBC : public ADIntegratedBC
{
public:
  ADHeatFluxBaseBC(const InputParameters & parameters);

protected:
  /// User object that computes the heat flux
  const ADHeatFluxFromHeatStructureBaseUserObject & _q_uo;
  /// Perimeter of a single unit of heat structure
  const Real _P_hs_unit;
  /// Number of units of heat structure
  const unsigned int _n_unit;
  /// Factor by which to scale term on the flow channel side for the heat structure side
  const Real _hs_scale;

public:
  static InputParameters validParams();
};
