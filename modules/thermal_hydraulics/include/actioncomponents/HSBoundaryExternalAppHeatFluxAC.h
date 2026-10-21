//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HSBoundaryAC.h"

/**
 * Cylindrical heat structure ActionComponent boundary condition applying a heat flux transferred
 * from an external application, along with the external application's heated perimeter
 * (transferred into the 'perimeter_ext' post-processor this component creates).
 *
 * The transferred heat flux is relative to the external perimeter, so the heat rate per unit
 * length entering the heat structure is q_ext * P_ext. Classic HSBoundaryExternalAppHeatFlux
 * applies q_ext * P_ext directly, since classic cylindrical heat structures do not weight their
 * boundary integrals by the circumference. The native COORD_RZ coordinate system of
 * HeatStructureCylindricalAC does weight them by the heat structure's own boundary perimeter P_hs,
 * so the boundary condition (and heat rate post-processor) here are scaled by P_ext / P_hs instead.
 */
class HSBoundaryExternalAppHeatFluxAC : public HSBoundaryAC
{
public:
  static InputParameters validParams();
  HSBoundaryExternalAppHeatFluxAC(const InputParameters & params);

protected:
  virtual void addAuxVariables() override;
  virtual void addFunctions() override;
  virtual void addBCs() override;
  virtual void addPostprocessors() override;

  /// Checks that the heat structure is cylindrical and that 'boundary' consists only of inner or
  /// only of outer boundaries, then returns the heat structure's unit perimeter on that side
  Real boundaryUnitPerimeter() const;

  /// Name of the heat flux variable transferred from the external application
  const VariableName _heat_flux_name;
  /// Name of the external application perimeter post-processor
  const PostprocessorName _perimeter_ext_pp_name;
  /// Name of the function scaling the transferred heat flux in the boundary condition
  const FunctionName _bc_coefficient_fn_name;
  /// Name of the function scaling the transferred heat flux in the heat rate post-processor
  const FunctionName _heat_rate_prefactor_fn_name;
};
