//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

/**
 *
 */
template <bool is_ad>
class HeatConductionFluxMaterialTempl : public Material
{
public:
  static InputParameters validParams();

  HeatConductionFluxMaterialTempl(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Temperature gradient
  const GenericVariableGradient<is_ad> & _grad_T;
  /// Thermal conductivity
  const GenericMaterialProperty<Real, is_ad> & _k;
  /// Heat flux
  GenericMaterialProperty<RealVectorValue, is_ad> & _heat_flux;
};

typedef HeatConductionFluxMaterialTempl<false> HeatConductionFluxMaterial;
typedef HeatConductionFluxMaterialTempl<true> ADHeatConductionFluxMaterial;
