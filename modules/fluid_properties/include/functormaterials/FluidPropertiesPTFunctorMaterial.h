//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorMaterial.h"

class SinglePhaseFluidProperties;

/**
 * Computes fluid properties from pressure and temperature functors.
 */
class FluidPropertiesPTFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  FluidPropertiesPTFunctorMaterial(const InputParameters & parameters);

protected:
  /// pressure
  const Moose::Functor<ADReal> & _pressure;
  /// temperature
  const Moose::Functor<ADReal> & _temperature;

  /// fluid properties user object
  const SinglePhaseFluidProperties & _fp;
};
