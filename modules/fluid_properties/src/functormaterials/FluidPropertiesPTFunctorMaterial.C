//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FluidPropertiesPTFunctorMaterial.h"
#include "SinglePhaseFluidProperties.h"

registerMooseObject("FluidPropertiesApp", FluidPropertiesPTFunctorMaterial);

InputParameters
FluidPropertiesPTFunctorMaterial::validParams()
{
  auto params = FunctorMaterial::validParams();

  params.addClassDescription("Computes fluid properties from pressure and temperature functors.");

  // std::vector<MooseEnum> fluid_properties{MooseEnum("rho e")};
  params.addRequiredParam<std::vector<std::string>>(
      "fluid_properties", "Fluid properties for which to create functor material properties");
  params.addRequiredParam<std::vector<std::string>>(
      "property_names", "Functor material property name to use for each selected fluid property");
  params.addRequiredParam<MooseFunctorName>("temperature", "temperature functor");
  params.addRequiredParam<MooseFunctorName>("pressure", "pressure functor");
  params.addRequiredParam<UserObjectName>("fluid_properties_object",
                                          "SinglePhaseFluidProperties object");

  return params;
}

FluidPropertiesPTFunctorMaterial::FluidPropertiesPTFunctorMaterial(
    const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _pressure(getFunctor<ADReal>("pressure")),
    _temperature(getFunctor<ADReal>("temperature")),
    _fp(getUserObject<SinglePhaseFluidProperties>("fluid_properties_object"))
{
  const auto & fluid_properties = getParam<std::vector<std::string>>("fluid_properties");
  const auto & property_names = getParam<std::vector<std::string>>("property_names");

  if (fluid_properties.size() != property_names.size())
    mooseError("The parameters 'fluid_properties' and 'property_names' must have the same size.");

  for (const auto i : index_range(fluid_properties))
  {
    std::function<ADReal(const ADReal &, const ADReal &)> property_function;
    if (fluid_properties[i] == "rho")
      property_function = [this](const ADReal & p, const ADReal & T)
      { return _fp.rho_from_p_T(p, T); };
    else if (fluid_properties[i] == "e")
      property_function = [this](const ADReal & p, const ADReal & T)
      { return _fp.e_from_p_T(p, T); };
    else
      mooseError("Invalid fluid property in 'fluid_properties'.");

    addFunctorProperty<ADReal>(property_names[i],
                               [this, property_function](const auto & r, const auto & t) -> ADReal
                               { return property_function(_pressure(r, t), _temperature(r, t)); });
  }
}
