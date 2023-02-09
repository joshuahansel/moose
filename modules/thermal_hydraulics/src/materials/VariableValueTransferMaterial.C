//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableValueTransferMaterial.h"
#include "Assembly.h"

#include "libmesh/elem.h"

registerMooseObject("ThermalHydraulicsApp", VariableValueTransferMaterial);

InputParameters
VariableValueTransferMaterial::validParams()
{
  InputParameters params = VariableTransferMaterialBase::validParams();
  params.addRequiredParam<MaterialPropertyName>(
      "property_name",
      "The name of the material property that will be "
      "declared that will represent the transferred variable.");
  params.addClassDescription("Creates an AD material property for a variable transferred from the "
                             "boundary of a 2D mesh onto a 1D mesh.");
  return params;
}

VariableValueTransferMaterial::VariableValueTransferMaterial(const InputParameters & parameters)
  : VariableTransferMaterialBase(parameters),
    _phi(_assembly.fePhi<Real>(FEType(FIRST, LAGRANGE))),
    _prop(declareADProperty<Real>(getParam<MaterialPropertyName>("property_name")))
{
}

void
VariableValueTransferMaterial::computeProperties()
{
  const auto nodal_values = getNodalVariableValues();

  for (const auto qp : make_range(_qrule->n_points()))
  {
    _prop[qp] = 0;
    for (const auto i : _current_elem->node_index_range())
      _prop[qp] += nodal_values[i] * _phi[i][qp];
  }
}
