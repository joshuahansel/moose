//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "MooseTypes.h"

class PenetrationLocator;
class NearestNodeLocator;
class SystemBase;
namespace libMesh
{
template <typename>
class NumericVector;
}

/**
 * Base class for Creating an AD material property for a variable transferred
 * from the boundary of a 2D mesh onto a 1D mesh.
 */
class VariableTransferMaterialBase : public Material
{
public:
  static InputParameters validParams();

  VariableTransferMaterialBase(const InputParameters & parameters);

protected:
  /**
   * Gets the nodal variable values from the primary boundary on the current element
   */
  std::vector<ADReal> getNodalVariableValues() const;

  /// Penetration locator
  PenetrationLocator & _penetration_locator;
  /// Nearest node locator
  NearestNodeLocator & _nearest_node;
  /// Nonlinear system
  const SystemBase & _nl_sys;
  /// Solution vector
  const NumericVector<Number> * const & _serialized_solution;
  /// Variable number of the variable to transfer
  unsigned int _paired_variable;
};
