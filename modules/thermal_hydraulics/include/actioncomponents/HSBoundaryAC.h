//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ActionComponent.h"
#include "NamingInterface.h"

class HeatStructureInterfaceAC;
class HeatStructureBaseAC;
class HeatStructureCylindricalAC;

/**
 * Base class for a boundary condition ActionComponent applied to boundaries of a heat structure
 * ActionComponent (HeatStructureCylindricalAC/HeatStructurePlateAC), addressed the same way as
 * classic THM's HSBoundary: through 'hs' (the heat structure name) and 'boundary'.
 *
 * Classic THM cylindrical heat structures live on a Cartesian mesh and integrate their own RZ
 * weighting, so classic HSBoundary components pick RZ-specific boundary conditions and
 * postprocessors (e.g. ADConvectionHeatTransferRZBC, HeatRateConvectionRZ) that multiply by the
 * local circumference. HeatStructureCylindricalAC instead uses MOOSE's native COORD_RZ coordinate
 * system, which already applies that weighting during assembly, so the derived classes use the
 * plain (Cartesian) boundary conditions and postprocessors for both geometries.
 */
class HSBoundaryAC : public virtual ActionComponent, public NamingInterface
{
public:
  static InputParameters validParams();
  HSBoundaryAC(const InputParameters & params);

protected:
  virtual void actOnAdditionalTasks() override;
  virtual void checkIntegrity() override;

  /// Adds auxiliary variables (on the 'add_aux_variable' task)
  virtual void addAuxVariables() {}
  /// Adds functions (on the 'add_function' task)
  virtual void addFunctions() {}
  /// Adds the boundary condition(s) (on the 'add_bc' task)
  virtual void addBCs() = 0;
  /// Adds postprocessors (on the 'add_postprocessor' task)
  virtual void addPostprocessors() {}

  /// The heat structure component, looked up by name lazily so that 'hs' may name a component
  /// declared before or after this one
  const HeatStructureInterfaceAC & hs() const;
  /// The heat structure as a 2D heat structure, or nullptr if it is not (e.g. loaded from a file)
  const HeatStructureBaseAC * hs2D() const;
  /// The heat structure as a cylindrical heat structure, or nullptr if it is not cylindrical
  const HeatStructureCylindricalAC * cylindricalHS() const;
  /// Name of the heat structure's temperature variable
  const VariableName & temperatureName() const;

  /// Requests an auxiliary variable on the heat structure's blocks through THMVariableCoordinator,
  /// so that several boundary components may share one variable name (as classic THM's
  /// Simulation::addSimVariable() allows)
  void requestAuxVariable(const VariableName & var_name,
                          const std::string & family,
                          const std::string & order);

  /// Boundaries to which this component applies
  const std::vector<BoundaryName> _boundary;
  /// Name of the heat structure component
  const std::string _hs_name;
};
