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
#include "HeatStructurePhysics.h"

/**
 * Interface for heat structure ActionComponents (HeatStructureBaseAC, HeatStructureFromFile3DAC),
 * mirroring the classic HeatStructureInterface. Components that apply boundary conditions or heat
 * sources to a heat structure only need what is declared here; the 2D-specific boundary and
 * perimeter queries are on HeatStructureBaseAC.
 *
 * This class creates and attaches the HeatStructurePhysics that defines the heat conduction
 * equation, the same way FlowChannel1PhaseAC creates and attaches its own FlowChannel1PhasePhysics
 * - so no separate [Physics] block is needed.
 */
class HeatStructureInterfaceAC : public virtual ActionComponent
{
public:
  static InputParameters validParams();
  HeatStructureInterfaceAC(const InputParameters & params);

  /// Returns this component's heat conduction physics
  const HeatStructurePhysics & getPhysics() const
  {
    mooseAssert(_physics, "The physics has not been created yet");
    return *_physics;
  }

  /// Whether \p region is one of this heat structure's regions
  virtual bool hasRegion(const std::string & region) const = 0;
  /// Number of units (rods) represented by this heat structure
  virtual Real numberOfUnits() const = 0;
  /// Length of the heat structure, used to normalize a uniform axial power shape
  virtual Real length() const = 0;

protected:
  /// Creates the heat conduction physics and registers it with the action warehouse
  void createPhysics(const InputParameters & params, bool second_order);

  virtual void addPhysics() override;

  /// The heat conduction physics this component creates and attaches to itself
  HeatStructurePhysics * _physics;
};
