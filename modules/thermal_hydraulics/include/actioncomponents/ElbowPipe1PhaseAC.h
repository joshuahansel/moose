//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FlowChannel1PhaseAC.h"

/**
 * Bent pipe for 1-phase flow, the native-ActionComponent counterpart of ElbowPipe1Phase.
 *
 * The straight mesh of a FlowChannel1PhaseAC is generated along the x axis with the length of the
 * arc, and is then mapped onto the arc of a circle of the given radius in the xy plane, starting at
 * 'start_angle' and ending at 'end_angle'. The element sizes along the arc stay as in the straight
 * mesh. The rotation and the translation of the flow channel are applied after the mapping.
 */
class ElbowPipe1PhaseAC : public FlowChannel1PhaseAC
{
public:
  static InputParameters validParams();
  ElbowPipe1PhaseAC(const InputParameters & params);

protected:
  virtual void addShapeMeshGenerators() override;

  /// Radius of the pipe
  const Real _radius;
  /// Angle at which the pipe starts [degrees]
  const Real _start_angle;
  /// Angle at which the pipe ends [degrees]
  Real _end_angle;
  /// Central angle covered by the pipe [degrees]
  Real _central_angle;
};
