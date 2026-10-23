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

class FlowChannel1PhasePhysics;

/**
 * Base class for prescribing a form loss over a FlowChannel1PhaseAC, the native-ActionComponent
 * counterpart of FormLoss1PhaseBase. It adds the momentum form loss kernel; derived classes provide
 * the form loss coefficient per unit length, as the material property 'K_prime'.
 */
class FormLoss1PhaseBaseAC : public virtual ActionComponent, public NamingInterface
{
public:
  static InputParameters validParams();
  FormLoss1PhaseBaseAC(const InputParameters & params);

protected:
  virtual void actOnAdditionalTasks() override;
  virtual void addMaterials() override;

  /// Adds the form loss coefficient per unit length (on the 'add_material' task)
  virtual void addFormLossMaterials() {}
  /// Adds variables of the form loss coefficient (on the 'add_aux_variable' task)
  virtual void addFormLossVariables() {}

  /// Physics of the flow channel, looked up by name lazily so that 'flow_channel' may name a
  /// component declared before or after this one
  const FlowChannel1PhasePhysics & connectedPhysics() const;
};
