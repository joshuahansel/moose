//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Component.h"

/**
 * Dirac flux source for a FlowChannel1Phase.
 */
class DiracJunction1Phase : public Component
{
public:
  static InputParameters validParams();

  DiracJunction1Phase(const InputParameters & params);

  virtual void addMooseObjects() override;

protected:
  virtual void check() const override;

  /// Adds the DiracJunction1PhaseKernel for the specified variable
  void addDiracJunction1PhaseKernel(const std::string & var);

  /// Flow channel name
  const std::string & _flow_channel_name;
  /// Source point
  const Point & _point;
  /// DiracJunction1PhaseUserObject name
  const std::string _uo_name;
};
