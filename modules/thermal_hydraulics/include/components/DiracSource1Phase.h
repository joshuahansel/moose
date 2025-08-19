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
 * Dirac mass and energy sources for a FlowChannel1Phase.
 */
class DiracSource1Phase : public Component
{
public:
  static InputParameters validParams();

  DiracSource1Phase(const InputParameters & params);

  virtual void addMooseObjects() override;

protected:
  virtual void check() const override;

  /// Adds a FunctorDiracKernel for the given variable and functor
  void addFunctorDiracKernel(const std::string & var, const MooseFunctorName & source_functor);

  /// Flow channel name
  const std::string & _flow_channel_name;
  /// Source point
  const Point & _point;
};
