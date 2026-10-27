//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FlowChannel1PhasePhysics.h"

/**
 * Physics for a 1D flow channel ActionComponent (FlowChannelGasMixAC) with a binary gas mixture,
 * the counterpart of FlowModelGasMix. It adds the mass equation of the secondary gas (the variable
 * xi*rho*A, with xi the mass fraction of the secondary gas) to the equations of the single-phase
 * physics, with mass diffusion in the species and energy equations, and computes the fluid
 * properties of the mixture from them.
 *
 * Passive transport variables are not supported for a gas mixture.
 */
class FlowChannel1PhaseGasMixPhysics : public FlowChannel1PhasePhysics
{
public:
  static InputParameters validParams();
  FlowChannel1PhaseGasMixPhysics(const InputParameters & params);

  virtual std::vector<VariableName> solutionVariableNames() const override;

protected:
  virtual void addSolverVariables() override;
  virtual void addAuxiliaryVariables() override;
  virtual void addInitialConditions() override;
  virtual void addFEKernels() override;
  virtual void addDGKernels() override;
  virtual void addAuxiliaryKernels() override;

  virtual void checkFluidProperties() const override;
  virtual void addRhoEAInitialCondition() override;
  virtual void addDensityInitialCondition() override;
  virtual void addAdvectionDGKernels() override;
  virtual void addFluidPropertiesMaterial() override;
  virtual void addSlopeReconstructionMaterial() override;
  virtual void addNumericalFluxUserObject() override;
  virtual void addPressureTemperatureAuxKernels() override;
};
