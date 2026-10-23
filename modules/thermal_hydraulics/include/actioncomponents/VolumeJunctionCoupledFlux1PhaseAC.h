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

class VolumeJunction1PhaseAC;

/**
 * Applies a flux between a VolumeJunction1PhaseAC and an external application, the
 * native-ActionComponent counterpart of VolumeJunctionCoupledFlux1Phase.
 *
 * The mass and energy flow rates to the junction are computed from the pressure and temperature of
 * the coupled volume, which come from post-processors (receivers) that the external application
 * writes into, and are applied to the junction equations by a kernel on each of its variables. If a
 * MultiApp is given, the flow rates are transferred to it and the pressure and temperature are
 * transferred back.
 */
class VolumeJunctionCoupledFlux1PhaseAC : public virtual ActionComponent, public NamingInterface
{
public:
  static InputParameters validParams();
  VolumeJunctionCoupledFlux1PhaseAC(const InputParameters & params);

protected:
  virtual void actOnAdditionalTasks() override;

  /// The volume junction, looked up by name lazily
  const VolumeJunction1PhaseAC & volumeJunction() const;
  /// Adds the kernel of the junction variable \p var, whose equation has index \p i
  void addVolumeJunctionKernel(const std::string & var, unsigned int i);
  /// Adds the post-processor computing the flow rate of an equation
  void addFluxPostprocessor(const std::string & equation);
  /// Adds the receiver post-processor for a property of the coupled volume
  void addPropertyPostprocessor(const std::string & property);
  /// Adds the transfer of a flow rate to the MultiApp
  void addFluxTransfer(const std::string & equation);
  /// Adds the transfer of a property from the MultiApp
  void addPropertyTransfer(const std::string & property);
  /// Appends the suffix of this component to a post-processor name
  PostprocessorName addPostprocessorSuffix(const std::string & base_name) const;

  /// Name of the volume junction
  const std::string _volume_junction_name;
  /// Direction pointing from the junction to the coupled volume, normalized
  const RealVectorValue _normal_from_junction;
  /// Suffix appended to the post-processor names
  const std::string _pp_suffix;
};
