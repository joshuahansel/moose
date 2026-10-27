//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PhysicsBase.h"
#include "GravityInterface.h"
#include "NamingInterface.h"
#include "THMVariableCoordinator.h"
#include "FlowChannelClosuresInterface.h"
#include "THMNames.h"

#include <optional>

class ClosuresBase;

/**
 * Single-phase flow (mass, momentum, energy) physics for a 1D flow channel ActionComponent (see
 * FlowChannel1PhaseAC), which creates and attaches one of these to itself in its constructor.
 *
 * This implements FlowChannelClosuresInterface so that the same ClosuresBase-derived objects
 * (Closures1PhaseSimple, Closures1PhaseTHM, ...) that classic THM's FlowChannelBase uses for wall
 * friction and wall heat transfer also work here - see FlowChannelClosuresInterface for why.
 *
 * A heat-transfer-to-flow-channel ActionComponent (HeatTransferFromHeatFlux1PhaseAC,
 * HeatTransferFromSpecifiedTemperature1PhaseAC) registers itself here via addHeatTransferName(),
 * the same way classic HeatTransferBase registers with FlowChannelBase. With a single connection
 * the bare THM:: wall-quantity names are used; with several, each connection's quantities get a
 * suffix (its index, or its name if 'name_multiple_ht_by_index' is false) and the flow channel
 * combines them into the bare names.
 *
 * This class is currently only meant to be attached to a single flow channel component; using it
 * on more than one component would incorrectly apply one component's fluid properties/area/initial
 * conditions to all of them.
 */
class FlowChannel1PhasePhysics : public PhysicsBase,
                                 public GravityInterface,
                                 public NamingInterface,
                                 public FlowChannelClosuresInterface
{
public:
  static InputParameters validParams();
  FlowChannel1PhasePhysics(const InputParameters & params);

  /// Name of the user object that defines fluid properties for this flow channel
  const UserObjectName & fluidPropertiesName() const { return _fp_name; }
  /// Name of the numerical flux user object this Physics creates, used by boundary conditions
  const UserObjectName & numericalFluxUserObjectName() const { return _numerical_flux_name; }
  /// Names of the passive transport solution variables, which are the passive transport variables
  /// multiplied by the area
  const std::vector<VariableName> & passiveTransportSolutionVariableNames() const
  {
    return _passives_times_area_names;
  }
  /// Names of all the variables solved for by the flow equations (mass, momentum, energy and the
  /// passive transport variables)
  virtual std::vector<VariableName> solutionVariableNames() const;
  /// Name of the function/constant defining the flow channel area, used to build a default heated
  /// perimeter when a heat-transfer component does not provide one
  const FunctionName & areaFunctionName() const { return _area_fn_name; }

  /// Registers a heat-transfer connection with this flow channel, mirroring
  /// FlowChannelBase::addHeatTransferName(). Const with mutable members because it is called on the
  /// const reference a heat-transfer ActionComponent gets from FlowChannel1PhaseAC::getPhysics().
  /// @param T_wall_name Name of the wall temperature variable of a temperature-type connection, if
  ///                    it is not the default one that this flow channel names
  void addHeatTransferName(const std::string & ht_name,
                           bool is_temperature_type,
                           const std::optional<VariableName> & T_wall_name = std::nullopt) const;

  /// Suffix of the names of the wall quantities of a heat transfer connection (empty if it is the
  /// only one), mirroring FlowChannelBase::getHeatTransferNamesSuffix()
  std::string heatTransferSuffix(const std::string & ht_name) const;

  /// The closures objects named by this flow channel's 'closures' parameter, resolved from the
  /// ClosuresRegistry - shared with any heat-transfer connection, exactly as classic THM's
  /// HeatTransferBase reuses its flow channel's closures objects
  std::vector<std::shared_ptr<ClosuresBase>> getClosuresObjects() const;

  // FlowChannelClosuresInterface implementation ----
  virtual const std::vector<SubdomainName> & getClosuresBlocks() const override { return _blocks; }
  virtual const std::string & getClosuresName() const override { return name(); }
  virtual bool getClosuresTemperatureMode() const override { return _temperature_mode; }
  virtual unsigned int getClosuresNumberOfHeatTransferConnections() const override
  {
    return _heat_transfers.size();
  }
  virtual std::vector<VariableName> getClosuresHeatedPerimeterNames() const override;
  virtual std::vector<VariableName> getClosuresWallTemperatureNames() const override;
  virtual std::vector<MaterialPropertyName> getClosuresWallHTCNames() const override;
  virtual HeatTransferGeometry getClosuresHeatTransferGeometry() const override
  {
    return HeatTransferGeometry::PIPE;
  }
  virtual PipeLocation getClosuresPipeLocation() const override { return PipeLocation::INTERIOR; }
  virtual bool hasClosuresWallFrictionFactorFunction() const override { return isParamValid("f"); }
  virtual const FunctionName & getClosuresWallFrictionFactorFunction() const override
  {
    return getParam<FunctionName>("f");
  }
  virtual Real getClosuresRoughness() const override { return getParam<Real>("roughness"); }
  virtual Real getClosuresPoD() const override { return getParam<Real>("PoD"); }
  virtual void connectClosuresObject(const InputParameters & obj_params,
                                     const std::string & obj_name,
                                     const std::string & param) const override;

protected:
  /// The checks on the fluid properties of the flow channel
  virtual void checkFluidProperties() const;

  virtual void addSolverVariables() override;
  virtual void addAuxiliaryVariables() override;
  virtual void addInitialConditions() override;
  virtual void addFEKernels() override;
  virtual void addDGKernels() override;
  virtual void addMaterials() override;
  virtual void addUserObjects() override;
  virtual void addAuxiliaryKernels() override;
  virtual void checkIntegrity() const override;

  /// Adds the initial condition of rhoEA
  virtual void addRhoEAInitialCondition();
  /// Adds the initial condition of the density
  virtual void addDensityInitialCondition();
  /// Adds the DG kernels of the advection terms of the solution variables
  virtual void addAdvectionDGKernels();
  /// Adds the material that computes the fluid properties from the solution variables
  virtual void addFluidPropertiesMaterial();
  /// Adds the material that computes the slope reconstruction (the rDG scheme)
  virtual void addSlopeReconstructionMaterial();
  /// Adds the user object of the numerical flux
  virtual void addNumericalFluxUserObject();
  /// Adds the auxiliary kernels of the pressure and the temperature
  virtual void addPressureTemperatureAuxKernels();

  virtual std::vector<UserObjectName> getSuppliedUserObjects() const override
  {
    return {_numerical_flux_name};
  }

protected:
  /// Adds the hydraulic diameter material, needed by the wall friction kernel/closures
  /// correlations - mirrors FlowChannel1PhaseBase::addHydraulicDiameterMaterial()
  void addHydraulicDiameterMaterial();
  /// Adds a nonlinear (nl = true) or auxiliary (nl = false) flow variable on this Physics' blocks
  void addFlowVariable(bool nl,
                       const VariableName & var_name,
                       const std::string & family = "MONOMIAL",
                       const std::string & order = "CONSTANT",
                       Real scaling_factor = 1.0);
  /// Adds a FunctionIC for a variable
  void addFunctionIC(const VariableName & var_name, const FunctionName & function_name);
  /// Seeds AREA/AREA_LINEAR's initial value so other ICs (e.g. rhoA = density * area) can read it
  void addAreaInitialCondition();
  /// Adds an ADTimeDerivative kernel for the given variable if the physics is transient
  void addTimeDerivativeKernelIfTransient(const VariableName & var_name);

protected:
  /// Name of the user object that defines fluid properties
  const UserObjectName _fp_name;
  /// Area of the flow channel; can be a constant or a function
  const FunctionName _area_fn_name;
  /// Initial pressure, temperature and velocity functions
  const FunctionName _initial_p_fn;
  const FunctionName _initial_T_fn;
  const FunctionName _initial_vel_fn;
  /// Scaling factors for rhoA, rhouA, rhoEA
  const std::vector<Real> _scaling_factors;
  /// Names of the passive transport solution variables (the passives multiplied by the area)
  std::vector<VariableName> _passives_times_area_names;
  /// Name of the numerical flux user object this Physics creates
  const UserObjectName _numerical_flux_name;
  /// Coordinator that merges this and other components' requests for shared, bare-named
  /// variables (rhoA, A, ...) into single MOOSE variables spanning the union of their blocks
  THMVariableCoordinator & _coordinator;

  /// A connected heat-transfer component
  struct HeatTransferConnection
  {
    /// Name of the heat-transfer component
    std::string _name;
    /// Whether it specifies a wall temperature (as opposed to a heat flux)
    bool _temperature_type;
    /// Name of its wall temperature variable, if it is not the one named by this flow channel
    std::optional<VariableName> _T_wall_name;
  };

private:
  /// Suffix of the names of the wall quantities of the i-th heat transfer connection
  std::string connectionSuffix(unsigned int i) const;
  /// Names of the wall quantities of all the heat transfer connections, with the connection suffix
  std::vector<std::string> connectionNames(const std::string & base_name) const;

  /// The connected heat-transfer components. Mutable because addHeatTransferName() is called
  /// through a const reference (see that method).
  mutable std::vector<HeatTransferConnection> _heat_transfers;
  /// Whether any connected heat-transfer component specifies a wall temperature. Mutable for the
  /// same reason as _heat_transfers.
  mutable bool _temperature_mode;
};
