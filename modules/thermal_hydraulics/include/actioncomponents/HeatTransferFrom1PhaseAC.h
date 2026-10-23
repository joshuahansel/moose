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
#include "HeatTransferClosuresInterface.h"
#include "MooseTypes.h"
#include "THMNames.h"

#include <optional>

class FlowChannel1PhasePhysics;

/// Registers a HeatTransferFrom1PhaseAC-derived component to all the tasks it acts on
#define registerHeatTransferFrom1PhaseACTasks(app_name, derived_name, alias)                       \
  registerMooseActionAliased(app_name, derived_name, alias, "THM:add_heat_transfers");             \
  registerMooseActionAliased(app_name, derived_name, alias, "add_function");                       \
  registerMooseActionAliased(app_name, derived_name, alias, "add_aux_variable");                   \
  registerMooseActionAliased(app_name, derived_name, alias, "add_ics_physics");                    \
  registerMooseActionAliased(app_name, derived_name, alias, "add_material");                       \
  registerMooseActionAliased(app_name, derived_name, alias, "add_kernel");                         \
  registerMooseActionAliased(app_name, derived_name, alias, "add_aux_kernel");                     \
  registerMooseActionAliased(app_name, derived_name, alias, "check_integrity");                    \
  registerActionComponentAliased(app_name, derived_name, alias)

/**
 * Base class for a heat-transfer connection ActionComponent attaching a heat source (a
 * function-specified wall heat flux or wall temperature) to one FlowChannel1PhaseAC, addressed the
 * same way classic THM's HeatTransferBase connects: through a 'flow_channel' parameter naming the
 * flow channel component.
 *
 * This is the native-ActionComponent counterpart of classic THM's HeatTransfer1PhaseBase. It
 * implements HeatTransferClosuresInterface (identically) so the same ClosuresBase-derived objects
 * the flow channel uses supply the wall heat transfer coefficient material. The flow channel's own
 * closures objects are reused (via FlowChannel1PhasePhysics::getClosuresObjects()), never
 * re-created here - matching classic HeatTransferBase.
 *
 * Derived leaves provide only the source-specific pieces through the addHeatTransfer* hooks: the
 * heat-flux leaf adds a q_wall material and an ADOneDEnergyWallHeatFlux kernel; the
 * specified-temperature leaf adds a T_wall variable/IC/aux kernel and an ADOneDEnergyWallHeating
 * kernel. The heated-perimeter (P_hf) variable/IC/aux kernel common to both is handled here.
 */
class HeatTransferFrom1PhaseAC : public virtual ActionComponent,
                                 public NamingInterface,
                                 public HeatTransferClosuresInterface
{
public:
  static InputParameters validParams();
  HeatTransferFrom1PhaseAC(const InputParameters & params);

  // HeatTransferClosuresInterface implementation (identical to HeatTransfer1PhaseBase) ----
  virtual const std::string & getClosuresName() const override { return name(); }
  virtual bool hasClosuresWallHeatTransferCoefficientFunction() const override
  {
    return isParamValid("Hw");
  }
  virtual const FunctionName & getClosuresWallHeatTransferCoefficientFunction() const override
  {
    return getParam<FunctionName>("Hw");
  }
  virtual const MaterialPropertyName &
  getClosuresWallHeatTransferCoefficient1PhaseName() const override;

protected:
  virtual void actOnAdditionalTasks() override;
  virtual void addMaterials() override;

  /// Whether this connection specifies a wall temperature (true) or a wall heat flux (false).
  /// Mirrors HeatTransferBase::isTemperatureType(); drives the flow channel's temperature mode.
  virtual bool isTemperatureType() const = 0;

  /// Name of the wall temperature variable of a temperature-type connection if it is not created by
  /// the flow channel's naming (e.g. it is created by another application)
  virtual std::optional<VariableName> specifiedWallTemperatureName() const { return std::nullopt; }

  /// Suffix of the names of this connection's wall quantities, which is not empty if there is more
  /// than one heat transfer connection to the flow channel. Only valid after the task
  /// THM:add_heat_transfers, once all the connections are known.
  std::string suffix() const;
  /// Name of the heated perimeter variable of this connection
  VariableName heatedPerimeterName() const { return THM::HEAT_FLUX_PERIMETER + suffix(); }
  /// Name of the wall temperature variable of this connection
  VariableName wallTemperatureName() const;
  /// Name of the wall heat flux of this connection
  std::string wallHeatFluxName() const { return THM::HEAT_FLUX_WALL + suffix(); }
  /// Name of the wall heat transfer coefficient of this connection
  MaterialPropertyName wallHeatTransferCoefficientName() const
  {
    return THM::HEAT_TRANSFER_COEFFICIENT_WALL + suffix();
  }

  /// Physics of the connected flow channel component, looked up by name lazily so that
  /// 'flow_channel' may name a component declared before or after this one (mirrors
  /// FlowBoundary1PhaseAC::connectedPhysics()).
  const FlowChannel1PhasePhysics & connectedPhysics() const;

  /// Physics of all the flow channels this component connects to. Derived classes that connect
  /// several flow channels override this; the tasks of this class are applied to each of them.
  virtual std::vector<const FlowChannel1PhasePhysics *> connectedPhysicsList() const
  {
    return {&connectedPhysics()};
  }
  /// Union of the flow channel blocks of all the connected flow channels
  std::vector<SubdomainName> closuresBlocks() const;

  // Source-specific hooks for leaves (all no-ops by default) ----
  /// Add source-specific variables (e.g. the wall temperature variable)
  virtual void addHeatTransferVariables() {}
  /// Add source-specific initial conditions
  virtual void addHeatTransferICs() {}
  /// Add source-specific materials (e.g. the wall heat flux material)
  virtual void addHeatTransferMaterials() {}
  /// Add the source's energy-equation kernel
  virtual void addHeatTransferKernels() {}
  /// Add source-specific auxiliary kernels (e.g. the wall temperature aux kernel)
  virtual void addHeatTransferAuxKernels() {}

  /// The heated perimeter function name: the 'P_hf' parameter if given, otherwise the name of a
  /// GeneralizedCircumference function this component builds from the flow channel area
  const FunctionName _P_hf_fn_name;
  /// Wall heat transfer coefficient material property name, updated when it is requested because
  /// the suffix is not known until all the connections are
  mutable MaterialPropertyName _Hw_name;
};
