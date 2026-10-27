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
#include "ComponentMeshTransformHelper.h"
#include "GravityInterface.h"
#include "ActionComponentConvergenceInterface.h"

class FlowChannel1PhasePhysics;
class ActionWarehouse;

/**
 * A 1D single-phase flow channel.
 *
 * The mesh is generated through the MeshGenerator system (see addMeshGenerators()). The governing
 * equations are defined by an internally-created FlowChannel1PhasePhysics, built and registered by
 * this component's constructor (rather than requiring the user to declare a separate [Physics]
 * block), so that an input file only needs to rename [Components] to [ActionComponents] and
 * FlowChannel1Phase to FlowChannel1PhaseAC to move a plain, unheated, frictionless flow channel
 * over to the native ActionComponent - see FlowChannel1PhasePhysics for what is and is not
 * currently covered.
 */
class FlowChannel1PhaseAC : public virtual ActionComponent,
                            public ComponentMeshTransformHelper,
                            public GravityInterface,
                            public ActionComponentConvergenceInterface
{
public:
  static InputParameters validParams();
  FlowChannel1PhaseAC(const InputParameters & params);

  /// The physics this component created for itself, e.g. for a boundary condition component to
  /// look up this flow channel's fluid properties/numerical flux user object names
  const FlowChannel1PhasePhysics & getPhysics() const;

  /// Looks up the flow channel that another component connects to by name, which is reported with
  /// the name of the connecting component \p from if there is no flow channel of that name
  static const FlowChannel1PhaseAC & connectedFlowChannel(const ActionWarehouse & awh,
                                                          const std::string & name,
                                                          const std::string & from);

  /// Length of the flow channel
  Real length() const { return _length; }
  /// Number of elements along the flow channel
  unsigned int numElems() const { return _n_elems; }

  /// Name of the nonlinear convergence of the flow channel
  virtual std::optional<ConvergenceName> nonlinearConvergenceName() const override;

protected:
  virtual void addMeshGenerators() override;
  virtual void addPhysics() override;
  virtual void actOnAdditionalTasks() override;

  /// Adds the objects that ComponentsConvergence uses: the step and residual of the pressure,
  /// temperature and velocity, and their combination into one Convergence object
  virtual void addConvergence();

  /// Coordinate (x, y or z) along which the vector postprocessor values are sorted: the dominant
  /// direction of the flow channel
  std::string sortBy() const;
  /// Adds the vector postprocessors that sample variables ('vpp_vars') and the numerical fluxes
  /// ('create_flux_vpp') along the flow channel
  void addVectorPostprocessors();

  /// Adds mesh generators that modify the straight mesh of the flow channel, which is generated
  /// along the x axis and has the name of the last mesh generator in _mg_names. Derived classes
  /// that are not straight override this; it is called before the rotation and the translation.
  virtual void addShapeMeshGenerators() {}

  /// Length of each axial section
  const std::vector<Real> _lengths;
  /// Number of elements of each axial section
  const std::vector<unsigned int> _n_elems_sections;
  /// Names of the axial sections (empty if the flow channel is one, unnamed section)
  const std::vector<std::string> _axial_region_names;
  /// Total length of the flow channel
  Real _length;
  /// Total number of elements along the flow channel
  const unsigned int _n_elems;

  /// The physics this component creates and owns; attaches itself to this component in addPhysics()
  FlowChannel1PhasePhysics * _physics;
};
