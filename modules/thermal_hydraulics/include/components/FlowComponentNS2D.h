//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Component2D.h"
#include "NSFVBase.h"
#include "THMMesh.h"

/**
 * Navier-Stokes flow component with a 2D generated mesh
 */
class FlowComponentNS2D : public NSFVBase<Component2D>
{
public:
  static InputParameters validParams();

  FlowComponentNS2D(const InputParameters & parameters);

  virtual void addRelationshipManagers(Moose::RelationshipManagerType input_rm_type) override;
  virtual void addVariables() override;
  virtual void addMooseObjects() override;

  virtual bool isCylindrical() const override { return _is_cylindrical; }
  virtual Real getDepth() const override { return _depth; }
  virtual Real getInnerRadius() const override { return _inner_radius; }

protected:
  virtual void init() override;
  virtual void check() const override;
  virtual bool usingSecondOrderMesh() const override { return false; }


  virtual std::vector<SubdomainName> getBlocks() const override { return getSubdomainNames(); }
  virtual Factory & getFactory() override { return getMooseApp().getFactory(); }
  virtual FEProblemBase & getProblem() override { return getMooseApp().feProblem(); }
  virtual const MooseMesh & getMesh() const override { return constMesh(); }
  virtual void addNSNonlinearVariable(const std::string & var_type,
                                      const std::string & var_name,
                                      InputParameters & params) override
  {
    getTHMProblem().addSimVariable(true, var_type, var_name, params);
  }
  virtual void addNSAuxVariable(const std::string & var_type,
                                const std::string & var_name,
                                InputParameters & params) override
  {
    getTHMProblem().addSimVariable(false, var_type, var_name, params);
  }
  virtual void addNSInitialCondition(const std::string & type,
                                     const std::string & name,
                                     InputParameters & params) override
  {
    getTHMProblem().addSimInitialCondition(type, name, params);
  }
  virtual std::string prefix() const override { return name() + ":"; }

  /// True if the component is cylindrical
  const bool _is_cylindrical;
  /// Depth if the component is not cylindrical
  const Real _depth;
  /// Depth if the component is cylindrical
  const Real _inner_radius;
};
