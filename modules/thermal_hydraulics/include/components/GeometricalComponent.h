//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Component.h"
#include "libmesh/enum_elem_type.h"

/**
 * Intermediate class for components that have mesh
 */
class GeometricalComponent : public Component
{
public:
  GeometricalComponent(const InputParameters & parameters);

  /**
   * Gets the subdomain names for this component
   *
   * @return vector of subdomain names for this component
   */
  virtual const std::vector<SubdomainName> & getSubdomainNames() const;

  /**
   * Gets the coordinate system types for this component
   *
   * @return vector of coordinate system types for this component
   */
  virtual const std::vector<Moose::CoordinateSystemType> & getCoordSysTypes() const;

  /**
   * Gets the name of each RZ subdomain in this component
   */
  const std::vector<SubdomainName> & getRZSubdomainNames() const;

  /**
   * Gets the RZ coordinate axis for each RZ subdomain in this component
   */
  const std::vector<std::pair<Point, RealVectorValue>> & getRZAxes() const;

  /**
   * Gets the node IDs corresponding to this component
   */
  const std::vector<dof_id_type> & getNodeIDs() const;

  /**
   * Gets the element IDs corresponding to this component
   */
  const std::vector<dof_id_type> & getElementIDs() const;

protected:
  Node * addNode(const Point & pt);

  Elem * addElement(libMesh::ElemType elem_type, const std::vector<dof_id_type> & node_ids);
  Elem * addElementEdge2(dof_id_type node0, dof_id_type node1);
  Elem * addElementEdge3(dof_id_type node0, dof_id_type node1, dof_id_type node2);
  Elem *
  addElementQuad4(dof_id_type node0, dof_id_type node1, dof_id_type node2, dof_id_type node3);
  Elem * addElementQuad9(dof_id_type node0,
                         dof_id_type node1,
                         dof_id_type node2,
                         dof_id_type node3,
                         dof_id_type node4,
                         dof_id_type node5,
                         dof_id_type node6,
                         dof_id_type node7,
                         dof_id_type node8);

  /**
   * Sets the subdomain ID, name, and coordinate system for an XYZ subdomain
   *
   * @param[in] subdomain_id  subdomain index
   * @param[in] subdomain_name  name of the new subdomain
   */
  virtual void setSubdomainInfoXYZ(SubdomainID subdomain_id, const std::string & subdomain_name);

  /**
   * Sets the subdomain ID, name, coordinate system, and RZ axis data for an RZ_GENERAL subdomain
   *
   * @param[in] subdomain_id  subdomain index
   * @param[in] subdomain_name  name of the new subdomain
   * @param[in] axis  RZ axis
   */
  virtual void setSubdomainInfoRZ(SubdomainID subdomain_id,
                                  const SubdomainName & subdomain_name,
                                  const std::pair<Point, RealVectorValue> & axis);

  /**
   * Makes a constant function parameter controllable and returns its name
   *
   * @param[in] fn_param_name   FunctionName parameter
   */
  const FunctionName & getVariableFn(const FunctionName & fn_param_name);

  /// List of subdomain IDs this components owns
  std::vector<SubdomainID> _subdomain_ids;
  /// List of subdomain names this components owns
  std::vector<SubdomainName> _subdomain_names;
  /// List of coordinate system for each subdomain
  std::vector<Moose::CoordinateSystemType> _coord_sys;

  /// List of RZ subdomain names
  std::vector<SubdomainName> _rz_subdomain_names;
  /// List of RZ axes
  std::vector<std::pair<Point, RealVectorValue>> _rz_axes;

  /// List of node IDs this components owns
  std::vector<dof_id_type> _node_ids;
  /// Elements ids of this component
  std::vector<dof_id_type> _elem_ids;

public:
  static InputParameters validParams();
};
