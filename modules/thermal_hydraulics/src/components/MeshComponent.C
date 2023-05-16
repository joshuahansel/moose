//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeshComponent.h"
#include "ConstantFunction.h"
#include "THMMesh.h"

InputParameters
MeshComponent::validParams()
{
  InputParameters params = Component::validParams();
  return params;
}

MeshComponent::MeshComponent(const InputParameters & parameters) : Component(parameters) {}

Node *
MeshComponent::addNode(const Point & pt)
{
  auto node = mesh().addNode(pt);
  _node_ids.push_back(node->id());
  return node;
}

Elem *
MeshComponent::addElement(libMesh::ElemType elem_type, const std::vector<dof_id_type> & node_ids)
{
  auto elem = mesh().addElement(elem_type, node_ids);
  _elem_ids.push_back(elem->id());
  return elem;
}

Elem *
MeshComponent::addElementEdge2(dof_id_type node0, dof_id_type node1)
{
  auto elem = mesh().addElementEdge2(node0, node1);
  _elem_ids.push_back(elem->id());
  return elem;
}

Elem *
MeshComponent::addElementEdge3(dof_id_type node0, dof_id_type node1, dof_id_type node2)
{
  auto elem = mesh().addElementEdge3(node0, node1, node2);
  _elem_ids.push_back(elem->id());
  return elem;
}

Elem *
MeshComponent::addElementQuad4(dof_id_type node0,
                               dof_id_type node1,
                               dof_id_type node2,
                               dof_id_type node3)
{
  auto elem = mesh().addElementQuad4(node0, node1, node2, node3);
  _elem_ids.push_back(elem->id());
  return elem;
}

Elem *
MeshComponent::addElementQuad9(dof_id_type node0,
                               dof_id_type node1,
                               dof_id_type node2,
                               dof_id_type node3,
                               dof_id_type node4,
                               dof_id_type node5,
                               dof_id_type node6,
                               dof_id_type node7,
                               dof_id_type node8)
{
  auto elem = mesh().addElementQuad9(node0, node1, node2, node3, node4, node5, node6, node7, node8);
  _elem_ids.push_back(elem->id());
  return elem;
}

const std::vector<SubdomainName> &
MeshComponent::getSubdomainNames() const
{
  checkSetupStatus(MESH_PREPARED);

  return _subdomain_names;
}

const std::vector<Moose::CoordinateSystemType> &
MeshComponent::getCoordSysTypes() const
{
  checkSetupStatus(MESH_PREPARED);

  return _coord_sys;
}

const std::vector<SubdomainName> &
MeshComponent::getRZSubdomainNames() const
{
  checkSetupStatus(MESH_PREPARED);

  return _rz_subdomain_names;
}

const std::vector<std::pair<Point, RealVectorValue>> &
MeshComponent::getRZAxes() const
{
  checkSetupStatus(MESH_PREPARED);

  return _rz_axes;
}

const FunctionName &
MeshComponent::getVariableFn(const FunctionName & fn_param_name)
{
  const FunctionName & fn_name = getParam<FunctionName>(fn_param_name);
  const Function & fn = getTHMProblem().getFunction(fn_name);

  if (dynamic_cast<const ConstantFunction *>(&fn) != nullptr)
  {
    connectObject(fn.parameters(), fn_name, fn_param_name, "value");
  }

  return fn_name;
}

void
MeshComponent::setSubdomainInfoXYZ(SubdomainID subdomain_id, const std::string & subdomain_name)
{
  _subdomain_ids.push_back(subdomain_id);
  _subdomain_names.push_back(subdomain_name);
  _coord_sys.push_back(Moose::COORD_XYZ);
  if (_parent)
  {
    MeshComponent * mesh_comp = dynamic_cast<MeshComponent *>(_parent);
    mesh_comp->_subdomain_ids.push_back(subdomain_id);
    mesh_comp->_subdomain_names.push_back(subdomain_name);
    mesh_comp->_coord_sys.push_back(Moose::COORD_XYZ);
  }
  mesh().setSubdomainName(subdomain_id, subdomain_name);
}

void
MeshComponent::setSubdomainInfoRZ(SubdomainID subdomain_id,
                                  const SubdomainName & subdomain_name,
                                  const std::pair<Point, RealVectorValue> & axis)
{
  _subdomain_ids.push_back(subdomain_id);
  _subdomain_names.push_back(subdomain_name);
  _coord_sys.push_back(Moose::COORD_RZ);
  _rz_subdomain_names.push_back(subdomain_name);
  _rz_axes.push_back(axis);
  if (_parent)
  {
    MeshComponent * mesh_comp = dynamic_cast<MeshComponent *>(_parent);
    mesh_comp->_subdomain_ids.push_back(subdomain_id);
    mesh_comp->_subdomain_names.push_back(subdomain_name);
    mesh_comp->_coord_sys.push_back(Moose::COORD_RZ);
    mesh_comp->_rz_subdomain_names.push_back(subdomain_name);
    mesh_comp->_rz_axes.push_back(axis);
  }
  mesh().setSubdomainName(subdomain_id, subdomain_name);
}

const std::vector<dof_id_type> &
MeshComponent::getNodeIDs() const
{
  checkSetupStatus(MESH_PREPARED);

  return _node_ids;
}

const std::vector<dof_id_type> &
MeshComponent::getElementIDs() const
{
  checkSetupStatus(MESH_PREPARED);

  return _elem_ids;
}
