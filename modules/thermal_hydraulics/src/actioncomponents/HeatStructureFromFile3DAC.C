//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatStructureFromFile3DAC.h"
#include "THMUtils.h"
#include "MooseUtils.h"
#include "libmesh/exodusII_io_helper.h"

// Registered under the classic Component's name ('HeatStructureFromFile3D'); see
// FlowChannel1PhaseAC.C for why.
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructureFromFile3DAC,
                           "HeatStructureFromFile3D",
                           "add_mesh_generator");
registerMooseActionAliased("ThermalHydraulicsApp",
                           HeatStructureFromFile3DAC,
                           "HeatStructureFromFile3D",
                           "init_component_physics");
registerActionComponentAliased("ThermalHydraulicsApp",
                               HeatStructureFromFile3DAC,
                               "HeatStructureFromFile3D");

InputParameters
HeatStructureFromFile3DAC::validParams()
{
  InputParameters params = HeatStructureInterfaceAC::validParams();
  params += ComponentMeshTransformHelper::validParams();
  params.addRequiredParam<FileName>("file", "The ExodusII mesh file name");
  params.renameParam("direction",
                     "orientation",
                     "Direction to orient the heat structure with, assuming it is initially "
                     "oriented along the X-axis (1, 0, 0).");
  params.addClassDescription("Heat structure component that loads a 3D mesh from an ExodusII file");
  return params;
}

HeatStructureFromFile3DAC::HeatStructureFromFile3DAC(const InputParameters & params)
  : ActionComponent(params), HeatStructureInterfaceAC(params), ComponentMeshTransformHelper(params)
{
  THM::informExperimental(*this);

  const auto & file_name = getParam<FileName>("file");
  // This is an error instead of a check at a later task because the variables are added to this
  // component's blocks before the integrity check, so an unreadable file would otherwise only
  // result in an error about variables being added to non-existent blocks.
  if (!MooseUtils::pathExists(file_name) || !MooseUtils::checkFileReadable(file_name, false, false))
    paramError("file",
               "The file '",
               file_name,
               "' could not be opened. Check that the file exists and that you have permission to "
               "open it.");

  libMesh::ExodusII_IO_Helper exio_helper(*this, false, true, false);
  exio_helper.open(file_name.c_str(), true);
  exio_helper.read_and_store_header_info();
  if (exio_helper.num_dim != 3)
    paramError("file", "File '", file_name, "' does not contain a 3D mesh.");

  // The block or sideset ID is used if the file does not give it a name
  exio_helper.read_block_info();
  for (const auto i : make_range(exio_helper.num_elem_blk))
  {
    auto block_name = exio_helper.get_block_name(i);
    if (block_name.empty())
      block_name = Moose::stringify(exio_helper.get_block_id(i));
    _region_names.push_back(block_name);
  }
  exio_helper.read_sideset_info();
  for (const auto i : make_range(exio_helper.num_side_sets))
  {
    auto sideset_name = exio_helper.get_side_set_name(i);
    if (sideset_name.empty())
      sideset_name = Moose::stringify(exio_helper.get_side_set_id(i));
    _file_boundary_names.push_back(sideset_name);
  }

  _dimension = 3;
  for (const auto & region_name : _region_names)
    _blocks.push_back(genName(name(), region_name));

  addRequiredTask("add_mesh_generator");
  addRequiredTask("init_component_physics");

  createPhysics(params, /*second_order=*/false);
}

bool
HeatStructureFromFile3DAC::hasRegion(const std::string & region) const
{
  return std::find(_region_names.begin(), _region_names.end(), region) != _region_names.end();
}

bool
HeatStructureFromFile3DAC::hasBoundary(const BoundaryName & boundary) const
{
  for (const auto & file_boundary_name : _file_boundary_names)
    if (boundary == genName(name(), file_boundary_name))
      return true;
  return false;
}

void
HeatStructureFromFile3DAC::addMeshGenerators()
{
  const std::string file_mg_name = genSafeName(name(), "file");
  {
    InputParameters params = _factory.getValidParams("FileMeshGenerator");
    params.set<MeshFileName>("file") = getParam<FileName>("file");
    _app.getMeshGeneratorSystem().addMeshGenerator("FileMeshGenerator", file_mg_name, params);
    _mg_names.push_back(file_mg_name);
  }

  const std::string block_mg_name = genSafeName(name(), "block_names");
  {
    InputParameters params = _factory.getValidParams("RenameBlockGenerator");
    params.set<MeshGeneratorName>("input") = file_mg_name;
    params.set<std::vector<SubdomainName>>("old_block") =
        std::vector<SubdomainName>(_region_names.begin(), _region_names.end());
    params.set<std::vector<SubdomainName>>("new_block") = _blocks;
    _app.getMeshGeneratorSystem().addMeshGenerator("RenameBlockGenerator", block_mg_name, params);
    _mg_names.push_back(block_mg_name);
  }

  const std::string boundary_mg_name = genSafeName(name(), "boundary_names");
  {
    std::vector<BoundaryName> old_boundary, new_boundary;
    for (const auto & file_boundary_name : _file_boundary_names)
    {
      old_boundary.push_back(file_boundary_name);
      new_boundary.push_back(genName(name(), file_boundary_name));
    }
    InputParameters params = _factory.getValidParams("RenameBoundaryGenerator");
    params.set<MeshGeneratorName>("input") = block_mg_name;
    params.set<std::vector<BoundaryName>>("old_boundary") = old_boundary;
    params.set<std::vector<BoundaryName>>("new_boundary") = new_boundary;
    _app.getMeshGeneratorSystem().addMeshGenerator(
        "RenameBoundaryGenerator", boundary_mg_name, params);
    _mg_names.push_back(boundary_mg_name);
  }

  _top_mg_name = boundary_mg_name;

  ComponentMeshTransformHelper::addMeshGenerators();
}
