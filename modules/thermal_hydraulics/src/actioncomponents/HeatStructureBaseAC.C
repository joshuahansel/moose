//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatStructureBaseAC.h"
#include "THMUtils.h"
#include "ActionFactory.h"

#include <numeric>

InputParameters
HeatStructureBaseAC::validParams()
{
  InputParameters params = HeatStructureInterfaceAC::validParams();
  params += ComponentMeshTransformHelper::validParams();
  params.addClassDescription(
      "Base class for 2D heat-conducting solid ActionComponents, with a mesh generated through "
      "the MeshGenerator system.");

  params.addRequiredParam<std::vector<std::string>>("names", "Name of each transverse region");
  params.addRequiredParam<std::vector<Real>>("widths", "Width of each transverse region [m]");
  params.addRequiredParam<std::vector<unsigned int>>(
      "n_part_elems", "Number of elements of each transverse region");
  params.addParam<std::vector<UserObjectName>>(
      "solid_properties", "Solid properties object name for each transverse region");
  params.addParam<std::vector<Real>>(
      "solid_properties_T_ref",
      {},
      "Density reference temperatures for each transverse region. This is required if "
      "'solid_properties' is provided. The density in each region will be a constant value "
      "computed by evaluating the density function at the reference temperature.");
  params.addParam<Real>("num_rods", 1.0, "Number of rods represented by this heat structure");

  params.addRequiredRangeCheckedParam<std::vector<Real>>(
      "length", "length > 0", "Length of each axial section of the heat structure [m]");
  params.addRequiredRangeCheckedParam<std::vector<unsigned int>>(
      "n_elems", "n_elems > 0", "Number of elements of each axial section of the heat structure");
  params.addParam<bool>("2nd_order_mesh", false, "Use 2nd order elements in the mesh");
  params.addParam<std::vector<std::string>>(
      "axial_region_names",
      {},
      "Names of the axial sections. This is required if there is more than 1 axial section.");

  params.renameParam("direction",
                     "orientation",
                     "Direction to orient the heat structure with, assuming it is initially "
                     "oriented along the X-axis (1, 0, 0).");

  return params;
}

HeatStructureBaseAC::HeatStructureBaseAC(const InputParameters & params)
  : ActionComponent(params),
    HeatStructureInterfaceAC(params),
    ComponentMeshTransformHelper(params),
    _region_names(getParam<std::vector<std::string>>("names")),
    _widths(getParam<std::vector<Real>>("widths")),
    _n_part_elems(getParam<std::vector<unsigned int>>("n_part_elems")),
    _num_rods(getParam<Real>("num_rods")),
    _lengths(getParam<std::vector<Real>>("length")),
    _n_elems(getParam<std::vector<unsigned int>>("n_elems")),
    _axial_region_names(getParam<std::vector<std::string>>("axial_region_names")),
    _length(std::accumulate(_lengths.begin(), _lengths.end(), 0.0)),
    _second_order(getParam<bool>("2nd_order_mesh")),
    _has_axial_regions(_lengths.size() > 1 && _axial_region_names.size() == _lengths.size())
{
  THM::informExperimental(*this);
  if (_widths.size() != _region_names.size() || _n_part_elems.size() != _region_names.size())
    paramError("names", "'names', 'widths' and 'n_part_elems' must all be the same size.");
  if (_n_elems.size() != _lengths.size())
    paramError("n_elems",
               "The number of entries in parameter 'n_elems' (",
               _n_elems.size(),
               ") must equal the number of entries of parameter 'length' (",
               _lengths.size(),
               ")");
  if (_axial_region_names.size())
  {
    if (_axial_region_names.size() != _lengths.size())
      paramError("axial_region_names",
                 "The number of entries in parameter 'axial_region_names' (",
                 _axial_region_names.size(),
                 ") must equal the number of entries of parameter 'length' (",
                 _lengths.size(),
                 ")");
  }
  else if (_lengths.size() > 1)
    paramError("length",
               "If there is more than 1 axial region, then the parameter 'axial_region_names' must "
               "be specified.");
  if (isParamValid("solid_properties"))
  {
    const auto sp_names = getParam<std::vector<UserObjectName>>("solid_properties");
    const auto T_ref = getParam<std::vector<Real>>("solid_properties_T_ref");
    if (sp_names.size() != _region_names.size() || T_ref.size() != _region_names.size())
      paramError("solid_properties",
                 "'solid_properties' and 'solid_properties_T_ref' must have as many entries as "
                 "'names'.");
  }

  _dimension = 2;
  for (const auto & region_name : _region_names)
    _blocks.push_back(genName(name(), region_name));

  addRequiredTask("add_mesh_generator");
  addRequiredTask("init_component_physics");
  addRequiredTask("add_material");
  addRequiredTask("check_integrity");

  createPhysics(params, _second_order);
}

void
HeatStructureBaseAC::addMeshGenerators()
{
  std::vector<MeshGeneratorName> region_mg_names;
  for (unsigned int i = 0; i < _region_names.size(); i++)
  {
    const std::string base_name = genSafeName(name(), _region_names[i], "base");
    InputParameters params = _factory.getValidParams("CartesianMeshGenerator");
    params.set<MooseEnum>("dim") = "2";
    params.set<std::vector<Real>>("dx") = _lengths;
    params.set<std::vector<unsigned int>>("ix") = _n_elems;
    params.set<std::vector<Real>>("dy") = {_widths[i]};
    params.set<std::vector<unsigned int>>("iy") = {_n_part_elems[i]};
    // StackGenerator (unlike CombinerGenerator) does not disambiguate colliding subdomain ids when
    // stitching meshes together, and CartesianMeshGenerator defaults every region to id 0 - so an
    // explicit, unique id per region (the same for each of its axial sections) is required here.
    params.set<std::vector<unsigned int>>("subdomain_id") =
        std::vector<unsigned int>(_lengths.size(), i);
    _app.getMeshGeneratorSystem().addMeshGenerator("CartesianMeshGenerator", base_name, params);
    _mg_names.push_back(base_name);
    region_mg_names.push_back(base_name);
  }

  // Stack the regions along the transverse ('y') direction. All regions were left with their
  // default 'top'/'bottom' boundary names, so StackGenerator can stitch consecutive regions
  // together automatically; the outermost region's 'top' and innermost region's 'bottom' survive
  // un-stitched. StackGenerator requires at least 2 inputs (it indexes into a "rest of the inputs"
  // vector unconditionally), so a single region is instead just translated by the same offset.
  std::string stack_name;
  if (region_mg_names.size() > 1)
  {
    stack_name = genSafeName(name(), "stack");
    InputParameters params = _factory.getValidParams("StackGenerator");
    params.set<MooseEnum>("dim") = "2";
    params.set<std::vector<MeshGeneratorName>>("inputs") = region_mg_names;
    params.set<Real>("bottom_height") = yOffset();
    _app.getMeshGeneratorSystem().addMeshGenerator("StackGenerator", stack_name, params);
  }
  else
  {
    stack_name = genSafeName(name(), "stack");
    InputParameters params = _factory.getValidParams("TransformGenerator");
    params.set<MeshGeneratorName>("input") = region_mg_names[0];
    params.set<MooseEnum>("transform") = "TRANSLATE";
    params.set<RealVectorValue>("vector_value") = RealVectorValue(0, yOffset(), 0);
    _app.getMeshGeneratorSystem().addMeshGenerator("TransformGenerator", stack_name, params);
  }
  _mg_names.push_back(stack_name);

  // stitch_meshes (used by StackGenerator) preserves each region's numeric subdomain id (made
  // unique above) but does not carry over their subdomain *names* from the meshes it consumes -
  // only the primary (first) region's name survives. Re-establish every region's name by its
  // now-guaranteed-unique numeric id.
  const std::string block_renamed_name = genSafeName(name(), "block_names");
  {
    std::vector<SubdomainName> old_block;
    for (unsigned int i = 0; i < _region_names.size(); i++)
      old_block.push_back(std::to_string(i));

    InputParameters params = _factory.getValidParams("RenameBlockGenerator");
    params.set<MeshGeneratorName>("input") = stack_name;
    params.set<std::vector<SubdomainName>>("old_block") = old_block;
    params.set<std::vector<SubdomainName>>("new_block") = _blocks;
    _app.getMeshGeneratorSystem().addMeshGenerator(
        "RenameBlockGenerator", block_renamed_name, params);
    _mg_names.push_back(block_renamed_name);
  }

  // Rename the surviving 'top'/'bottom' boundaries to the aggregate inner/outer names.
  const std::string renamed_name = genSafeName(name(), "outer_inner");
  {
    InputParameters params = _factory.getValidParams("RenameBoundaryGenerator");
    params.set<MeshGeneratorName>("input") = block_renamed_name;
    params.set<std::vector<BoundaryName>>("old_boundary") = {"top", "bottom"};
    params.set<std::vector<BoundaryName>>("new_boundary") = {genName(name(), "outer"),
                                                             genName(name(), "inner")};
    _app.getMeshGeneratorSystem().addMeshGenerator("RenameBoundaryGenerator", renamed_name, params);
    _mg_names.push_back(renamed_name);
  }

  // Tag the start/end faces, per region and in aggregate (spanning every region), directly from
  // geometry on the now-stacked mesh - SideSetsAroundSubdomainGenerator before stacking would only
  // have survived for the primary (first) region, the same way subdomain names did not survive
  // for the others (see the RenameBlockGenerator step above).
  std::string current = renamed_name;
  std::vector<std::pair<std::string, std::vector<SubdomainName>>> start_end_blocks;
  for (unsigned int i = 0; i < _region_names.size(); i++)
    start_end_blocks.emplace_back(_region_names[i], std::vector<SubdomainName>{_blocks[i]});
  start_end_blocks.emplace_back("", _blocks);
  for (const auto & name_and_blocks : start_end_blocks)
    for (const auto & side_normal :
         {std::make_pair(std::string("start"), RealVectorValue(-1, 0, 0)),
          std::make_pair(std::string("end"), RealVectorValue(1, 0, 0))})
    {
      const std::string step_name = genSafeName(
          name(), name_and_blocks.first.empty() ? "agg" : name_and_blocks.first, side_normal.first);
      InputParameters params = _factory.getValidParams("SideSetsAroundSubdomainGenerator");
      params.set<MeshGeneratorName>("input") = current;
      params.set<std::vector<SubdomainName>>("block") = name_and_blocks.second;
      params.set<Point>("normal") = side_normal.second;
      params.set<std::vector<BoundaryName>>("new_boundary") = {
          name_and_blocks.first.empty()
              ? genName(name(), side_normal.first)
              : genName(name(), name_and_blocks.first, side_normal.first)};
      _app.getMeshGeneratorSystem().addMeshGenerator(
          "SideSetsAroundSubdomainGenerator", step_name, params);
      _mg_names.push_back(step_name);
      current = step_name;
    }

  // Tag the interior interfaces between adjacent regions.
  for (unsigned int i = 0; i + 1 < _region_names.size(); i++)
  {
    const std::string step_name = genSafeName(name(), "interface", std::to_string(i));
    InputParameters params = _factory.getValidParams("SideSetsBetweenSubdomainsGenerator");
    params.set<MeshGeneratorName>("input") = current;
    params.set<std::vector<SubdomainName>>("primary_block") = {_blocks[i]};
    params.set<std::vector<SubdomainName>>("paired_block") = {_blocks[i + 1]};
    params.set<std::vector<BoundaryName>>("new_boundary") = {
        genName(name(), _region_names[i], _region_names[i + 1])};
    _app.getMeshGeneratorSystem().addMeshGenerator(
        "SideSetsBetweenSubdomainsGenerator", step_name, params);
    _mg_names.push_back(step_name);
    current = step_name;
  }

  if (_has_axial_regions)
  {
    // Axial sections' end positions along the (not yet transformed) axial direction
    std::vector<Real> x_ends;
    Real x = 0;
    for (const auto length : _lengths)
      x_ends.push_back(x += length);

    // Sides of a ParsedGenerateSideset are selected by their centroid, so this tolerance only needs
    // to be small compared to the element size; a relative tolerance on the total length is used.
    const Real x_tol = 1e-8 * _length;

    // Per axial region, the inner and outer boundaries
    for (const auto i : index_range(_axial_region_names))
      for (const auto & side : {"inner", "outer"})
      {
        const std::string step_name = genSafeName(name(), _axial_region_names[i], side);
        InputParameters params = _factory.getValidParams("ParsedGenerateSideset");
        params.set<MeshGeneratorName>("input") = current;
        params.set<std::vector<BoundaryName>>("included_boundaries") = {genName(name(), side)};
        params.set<std::string>("combinatorial_geometry") =
            "x > " + Moose::stringifyExact(x_ends[i] - _lengths[i]) + " & x < " +
            Moose::stringifyExact(x_ends[i]);
        params.set<BoundaryName>("new_sideset_name") =
            genName(name(), _axial_region_names[i], side);
        _app.getMeshGeneratorSystem().addMeshGenerator("ParsedGenerateSideset", step_name, params);
        _mg_names.push_back(step_name);
        current = step_name;
      }

    // Per radial region, the interior axial interfaces between consecutive axial regions. The
    // faces are tagged from the side of the element on the lower-x side of the interface.
    for (const auto i : make_range(_axial_region_names.size() - 1))
      for (const auto j : index_range(_region_names))
      {
        const std::string interface_name =
            _axial_region_names[i] + ":" + _axial_region_names[i + 1];
        const std::string step_name = genSafeName(name(), _region_names[j], interface_name);
        InputParameters params = _factory.getValidParams("ParsedGenerateSideset");
        params.set<MeshGeneratorName>("input") = current;
        params.set<std::vector<SubdomainName>>("included_subdomains") = {_blocks[j]};
        params.set<Point>("normal") = Point(1, 0, 0);
        params.set<std::string>("combinatorial_geometry") =
            "abs(x - " + Moose::stringifyExact(x_ends[i]) + ") < " + Moose::stringifyExact(x_tol);
        params.set<BoundaryName>("new_sideset_name") =
            genName(name(), _region_names[j], interface_name);
        _app.getMeshGeneratorSystem().addMeshGenerator("ParsedGenerateSideset", step_name, params);
        _mg_names.push_back(step_name);
        current = step_name;
      }
  }

  // Every region's own default 'left'/'right' boundary tags (never renamed - see above) are left
  // dangling once the aggregate/per-region start/end boundaries have been derived from geometry.
  // They are redundant (each spans the same faces as one of the boundaries just derived) and,
  // being unprefixed, collide by name with every other heat structure's own leftover 'left'/
  // 'right' tags once multiple instances are combined into one simulation - so delete them here.
  {
    const std::string step_name = genSafeName(name(), "cleanup");
    InputParameters params = _factory.getValidParams("BoundaryDeletionGenerator");
    params.set<MeshGeneratorName>("input") = current;
    params.set<std::vector<BoundaryName>>("boundary_names") = {"left", "right"};
    _app.getMeshGeneratorSystem().addMeshGenerator("BoundaryDeletionGenerator", step_name, params);
    _mg_names.push_back(step_name);
    current = step_name;
  }

  if (_second_order)
  {
    const std::string step_name = genSafeName(name(), "2nd_order");
    InputParameters params = _factory.getValidParams("ElementOrderConversionGenerator");
    params.set<MeshGeneratorName>("input") = current;
    params.set<MooseEnum>("conversion_type") = "SECOND_ORDER";
    _app.getMeshGeneratorSystem().addMeshGenerator(
        "ElementOrderConversionGenerator", step_name, params);
    _mg_names.push_back(step_name);
    current = step_name;
  }

  _top_mg_name = current;

  ComponentMeshTransformHelper::addMeshGenerators();
}

void
HeatStructureBaseAC::addConstantDensitySolidPropertiesMaterial(const UserObjectName & sp_name,
                                                               Real T_ref,
                                                               const SubdomainName & block)
{
  const std::string class_name = "ADConstantDensityThermalSolidPropertiesMaterial";
  InputParameters params = _factory.getValidParams(class_name);
  params.set<std::vector<SubdomainName>>("block") = {block};
  params.set<std::vector<VariableName>>("temperature") = {"T_solid"};
  params.set<UserObjectName>("sp") = sp_name;
  params.set<Real>("T_ref") = T_ref;
  getProblem().addMaterial(class_name, genSafeName(name(), class_name, block), params);
}

void
HeatStructureBaseAC::addMaterials()
{
  if (!isParamValid("solid_properties"))
    return;

  const auto sp_names = getParam<std::vector<UserObjectName>>("solid_properties");
  const auto T_ref = getParam<std::vector<Real>>("solid_properties_T_ref");
  for (unsigned int i = 0; i < sp_names.size(); i++)
    addConstantDensitySolidPropertiesMaterial(sp_names[i], T_ref[i], _blocks[i]);
}

void
HeatStructureBaseAC::checkIntegrity()
{
  // The TRAP quadrature rule cannot integrate 2nd order elements
  if (_second_order)
    for (const auto * action : _awh.getActionListByName("setup_quadrature"))
      // Not every action on this task (e.g. THMSetupQuadratureAction) selects a quadrature type
      if (action->parameters().have_parameter<MooseEnum>("type") &&
          action->getParam<MooseEnum>("type") == "TRAP")
        mooseError(name(),
                   ": Cannot use TRAP quadrature rule with 2nd order elements.  Use SIMPSON or "
                   "GAUSS instead.");
}

bool
HeatStructureBaseAC::hasRegion(const std::string & region) const
{
  return std::find(_region_names.begin(), _region_names.end(), region) != _region_names.end();
}

bool
HeatStructureBaseAC::hasExternalBoundary(const BoundaryName & boundary) const
{
  for (const auto & side : {"inner", "outer", "start", "end"})
    if (boundary == genName(name(), side))
      return true;
  for (const auto & region_name : _region_names)
    for (const auto & side : {"start", "end"})
      if (boundary == genName(name(), region_name, side))
        return true;
  if (_has_axial_regions)
    for (const auto & axial_region_name : _axial_region_names)
      for (const auto & side : {"inner", "outer"})
        if (boundary == genName(name(), axial_region_name, side))
          return true;
  return false;
}

HeatStructureBaseAC::ExternalBoundaryType
HeatStructureBaseAC::getExternalBoundaryType(const BoundaryName & boundary) const
{
  if (boundary == genName(name(), "inner"))
    return ExternalBoundaryType::INNER;
  if (boundary == genName(name(), "outer"))
    return ExternalBoundaryType::OUTER;
  if (boundary == genName(name(), "start"))
    return ExternalBoundaryType::START;
  if (boundary == genName(name(), "end"))
    return ExternalBoundaryType::END;
  for (const auto & region_name : _region_names)
  {
    if (boundary == genName(name(), region_name, "start"))
      return ExternalBoundaryType::START;
    if (boundary == genName(name(), region_name, "end"))
      return ExternalBoundaryType::END;
  }
  if (_has_axial_regions)
    for (const auto & axial_region_name : _axial_region_names)
    {
      if (boundary == genName(name(), axial_region_name, "inner"))
        return ExternalBoundaryType::INNER;
      if (boundary == genName(name(), axial_region_name, "outer"))
        return ExternalBoundaryType::OUTER;
    }
  mooseError(name(), ": '", boundary, "' is not an external boundary of this heat structure.");
}
