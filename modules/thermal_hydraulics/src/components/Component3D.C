//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Component3D.h"
#include "THMMesh.h"
#include "THMEnums.h"

// const std::map<std::string, Component3D::ExternalBoundaryType>
//     Component3D::_external_boundary_type_to_enum{{"INNER", ExternalBoundaryType::INNER},
//                                                  {"OUTER", ExternalBoundaryType::OUTER},
//                                                  {"START", ExternalBoundaryType::START},
//                                                  {"END", ExternalBoundaryType::END}};

// MooseEnum
// Component3D::getExternalBoundaryTypeMooseEnum(const std::string & name)
// {
//   return THM::getMooseEnum<ExternalBoundaryType>(name, _external_boundary_type_to_enum);
// }

// template <>
// Component3D::ExternalBoundaryType
// THM::stringToEnum(const std::string & s)
// {
//   return stringToEnum<Component3D::ExternalBoundaryType>(
//       s, Component3D::_external_boundary_type_to_enum);
// }

InputParameters
Component3D::validParams()
{
  InputParameters params = GeneratedMeshComponent::validParams();

  params.addRequiredParam<std::vector<std::string>>("radial_region_names", "Name of each radial region");
  params.addRequiredParam<std::vector<Real>>("radial_region_widths", "Width of each radial region [m]");
  params.addRequiredParam<std::vector<unsigned int>>("n_radial_elems",
                                                     "Number of elements of each radial region");
  params.addParam<Real>("inner_radius", 0., "Inner radius of the componet [m]");
  params.addRequiredParam<unsigned int>("n_azimuthal_elems", "Number of azimuthal divisions");
  return params;
}

Component3D::Component3D(const InputParameters & params)
  : GeneratedMeshComponent(params), _n_radial_regions(0), _total_radial_elems(0), _axial_offset(0.0), _n_azimuthal_elems(getParam<unsigned int>("n_azimuthal_elems")), _azimuthal_elem_width(2 * libMesh::pi / _n_azimuthal_elems)
{
  _radial_region_names = getParam<std::vector<std::string>>("radial_region_names");
  _n_radial_regions = _radial_region_names.size();
  // for (unsigned int i = 0; i < _radial_region_names.size(); i++)
  //   _name_index[_radial_region_names[i]] = i;

  _radial_region_widths = getParam<std::vector<Real>>("radial_region_widths");
  _total_width = std::accumulate(_radial_region_widths.begin(), _radial_region_widths.end(), 0.0);

  _n_radial_elems = getParam<std::vector<unsigned int>>("n_radial_elems");
  for (unsigned int j_section = 0; j_section < _n_radial_elems.size(); j_section++)
    _total_radial_elems += _n_radial_elems[j_section];

  _inner_radius = getParam<Real>("inner_radius");
  _axial_offset = _inner_radius;

  // store radial positions
  _radial_positions.resize(_total_radial_elems + 1);
  _radial_positions[0] = _inner_radius;
  unsigned int j = 0;
  for (unsigned int j_section = 0; j_section < _n_radial_elems.size(); j_section++)
  {
    const Real region_elem_width = _radial_region_widths[j_section] / _n_radial_elems[j_section];
    for (unsigned int j_local = 0; j_local < _n_radial_elems[j_section]; j_local++)
    {
      _radial_positions[j + 1] = _radial_positions[j] + region_elem_width;
      j++;
    }
  }

  if (_radial_region_widths.size() == _n_radial_regions)
  {
    std::vector<Real> r(_n_radial_regions + 1, _inner_radius);
    for (unsigned int i = 0; i < _n_radial_regions; i++)
    {
      r[i + 1] = r[i] + _radial_region_widths[i];
      _radial_region_volumes.push_back(M_PI * (r[i + 1] * r[i + 1] - r[i] * r[i]) * _length);
    }
  }
}

void
Component3D::check() const
{
  GeneratedMeshComponent::check();

  if (getParam<std::vector<std::string>>("axial_region_names").size())
    checkEqualSize<std::string, Real>("axial_region_names", "length");
  else if (_n_sections > 1)
    logError("If there is more than 1 axial region, then the parameter 'axial_region_names' must "
             "be specified.");
}

bool
Component3D::hasRadialRegion(const std::string & name) const
{
  return std::find(_radial_region_names.begin(), _radial_region_names.end(), name) != _radial_region_names.end();
}

void
Component3D::build3DMesh()
{
  const auto node_ids = build3DMeshNodes();
  build3DMeshElems(node_ids);

  // // auto & boundary_info = mesh().getMesh().get_boundary_info();

  // // create elements from nodes
  // unsigned int i = 0;
  // for (unsigned int i_section = 0; i_section < _n_sections; i_section++)
  // {
  //   // element axial index for end of axial section
  //   unsigned int i_section_end = 0;
  //   for (unsigned int ii_section = 0; ii_section <= i_section; ++ii_section)
  //     i_section_end += _n_elems[ii_section];
  //   i_section_end -= 1;

  //   for (unsigned int i_local = 0; i_local < _n_elems[i_section]; i_local++)
  //   {
  //     unsigned int j = 0;
  //     for (unsigned int j_section = 0; j_section < _n_radial_regions; j_section++)
  //       for (unsigned int j_local = 0; j_local < _n_radial_elems[j_section]; j_local++)
  //       {
  //         for (unsigned int k = 0; k < _n_azimuthal_elems; k++)
  //         {
  //           Elem * elem = addElementPrism6(
  //               node_ids[i][j + 1][k], node_ids[i][j][k], node_ids[i + 1][j], node_ids[i + 1][j + 1]);
  //           elem->subdomain_id() = _subdomain_ids[j_section];
  //         }

  //         // exterior axial boundaries (all radial sections)
  //         // if (i == 0)
  //         // {
  //         //   boundary_info.add_side(elem, 0, _start_bc_id);
  //         //   _boundary_info[_boundary_name_start].push_back(
  //         //       std::tuple<dof_id_type, unsigned short int>(elem->id(), 0));
  //         // }
  //         // if (i == _n_elem - 1)
  //         // {
  //         //   boundary_info.add_side(elem, 2, _end_bc_id);
  //         //   _boundary_info[_boundary_name_end].push_back(
  //         //       std::tuple<dof_id_type, unsigned short int>(elem->id(), 2));
  //         // }

  //         // exterior axial boundaries (per radial section)
  //         // if (_radial_regions_names.size() > 1)
  //         // {
  //         //   if (i == 0)
  //         //   {
  //         //     boundary_info.add_side(elem, 0, _radial_start_bc_id[j_section]);
  //         //     _boundary_info[_boundary_names_radial_start[j_section]].push_back(
  //         //         std::tuple<dof_id_type, unsigned short int>(elem->id(), 0));
  //         //   }
  //         //   if (i == _n_elem - 1)
  //         //   {
  //         //     boundary_info.add_side(elem, 2, _radial_end_bc_id[j_section]);
  //         //     _boundary_info[_boundary_names_radial_end[j_section]].push_back(
  //         //         std::tuple<dof_id_type, unsigned short int>(elem->id(), 2));
  //         //   }
  //         // }

  //         // interior axial boundaries (per radial section)
  //         // if (_n_sections > 1 && _axial_region_names.size() == _n_sections &&
  //         //     i_section != _n_sections - 1 && i == i_section_end)
  //         // {
  //         //   const unsigned int k = i_section * _n_radial_regions + j_section;
  //         //   boundary_info.add_side(elem, 2, _interior_axial_per_radial_section_bc_id[k]);
  //         //   _boundary_info[_boundary_names_interior_axial_per_radial_section[k]].push_back(
  //         //       std::tuple<dof_id_type, unsigned short int>(elem->id(), 2));
  //         // }

  //         // exterior radial boundaries (all axial sections)
  //         // if (j == 0)
  //         // {
  //         //   boundary_info.add_side(elem, 1, _inner_bc_id);
  //         //   _boundary_info[_boundary_name_inner].push_back(
  //         //       std::tuple<dof_id_type, unsigned short int>(elem->id(), 1));
  //         // }
  //         // if (j == _total_radial_elems - 1)
  //         // {
  //         //   boundary_info.add_side(elem, 3, _outer_bc_id);
  //         //   _boundary_info[_boundary_name_outer].push_back(
  //         //       std::tuple<dof_id_type, unsigned short int>(elem->id(), 3));
  //         // }

  //         // exterior radial boundaries (per axial section)
  //         // if (_n_sections > 1 && _axial_region_names.size() == _n_sections)
  //         // {
  //         //   if (j == 0)
  //         //   {
  //         //     boundary_info.add_side(elem, 1, _axial_inner_bc_id[i_section]);
  //         //     _boundary_info[_boundary_names_axial_inner[i_section]].push_back(
  //         //         std::tuple<dof_id_type, unsigned short int>(elem->id(), 1));
  //         //   }
  //         //   if (j == _total_radial_elems - 1)
  //         //   {
  //         //     boundary_info.add_side(elem, 3, _axial_outer_bc_id[i_section]);
  //         //     _boundary_info[_boundary_names_axial_outer[i_section]].push_back(
  //         //         std::tuple<dof_id_type, unsigned short int>(elem->id(), 3));
  //         //   }
  //         // }

  //         // interior radial boundaries (all axial sections)
  //         // if (_n_radial_regions > 1 && _radial_regions_names.size() == _n_radial_regions && j_section != 0)
  //         // {
  //         //   unsigned int j_section_begin = 0;
  //         //   for (unsigned int jj_section = 0; jj_section < j_section; ++jj_section)
  //         //     j_section_begin += _n_radial_elems[jj_section];

  //         //   if (j == j_section_begin)
  //         //   {
  //         //     boundary_info.add_side(elem, 1, _inner_radial_bc_id[j_section - 1]);
  //         //     _boundary_info[_boundary_names_inner_radial[j_section - 1]].push_back(
  //         //         std::tuple<dof_id_type, unsigned short int>(elem->id(), 1));
  //         //   }
  //         // }

  //         j++;
  //       }

  //     i++;
  //   }
  // }
}

std::vector<std::vector<std::vector<unsigned int>>>
Component3D::build3DMeshNodes()
{
  unsigned int n_axial_positions = _node_locations.size();
  std::vector<std::vector<std::vector<unsigned int>>> node_ids(n_axial_positions, std::vector<std::vector<unsigned int>>(_total_radial_elems + 1, std::vector<unsigned int>(_n_azimuthal_elems)));

  // loop over axial positions
  for (unsigned int i = 0; i < n_axial_positions; i++)
  {
    // loop over radial positions
    for (unsigned int j = 0; j < _radial_positions.size(); j++)
    {
        const Real r = _radial_positions[j];

        // loop over the azimuthal elems
        for (unsigned int k = 0; k < _n_azimuthal_elems; k++)
        {
          const Real theta = k * _azimuthal_elem_width;
          const Point p(_node_locations[i], r * cos(theta), r * sin(theta));
          Node * nd = addNode(p);
          node_ids[i][j][k] = nd->id();
          std::cout<<"Adding "<<i<<","<<j<<","<<k<<":"<<p<<std::endl;
        }
      }
  }

  return node_ids;
}

void
Component3D::build3DMeshElems(const std::vector<std::vector<std::vector<unsigned int>>> & node_ids)
{
  unsigned int i = 0;
  for (unsigned int i_section = 0; i_section < _n_sections; i_section++)
  {
    for (unsigned int i_local = 0; i_local < _n_elems[i_section]; i_local++)
    {
      unsigned int j = 0;
      for (unsigned int j_section = 0; j_section < _n_radial_regions; j_section++)
        for (unsigned int j_local = 0; j_local < _n_radial_elems[j_section]; j_local++)
        {
          for (unsigned int k = 0; k < _n_azimuthal_elems; k++)
          {
            const auto k_next = k == _n_azimuthal_elems - 1 ? 0 : k+1;
            std::cout<<i<<","<<j<<","<<k<<","<<k_next<<std::endl;
            Elem * elem = addElementHex8(
                node_ids[i][j][k], node_ids[i+1][j][k], node_ids[i+1][j+1][k], node_ids[i][j+1][k], node_ids[i][j][k_next], node_ids[i+1][j][k_next], node_ids[i+1][j+1][k_next], node_ids[i][j+1][k_next]);
            elem->subdomain_id() = _subdomain_ids[j_section];
          }
          j++;
        }
      i++;
    }
  }
}

// void
// Component3D::build2DMesh2ndOrder()
// {
//   unsigned int n_axial_positions = _node_locations.size();
//   std::vector<std::vector<unsigned int>> node_ids(
//       n_axial_positions, std::vector<unsigned int>(2 * _total_radial_elems + 1));

//   // loop over axial positions
//   for (unsigned int i = 0; i < n_axial_positions; i++)
//   {
//     Point p(_node_locations[i], _axial_offset, 0);

//     const Node * nd = addNode(p);
//     node_ids[i][0] = nd->id();

//     // loop over regions
//     unsigned int l = 1;
//     for (unsigned int j = 0; j < _n_radial_regions; j++)
//     {
//       Real elem_length = _radial_region_widths[j] / (2. * _n_radial_elems[j]);
//       for (unsigned int k = 0; k < 2. * _n_radial_elems[j]; k++, l++)
//       {
//         p(1) += elem_length;
//         nd = addNode(p);
//         node_ids[i][l] = nd->id();
//       }
//     }
//   }

//   auto & boundary_info = mesh().getMesh().get_boundary_info();

//   // create elements from nodes
//   unsigned int i = 0;
//   for (unsigned int i_section = 0; i_section < _n_sections; i_section++)
//     for (unsigned int i_local = 0; i_local < _n_elems[i_section]; i_local++)
//     {
//       unsigned int j = 0;
//       for (unsigned int j_section = 0; j_section < _n_radial_regions; j_section++)
//         for (unsigned int j_local = 0; j_local < _n_radial_elems[j_section]; j_local++)
//         {
//           Elem * elem = addElementQuad9(node_ids[2 * i][2 * j],
//                                         node_ids[2 * i][2 * (j + 1)],
//                                         node_ids[2 * (i + 1)][2 * (j + 1)],
//                                         node_ids[2 * (i + 1)][2 * j],
//                                         node_ids[2 * i][(2 * j) + 1],
//                                         node_ids[(2 * i) + 1][2 * (j + 1)],
//                                         node_ids[2 * (i + 1)][(2 * j) + 1],
//                                         node_ids[(2 * i) + 1][(2 * j)],
//                                         node_ids[(2 * i) + 1][(2 * j) + 1]);
//           elem->subdomain_id() = _subdomain_ids[j_section];

//           if (i == 0)
//           {
//             boundary_info.add_side(elem, 0, _start_bc_id);
//             _boundary_info[_boundary_name_start].push_back(
//                 std::tuple<dof_id_type, unsigned short int>(elem->id(), 0));
//           }
//           if (i == _n_elem - 1)
//           {
//             boundary_info.add_side(elem, 2, _end_bc_id);
//             _boundary_info[_boundary_name_end].push_back(
//                 std::tuple<dof_id_type, unsigned short int>(elem->id(), 2));
//           }
//           if (_radial_regions_names.size() > 1)
//           {
//             if (i == 0)
//             {
//               boundary_info.add_side(elem, 0, _radial_start_bc_id[j_section]);
//               _boundary_info[_boundary_names_radial_start[j_section]].push_back(
//                   std::tuple<dof_id_type, unsigned short int>(elem->id(), 0));
//             }
//             if (i == _n_elem - 1)
//             {
//               boundary_info.add_side(elem, 2, _radial_end_bc_id[j_section]);
//               _boundary_info[_boundary_names_radial_end[j_section]].push_back(
//                   std::tuple<dof_id_type, unsigned short int>(elem->id(), 2));
//             }
//           }

//           if (j == 0)
//           {
//             boundary_info.add_side(elem, 3, _inner_bc_id);
//             _boundary_info[_boundary_name_inner].push_back(
//                 std::tuple<dof_id_type, unsigned short int>(elem->id(), 3));
//           }
//           if (j == _total_radial_elems - 1)
//           {
//             boundary_info.add_side(elem, 1, _outer_bc_id);
//             _boundary_info[_boundary_name_outer].push_back(
//                 std::tuple<dof_id_type, unsigned short int>(elem->id(), 1));
//           }

//           if (_n_sections > 1 && _axial_region_names.size() == _n_sections)
//           {
//             if (j == 0)
//             {
//               boundary_info.add_side(elem, 1, _axial_inner_bc_id[i_section]);
//               _boundary_info[_boundary_names_axial_inner[i_section]].push_back(
//                   std::tuple<dof_id_type, unsigned short int>(elem->id(), 1));
//             }
//             if (j == _total_radial_elems - 1)
//             {
//               boundary_info.add_side(elem, 3, _axial_outer_bc_id[i_section]);
//               _boundary_info[_boundary_names_axial_outer[i_section]].push_back(
//                   std::tuple<dof_id_type, unsigned short int>(elem->id(), 3));
//             }
//           }

//           // interior radial boundaries
//           if (_n_radial_regions > 1 && _radial_regions_names.size() == _n_radial_regions && j_section != 0)
//           {
//             unsigned int j_section_begin = 0;
//             for (unsigned int jj_section = 0; jj_section < j_section; ++jj_section)
//               j_section_begin += _n_radial_elems[jj_section];

//             if (j == j_section_begin)
//             {
//               boundary_info.add_side(elem, 1, _inner_radial_bc_id[j_section - 1]);
//               _boundary_info[_boundary_names_inner_radial[j_section - 1]].push_back(
//                   std::tuple<dof_id_type, unsigned short int>(elem->id(), 1));
//             }
//           }

//           j++;
//         }

//       i++;
//     }
// }

void
Component3D::buildMesh()
{
  if (_n_radial_elems.size() != _n_radial_regions || _radial_region_widths.size() != _n_radial_regions)
    return;

  // Assign subdomain to each transverse region
  for (unsigned int i = 0; i < _n_radial_regions; i++)
    setSubdomainInfo(mesh().getNextSubdomainId(), genName(_name, _radial_region_names[i]), Moose::COORD_XYZ);

  // // Create boundary IDs and associated boundary names
  // _inner_bc_id = mesh().getNextBoundaryId();
  // _outer_bc_id = mesh().getNextBoundaryId();
  // _boundary_name_inner = genName(name(), "inner");
  // _boundary_name_outer = genName(name(), "outer");
  // _boundary_name_to_area[_boundary_name_inner] = computeRadialBoundaryArea(_length, 0.0);
  // _boundary_name_to_area[_boundary_name_outer] =
  //     computeRadialBoundaryArea(_length, getTotalWidth());
  // if (_n_sections > 1 && _axial_region_names.size() == _n_sections)
  //   for (unsigned int i = 0; i < _n_sections; i++)
  //   {
  //     _axial_inner_bc_id.push_back(mesh().getNextBoundaryId());
  //     _axial_outer_bc_id.push_back(mesh().getNextBoundaryId());
  //     const BoundaryName boundary_name_axial_inner =
  //         genName(name(), _axial_region_names[i], "inner");
  //     const BoundaryName boundary_name_axial_outer =
  //         genName(name(), _axial_region_names[i], "outer");
  //     _boundary_names_axial_inner.push_back(boundary_name_axial_inner);
  //     _boundary_names_axial_outer.push_back(boundary_name_axial_outer);
  //     _boundary_name_to_area[boundary_name_axial_inner] =
  //         computeRadialBoundaryArea(_lengths[i], 0.0);
  //     _boundary_name_to_area[boundary_name_axial_outer] =
  //         computeRadialBoundaryArea(_lengths[i], getTotalWidth());
  //   }

  // // exterior axial boundaries
  // _start_bc_id = mesh().getNextBoundaryId();
  // _end_bc_id = mesh().getNextBoundaryId();
  // _boundary_name_start = genName(name(), "start");
  // _boundary_name_end = genName(name(), "end");
  // _boundary_name_to_area[_boundary_name_start] = computeAxialBoundaryArea(0.0, getTotalWidth());
  // _boundary_name_to_area[_boundary_name_end] = computeAxialBoundaryArea(0.0, getTotalWidth());
  // if (_radial_regions_names.size() > 1)
  // {
  //   Real y1 = 0.0;
  //   for (unsigned int i = 0; i < _radial_regions_names.size(); i++)
  //   {
  //     const Real y2 = y1 + _radial_region_widths[i];

  //     _radial_start_bc_id.push_back(mesh().getNextBoundaryId());
  //     _radial_end_bc_id.push_back(mesh().getNextBoundaryId());
  //     const BoundaryName boundary_name_radial_start = genName(name(), _radial_regions_names[i], "start");
  //     const BoundaryName boundary_name_radial_end = genName(name(), _radial_regions_names[i], "end");
  //     _boundary_names_radial_start.push_back(boundary_name_radial_start);
  //     _boundary_names_radial_end.push_back(boundary_name_radial_end);
  //     _boundary_name_to_area[boundary_name_radial_start] = computeAxialBoundaryArea(y1, y2);
  //     _boundary_name_to_area[boundary_name_radial_end] = computeAxialBoundaryArea(y1, y2);
  //     if (i != _radial_regions_names.size() - 1)
  //     {
  //       _inner_radial_bc_id.push_back(mesh().getNextBoundaryId());
  //       const BoundaryName boundary_name_inner_radial = genName(name(), _radial_regions_names[i], _radial_regions_names[i + 1]);
  //       _boundary_names_inner_radial.push_back(boundary_name_inner_radial);
  //       _boundary_name_to_area[boundary_name_inner_radial] = computeRadialBoundaryArea(_length, y2);
  //     }
  //     y1 = y2;
  //   }
  // }

  // // interior axial boundaries
  // if (_n_sections > 1 && _axial_region_names.size() == _n_sections)
  //   for (unsigned int i = 0; i < _n_sections - 1; i++)
  //   {
  //     Real y1 = 0.0;
  //     for (unsigned int j = 0; j < _radial_regions_names.size(); j++)
  //     {
  //       const Real y2 = y1 + _radial_region_widths[j];

  //       _interior_axial_per_radial_section_bc_id.push_back(mesh().getNextBoundaryId());
  //       const BoundaryName boundary_name_interior_axial_per_radial_section =
  //           genName(name(), _radial_regions_names[j], _axial_region_names[i] + ":" + _axial_region_names[i + 1]);
  //       _boundary_names_interior_axial_per_radial_section.push_back(
  //           boundary_name_interior_axial_per_radial_section);
  //       _boundary_name_to_area[boundary_name_interior_axial_per_radial_section] =
  //           computeAxialBoundaryArea(y1, y2);
  //       y1 = y2;
  //     }
  //   }

  if (usingSecondOrderMesh())
    mooseError("Not implemented");

  // Build the mesh
  build3DMesh();

  // // Set boundary names
  // auto & binfo = mesh().getMesh().get_boundary_info();
  // binfo.sideset_name(_inner_bc_id) = _boundary_name_inner;
  // binfo.sideset_name(_outer_bc_id) = _boundary_name_outer;
  // if (_n_sections > 1 && _axial_region_names.size() == _n_sections)
  //   for (unsigned int i = 0; i < _n_sections; i++)
  //   {
  //     binfo.sideset_name(_axial_inner_bc_id[i]) = _boundary_names_axial_inner[i];
  //     binfo.sideset_name(_axial_outer_bc_id[i]) = _boundary_names_axial_outer[i];
  //   }
  // binfo.sideset_name(_start_bc_id) = _boundary_name_start;
  // binfo.sideset_name(_end_bc_id) = _boundary_name_end;
  // if (_radial_regions_names.size() > 1)
  //   for (unsigned int i = 0; i < _radial_regions_names.size(); i++)
  //   {
  //     binfo.sideset_name(_radial_start_bc_id[i]) = _boundary_names_radial_start[i];
  //     binfo.sideset_name(_radial_end_bc_id[i]) = _boundary_names_radial_end[i];
  //     if (i != _radial_regions_names.size() - 1)
  //       binfo.sideset_name(_inner_radial_bc_id[i]) = _boundary_names_inner_radial[i];
  //   }
  // for (unsigned int k = 0; k < _interior_axial_per_radial_section_bc_id.size(); k++)
  //   binfo.sideset_name(_interior_axial_per_radial_section_bc_id[k]) =
  //       _boundary_names_interior_axial_per_radial_section[k];
}

// bool
// Component3D::isBoundaryInVector(const BoundaryName & boundary_name,
//                                 const std::vector<BoundaryName> & boundary_name_vector) const
// {
//   return std::find(boundary_name_vector.begin(), boundary_name_vector.end(), boundary_name) !=
//          boundary_name_vector.end();
// }

// bool
// Component3D::hasBoundary(const BoundaryName & boundary_name) const
// {
//   checkSetupStatus(MESH_PREPARED);

//   return hasExternalBoundary(boundary_name) ||
//          isBoundaryInVector(boundary_name, _boundary_names_interior_axial_per_radial_section) ||
//          isBoundaryInVector(boundary_name, _boundary_names_inner_radial);
// }

// bool
// Component3D::hasExternalBoundary(const BoundaryName & boundary_name) const
// {
//   checkSetupStatus(MESH_PREPARED);

//   return boundary_name == _boundary_name_inner || boundary_name == _boundary_name_outer ||
//          boundary_name == _boundary_name_start || boundary_name == _boundary_name_end ||
//          isBoundaryInVector(boundary_name, _boundary_names_axial_inner) ||
//          isBoundaryInVector(boundary_name, _boundary_names_axial_outer) ||
//          isBoundaryInVector(boundary_name, _boundary_names_radial_start) ||
//          isBoundaryInVector(boundary_name, _boundary_names_radial_end);
// }

// Component3D::ExternalBoundaryType
// Component3D::getExternalBoundaryType(const BoundaryName & boundary_name) const
// {
//   checkSetupStatus(MESH_PREPARED);

//   if (boundary_name == _boundary_name_inner ||
//       isBoundaryInVector(boundary_name, _boundary_names_axial_inner))
//     return ExternalBoundaryType::INNER;
//   else if (boundary_name == _boundary_name_outer ||
//            isBoundaryInVector(boundary_name, _boundary_names_axial_outer))
//     return ExternalBoundaryType::OUTER;
//   else if (boundary_name == _boundary_name_start ||
//            isBoundaryInVector(boundary_name, _boundary_names_radial_start))
//     return ExternalBoundaryType::START;
//   else if (boundary_name == _boundary_name_end ||
//            isBoundaryInVector(boundary_name, _boundary_names_radial_end))
//     return ExternalBoundaryType::END;
//   else if (hasBoundary(boundary_name))
//     mooseError(name(), ": The boundary '", boundary_name, "' is an interior boundary.");
//   else
//     mooseError(name(), ": The boundary '", boundary_name, "' does not exist on this component.");
// }

// const std::vector<std::tuple<dof_id_type, unsigned short int>> &
// Component3D::getBoundaryInfo(const BoundaryName & boundary_name) const
// {
//   checkSetupStatus(MESH_PREPARED);

//   if (_boundary_info.find(boundary_name) != _boundary_info.end())
//     return _boundary_info.at(boundary_name);
//   else
//     mooseError(name(), ": The boundary '", boundary_name, "' does not exist on this component.");
// }

// const std::vector<std::tuple<dof_id_type, unsigned short int>> &
// Component3D::getBoundaryInfo(const ExternalBoundaryType & boundary_type) const
// {
//   checkSetupStatus(MESH_PREPARED);

//   switch (boundary_type)
//   {
//     case ExternalBoundaryType::INNER:
//       return getBoundaryInfo(_boundary_name_inner);
//     case ExternalBoundaryType::OUTER:
//       return getBoundaryInfo(_boundary_name_outer);
//     case ExternalBoundaryType::START:
//       return getBoundaryInfo(_boundary_name_start);
//     case ExternalBoundaryType::END:
//       return getBoundaryInfo(_boundary_name_end);
//     default:
//       mooseError(name(), ": Invalid external boundary type.");
//   }
// }

// const BoundaryName &
// Component3D::getExternalBoundaryName(const ExternalBoundaryType & boundary_type) const
// {
//   checkSetupStatus(MESH_PREPARED);

//   switch (boundary_type)
//   {
//     case ExternalBoundaryType::OUTER:
//       return _boundary_name_outer;
//     case ExternalBoundaryType::INNER:
//       return _boundary_name_inner;
//     case ExternalBoundaryType::START:
//       return _boundary_name_start;
//     case ExternalBoundaryType::END:
//       return _boundary_name_end;
//     default:
//       mooseError(name(), ": Invalid external boundary type.");
//   }
// }

// const Real &
// Component3D::getBoundaryArea(const BoundaryName & boundary_name) const
// {
//   checkSetupStatus(MESH_PREPARED);

//   if (_boundary_name_to_area.find(boundary_name) != _boundary_name_to_area.end())
//     return _boundary_name_to_area.at(boundary_name);
//   else
//     mooseError(name(), ": The boundary '", boundary_name, "' does not exist on this component.");
// }
