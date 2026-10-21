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
#include "NamingInterface.h"
#include "HeatStructureInterfaceAC.h"

/**
 * Base class for 2D heat-conducting solid ActionComponents (see HeatStructureCylindricalAC/
 * HeatStructurePlateAC), with a mesh built through the MeshGenerator system rather than classic
 * THM's hand-built Component2D mesh.
 *
 * The mesh is divided along the transverse ('y') direction into named regions (materials), each
 * built as its own GeneratedMeshGenerator, then combined with StackGenerator (which stitches
 * consecutive regions' shared top/bottom boundaries). Boundaries are named to match classic THM's
 * convention: '<name>:inner'/'<name>:outer' (aggregate transverse faces), '<name>:start'/
 * '<name>:end' (aggregate axial faces), '<name>:<region>:start'/'<name>:<region>:end' (per-region
 * axial faces), and '<name>:<region_i>:<region_j>' (interior region interfaces). Unlike classic
 * Component2D, each region is made of one or more axial sections ('length'/'n_elems' list one entry
 * per section). If there is more than one section, the sections are named by 'axial_region_names'
 * and additional boundaries are created, as in classic THM: '<name>:<axial_region>:inner'/
 * '<name>:<axial_region>:outer' and '<name>:<region>:<axial_region_i>:<axial_region_j>' (the
 * interior axial interfaces of each radial region).
 *
 * The heat conduction physics is created by HeatStructureInterfaceAC.
 */
class HeatStructureBaseAC : public virtual ActionComponent,
                            public HeatStructureInterfaceAC,
                            public ComponentMeshTransformHelper,
                            public NamingInterface
{
public:
  static InputParameters validParams();
  HeatStructureBaseAC(const InputParameters & params);

  /// Type of an external boundary, mirroring classic THM's Component2D::ExternalBoundaryType
  enum class ExternalBoundaryType
  {
    INNER,
    OUTER,
    START,
    END
  };

  /// Number of rods (units) represented by this heat structure
  virtual Real numberOfUnits() const override { return _num_rods; }
  /// Total length of the heat structure (summed over its axial sections)
  virtual Real length() const override { return _length; }
  /// Whether \p region is one of this heat structure's transverse regions
  virtual bool hasRegion(const std::string & region) const override;

  /// Whether \p boundary is an external boundary of this heat structure: the aggregate
  /// inner/outer/start/end boundaries, a per-region start/end boundary, or a per-axial-region
  /// inner/outer boundary
  bool hasExternalBoundary(const BoundaryName & boundary) const;
  /// Type of the external boundary \p boundary; errors if it is not an external boundary
  ExternalBoundaryType getExternalBoundaryType(const BoundaryName & boundary) const;
  /// Perimeter of one unit of the heat structure on the given side (NaN for the start/end sides,
  /// which have no well-defined perimeter), mirroring HeatStructureBase::getUnitPerimeter()
  virtual Real getUnitPerimeter(ExternalBoundaryType side) const = 0;

protected:
  virtual void addMeshGenerators() override;
  virtual void addMaterials() override;
  virtual void checkIntegrity() override;

  /// Returns the transverse ('y') offset of the innermost region's inner face from the component's
  /// own axis - the inner radius for a cylindrical heat structure, or 0 for a plate
  virtual Real yOffset() const = 0;

  /// Adds a single region's material properties from the 'solid_properties'/'solid_properties_T_ref'
  /// convenience parameters
  void addConstantDensitySolidPropertiesMaterial(const UserObjectName & sp_name,
                                                 Real T_ref,
                                                 const SubdomainName & block);

  /// Names of the transverse regions
  std::vector<std::string> _region_names;
  /// Width of each transverse region
  std::vector<Real> _widths;
  /// Number of elements of each transverse region
  std::vector<unsigned int> _n_part_elems;

  /// Number of rods (units) represented by this heat structure
  const Real _num_rods;

  /// Length of each axial section
  const std::vector<Real> _lengths;
  /// Number of axial elements of each axial section
  const std::vector<unsigned int> _n_elems;
  /// Names of the axial sections (empty if there is a single, unnamed section)
  const std::vector<std::string> _axial_region_names;
  /// Total length of the heat structure
  const Real _length;
  /// Whether 2nd order elements are used
  const bool _second_order;
  /// Whether the axial sections are named, so that per-section boundaries exist
  const bool _has_axial_regions;
};
