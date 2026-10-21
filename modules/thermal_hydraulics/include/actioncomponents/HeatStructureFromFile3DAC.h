//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HeatStructureInterfaceAC.h"
#include "ComponentMeshTransformHelper.h"
#include "NamingInterface.h"

/**
 * Heat structure ActionComponent whose 3D mesh is loaded from an ExodusII file.
 *
 * The mesh is read with a FileMeshGenerator. To match the naming convention of the classic
 * HeatStructureFromFile3D, every block and sideset of the file is renamed to '<name>:<file name>'
 * (the block or sideset ID is used when the file gives it no name). The materials of the regions
 * are not created by this component and must be supplied by the user.
 */
class HeatStructureFromFile3DAC : public virtual ActionComponent,
                                  public HeatStructureInterfaceAC,
                                  public ComponentMeshTransformHelper,
                                  public NamingInterface
{
public:
  static InputParameters validParams();
  HeatStructureFromFile3DAC(const InputParameters & params);

  virtual bool hasRegion(const std::string & region) const override;
  /// Whether \p boundary is one of the sidesets of the mesh file (by its name in the simulation)
  bool hasBoundary(const BoundaryName & boundary) const;
  /// A mesh from a file represents a single unit
  virtual Real numberOfUnits() const override { return 1.0; }
  /// Without a notion of length, a unit length is used to normalize a uniform power shape
  virtual Real length() const override { return 1.0; }

protected:
  virtual void addMeshGenerators() override;

  /// Names of the file's blocks (regions)
  std::vector<std::string> _region_names;
  /// Names of the file's sidesets
  std::vector<std::string> _file_boundary_names;
};
