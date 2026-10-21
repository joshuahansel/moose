[hs3d]
  type = HeatStructureFromFile3D
  file = mesh/mesh.e
  position = '0 0 0'
  initial_T = initial_T_matrix_fn
[]

[hs2d]
  type = HeatStructureCylindrical
  orientation = '0 0 1'
  position = '0 0 0'
  length = '${length_matrix} ${length_extend}'
  n_elems = '${n_elems_axial_matrix} ${n_elems_axial_extend}'
  axial_region_names = 'matrix extend'

  inner_radius = 0
  widths = '${R_pipe}'
  n_part_elems = '${n_elems_radial}'
  names = 'pipe'

  initial_T = 300
[]

[hs_coupler]
  type = HSCoupler2D3D
  heat_structure_2d = hs2d
  heat_structure_3d = hs3d
  boundary_2d = hs2d:matrix:outer
  boundary_3d = hs3d:rmin

  include_radiation = false
  gap_thickness = 0.00001
  gap_thermal_conductivity = 0.05
[]
