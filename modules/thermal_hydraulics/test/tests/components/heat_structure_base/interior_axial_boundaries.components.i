[hs]
  type = HeatStructurePlate

  position = '0 0 0'
  orientation = '1 0 0'
  length = '20 20 20'
  n_elems = '2 2 2'
  axial_region_names = 'axial1 axial2 axial3'

  names = 'radial1 radial2 radial3'
  widths = '0.5 0.5 0.5'
  n_part_elems = '2 2 2'
  solid_properties = 'hs_mat hs_mat hs_mat'
  solid_properties_T_ref = '300 300 300'

  depth = 1.0

  initial_T = initial_T_fn
[]
