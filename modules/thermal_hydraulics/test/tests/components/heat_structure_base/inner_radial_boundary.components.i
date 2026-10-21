[hs]
  type = HeatStructureCylindrical

  position = '0 0 0'
  orientation = '1 0 0'
  length = 10.0
  n_elems = 20

  names = 'region1 region2 region3'
  widths = '1.0 3.0 2.0'
  n_part_elems = '2 6 8'
  solid_properties = 'hs_mat hs_mat hs_mat'
  solid_properties_T_ref = '300 300 300'

  initial_T = initial_T_fn
[]
