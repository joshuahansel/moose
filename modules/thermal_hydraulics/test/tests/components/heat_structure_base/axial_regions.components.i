[hs]
  type = HeatStructurePlate
  position = '1 2 3'
  orientation = '0 0 1'
  depth = 1

  length = '3 2 1'
  n_elems = '2 4 3'

  names = 'radialregion'
  widths = '0.5'
  n_part_elems = '3'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '300'

  initial_T = 300
[]

[hs_boundary]
  type = HSBoundaryHeatFlux
  boundary = 'hs:region1:outer hs:region2:outer'
  hs = hs
  q = 1000
[]
