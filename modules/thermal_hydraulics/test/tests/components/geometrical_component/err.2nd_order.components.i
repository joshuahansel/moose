[hs]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '1 0 0'

  length = 1
  n_elems = 2
  names = 'blk'
  widths = '1'
  n_part_elems = '2'
  solid_properties = 'hs-mat'
  solid_properties_T_ref = '300'

  initial_T = 350
[]

[start]
  type = HSBoundarySpecifiedTemperature
  hs = hs
  boundary = hs:start
  T = 300
[]
