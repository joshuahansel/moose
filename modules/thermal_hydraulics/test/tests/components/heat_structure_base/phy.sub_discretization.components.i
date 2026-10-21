[hs]
  type = HeatStructureCylindrical
  position = '0 0 1'
  orientation = '1 0 0'

  axial_region_names = 'reg1 reg2'
  length = '2.0 1.6576'
  n_elems = '7   4'

  names = 'FUEL GAP CLAD'
  widths = '0.0046955  0.0000955  0.000673'
  n_part_elems = '10 3 3'
  solid_properties = 'fuel-mat gap-mat clad-mat'
  solid_properties_T_ref = '300 300 300'

  initial_T = 300
[]

[temp_outside]
  type = HSBoundarySpecifiedTemperature
  hs = hs
  boundary = hs:outer
  T = 300
[]
