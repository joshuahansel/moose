[reactor]
  type = TotalPower
  power = 296153.84615384615385
[]

[hs]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '1 0 0'

  length = 1
  n_elems = 1
  names = 'FUEL GAP CLAD'
  widths = '0.0046955  0.0000955  0.000673'
  n_part_elems = '1 1 1'
  solid_properties = 'fuel-mat gap-mat clad-mat'
  solid_properties_T_ref = '300 300 300'
[]

[temp_outside]
  type = HSBoundarySpecifiedTemperature
  hs = hs
  boundary = hs:outer
  T = 600
[]
