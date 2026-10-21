[total_power]
  type = TotalPower
  power = 3.0e4
[]

[CH1:solid]
  type = HeatStructureCylindrical
  position = '0 -0.024 0'
  orientation = '0 0 1'
  length = 0.8
  n_elems = 16

  initial_T = 628.15

  names = 'fuel gap clad'
  widths = '0.003015 0.000465  0.00052'
  n_part_elems = '20 2 2'
  solid_properties = 'fuel-mat gap-mat clad-mat'
  solid_properties_T_ref = '300 300 300'
[]

[CH1:hgen]
  type = HeatSourceFromPowerDensity
  hs = CH1:solid
  regions = 'fuel'
  power_density = power_density
[]
