[heat_structure]
  type = HeatStructureCylindrical
  num_rods = ${n_units}

  position = '0 1 0'
  orientation = '1 0 0'
  length = 0.8
  n_elems = 100

  names = 'rgn1 rgn2 rgn3'
  solid_properties = 'main-material main-material main-material'
  solid_properties_T_ref = '300 300 300'
  widths = '0.4 0.1 0.5'
  n_part_elems = '2 2 2'

  initial_T = 300
[]
[heat_generation]
  type = HeatSourceFromTotalPower
  hs = heat_structure
  regions = 'rgn1 rgn2'
  power = total_power
  power_fraction = ${power_fraction}
[]
[total_power]
  type = TotalPower
  power = ${power}
[]
