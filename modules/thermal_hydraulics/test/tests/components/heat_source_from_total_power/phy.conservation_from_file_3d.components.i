[heat_structure]
  type = HeatStructureFromFile3D
  file = box.e
  position = '0 0 0'
  initial_T = 300
[]
[heat_generation]
  type = HeatSourceFromTotalPower
  hs = heat_structure
  regions = 'rgn1'
  power = total_power
  power_fraction = ${power_fraction}
[]
[total_power]
  type = TotalPower
  power = ${power}
[]
