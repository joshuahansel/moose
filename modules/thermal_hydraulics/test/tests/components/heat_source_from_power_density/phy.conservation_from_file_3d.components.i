[heat_structure]
  type = HeatStructureFromFile3D
  file = box.e
  position = '0 0 0'
  initial_T = 300
[]
[heat_generation]
  type = HeatSourceFromPowerDensity
  hs = heat_structure
  regions = 'rgn1'
  power_density = power_density
[]
