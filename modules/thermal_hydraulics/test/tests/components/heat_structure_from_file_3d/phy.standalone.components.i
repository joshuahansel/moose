[blk]
  type = HeatStructureFromFile3D
  file = box.e
  position = '0 0 0'
  initial_T = 300
[]

[left_bnd]
  type = HSBoundarySpecifiedTemperature
  hs = blk
  boundary = blk:left
  T = 300
[]

[right_bnd]
  type = HSBoundarySpecifiedTemperature
  hs = blk
  boundary = blk:right
  T = 310
[]
