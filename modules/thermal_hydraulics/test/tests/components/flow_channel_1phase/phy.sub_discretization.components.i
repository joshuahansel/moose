[pipe]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  axial_region_names = 'r1 r2'
  length = '1 1'
  n_elems = '1 2'
  A = 1
  f = 0
  fp = eos
[]

[inlet]
  type = SolidWall1Phase
  input = 'pipe:in'
[]

[outlet]
  type = SolidWall1Phase
  input = 'pipe:out'
[]
