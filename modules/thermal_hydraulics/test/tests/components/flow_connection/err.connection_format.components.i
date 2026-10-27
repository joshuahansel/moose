[left_wall]
  type = SolidWall1Phase
[]

[pipe]
  type = FlowChannel1Phase
  fp = fp
  closures = simple_closures

  position = '0 0 0'
  orientation = '1 0 0'
  length = 1.0
  n_elems = 5
  A = 1.0

  initial_T = 300
  initial_p = 1e5
  initial_vel = 0

  f = 0
[]

[right_wall]
  type = SolidWall1Phase
  input = 'pipe:out'
[]
