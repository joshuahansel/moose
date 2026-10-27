[pipe]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  gravity_vector = '0 0 0'
  length = 1.0
  n_elems = 10
  A = 0.1

  initial_T = T_ic_fn
  initial_p = 1e5
  initial_vel = 0

  enable_heat_conduction = true

  fp = fp
  closures = simple_closures
  f = 0

  scaling_factor_1phase = '1 1 1e-5'
[]
[left_wall]
  type = SolidWall1Phase
  input = 'pipe:in'
[]
[right_wall]
  type = SolidWall1Phase
  input = 'pipe:out'
[]
