[pipe]
  type = FlowChannel1Phase
  gravity_vector = '0 0 0'
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1e4
  n_elems = 10
  A = 1e3

  initial_T = 300
  initial_p = 1e5
  initial_vel = 0

  fp = fp
  closures = simple_closures
  f = 0

  scaling_factor_1phase = '1 1 1'
[]

[inlet]
  type = InletStagnationPressureTemperature1Phase
  input = 'pipe:in'
  p0 = 1e5
  T0 = 300
[]

[outlet]
  type = Outlet1Phase
  input = 'pipe:out'
  p = 1.1e5
[]
