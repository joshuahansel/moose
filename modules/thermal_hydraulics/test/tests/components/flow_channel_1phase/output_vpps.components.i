[inlet]
  type = InletStagnationPressureTemperature1Phase
  input = 'pipe:in'
  p0 = 1e5
  T0 = 300
[]
[pipe]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1.0
  n_elems = 5
  A = 1.0
  initial_T = 300
  initial_p = 1e5
  initial_vel = 0
  fp = fp
  closures = simple_closures
  f = 0
  vpp_vars = 'p T'
  create_flux_vpp = true
[]
[outlet]
  type = Outlet1Phase
  input = 'pipe:out'
  p = 0.9e5
[]
