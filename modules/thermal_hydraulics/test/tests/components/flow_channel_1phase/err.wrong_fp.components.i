[pipe]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  A = 1.
  D_h = 1.12837916709551
  f = 0.01
  length = 1
  n_elems = 100
  fp = fp_2phase   # this is wrong
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
  p = 9.5e4
[]
