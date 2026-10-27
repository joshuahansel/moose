[pipe1]
  type = FlowChannel1Phase
  fp = water
  position = '0 0 0'
  orientation = '0 1 0'
  length = 1
  n_elems = 2
  A = 1e-4
  f = 0
[]
[inlet_1p]
  type = InletMassFlowRateTemperature1Phase
  input = 'pipe:in'
  m_dot = 1
  T = 300
[]
[outlet_1p]
  type = Outlet1Phase
  input = 'pipe1:out'
  p = 1e5
[]
