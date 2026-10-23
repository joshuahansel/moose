[inlet]
  type = InletMassFlowRateTemperature1Phase
  input = 'pipe1:in'
  m_dot = 10
  T = 517
[]

[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 10
  A = 1
[]

[turbine]
  type = SimpleTurbine1Phase
  connections = 'pipe1:out pipe2:in'
  position = '1 0 0'
  volume = 1
  on = true
  power = 1000
[]

[pipe2]
  type = FlowChannel1Phase
  position = '1. 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 10
  A = 1
[]

[outlet]
  type = Outlet1Phase
  input = 'pipe2:out'
  p = 1e6
[]
