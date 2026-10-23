[inlet]
  type = InletMassFlowRateTemperature1Phase
  input = 'pipe1:in'
  m_dot = 10
  T = 250
[]

[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  A = 1
  n_elems = 20
  initial_vel = 20
[]

[junction]
  type = JunctionParallelChannels1Phase
  connections = 'pipe1:out pipe2:in'
  scaling_factor_rhouV = 1e-4
  scaling_factor_rhoEV = 1e-5
  position = '1 0 0'
  volume = 1e-8
[]

[pipe2]
  type = FlowChannel1Phase
  position = '1 0 0'
  orientation = '1 0 0'
  length = 1
  A = 0.5
  n_elems = 20
  initial_vel = 15
[]

[outlet]
  type = Outlet1Phase
  input = 'pipe2:out'
  p = 1e5
[]
