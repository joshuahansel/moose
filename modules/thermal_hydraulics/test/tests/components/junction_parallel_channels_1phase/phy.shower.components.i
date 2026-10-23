[inlet_hot]
  type = InletDensityVelocity1Phase
  input = 'pipe_hot:in'
  # rho @ (p = 1e5, T = 310 K)
  rho = 1315.9279785683
  vel = 1
[]
[inlet_cold]
  type = InletDensityVelocity1Phase
  input = 'pipe_cold:in'
  # rho @ (p = 1e5, T = 280 K)
  rho = 1456.9202619863
  vel = 1
[]
[outlet]
  type = Outlet1Phase
  input = 'pipe_warm:out'
  p = 1e5
[]

[pipe_hot]
  type = FlowChannel1Phase
  position = '0 1 0'
  A = 1
[]
[pipe_cold]
  type = FlowChannel1Phase
  position = '0 0 0'
  A = 1
[]
[pipe_warm]
  type = FlowChannel1Phase
  position = '1 0.5 0'
  A = 2
  initial_vel = 0.5
[]

[junction]
  type = JunctionParallelChannels1Phase
  connections = 'pipe_cold:out pipe_hot:out pipe_warm:in'
  position = '1 0.5 0'
  volume = 1e-8
[]
