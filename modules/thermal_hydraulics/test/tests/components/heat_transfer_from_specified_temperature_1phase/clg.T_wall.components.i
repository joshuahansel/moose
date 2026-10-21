[pipe1]
  type = FlowChannel1Phase
  fp = eos
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1.0
  n_elems = 50

  A = 3.14e-2
  f = 0.1
[]

[ht_pipe1]
  type = HeatTransferFromSpecifiedTemperature1Phase
  flow_channel = pipe1
  T_wall = T_wall_fn
  Hw = 0
[]

[inlet1]
  type = InletDensityVelocity1Phase
  input = 'pipe1:in'
  rho = 996.557482499661660
  vel = 1
[]

[outlet]
  type = Outlet1Phase
  input = 'pipe1:out'
  p = 0.1e6
[]
