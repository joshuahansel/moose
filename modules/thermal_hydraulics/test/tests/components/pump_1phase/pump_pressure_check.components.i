[inlet]
  type = InletMassFlowRateTemperature1Phase
  input = 'pipe1:in'
  m_dot = 20
  T = 393.15
[]

[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  initial_p = 1.318964e+07
  n_elems = 10
[]

[pump]
  type = Pump1Phase
  connections = 'pipe1:out pipe2:in'
  position = '1.02 0 0'
  initial_p = 1.318964e+07
  scaling_factor_rhoEV = 1e-5
  head = ${head}
  volume = ${volume}
  A_ref = 0.567
  initial_vel_x = 1
  initial_vel_y = 1
  initial_vel_z = 0
[]

[pipe2]
  type = FlowChannel1Phase
  position = '1.04 0 0'
  orientation = '0 2 0'
  length = 0.96
  initial_p = 1.4072E+07
  n_elems = 10
[]

[outlet]
  type = Outlet1Phase
  input = 'pipe2:out'
  p = 1.4072E+07
[]
