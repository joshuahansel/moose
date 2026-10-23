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
  n_elems = 10
  A = 0.567

  initial_p = 1.318964e+07
[]

[pump]
  type = Pump1Phase
  connections = 'pipe1:out pipe2:in'
  position = '1.02 0 0'
  head = 0
  volume = 0.567
  A_ref = 0.567

  initial_p = 1.318964e+07
  initial_vel_x = 0.0372
  initial_vel_y = 0
  initial_vel_z = 0
  scaling_factor_rhoV  = 1
  scaling_factor_rhouV = 1
  scaling_factor_rhoEV = 1e-5
[]

[pipe2]
  type = FlowChannel1Phase
  position = '1.04 0 0'
  orientation = '1 0 0'
  length = 0.96
  n_elems = 10
  A = 0.567

  initial_p = 1.4072e+07
[]

[outlet]
  type = Outlet1Phase
  input = 'pipe2:out'
  p = 1.4072e+07
[]
