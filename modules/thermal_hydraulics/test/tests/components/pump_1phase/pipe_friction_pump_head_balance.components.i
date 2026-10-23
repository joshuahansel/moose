[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 10
  initial_p = 1.35e+07
  n_elems = 20
  f = 5000
  gravity_vector = '0 0 0'
[]
[pump]
  type = Pump1Phase
  connections = 'pipe1:out pipe1:in'
  position = '1.02 0 0'
  initial_p = 1.3e+07
  initial_vel_x = 1
  initial_vel_y = 0
  initial_vel_z = 0
  scaling_factor_rhoV  = 1
  scaling_factor_rhouV = 1
  scaling_factor_rhoEV = 1e-5
  head = 8
  volume = 0.567
  A_ref = 0.567
[]
