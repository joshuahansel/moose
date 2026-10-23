[wall_in]
  type = SolidWall1Phase
  input = 'pipe1:in'
[]

[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  initial_p = 1.7E+07
  n_elems = 10
  gravity_vector = '0 0 0'
[]

[pump]
  type = Pump1Phase
  connections = 'pipe1:out pipe2:in'
  position = '1.02 0 0'
  initial_p = 1.3e+07
  scaling_factor_rhoEV = 1e-5
  head = ${head}
  A_ref = ${A}
  volume = ${volume}
  initial_vel_x = 0
  initial_vel_y = 0
  initial_vel_z = 0
[]

[pipe2]
  type = FlowChannel1Phase
  position = '1.04 0 0'
  orientation = '1 0 0'
  length = 1
  initial_p = 1.3e+07
  n_elems = 10
  gravity_vector = '0 0 0'
[]

[wall_out]
  type = SolidWall1Phase
  input = 'pipe2:out'
[]
