[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 2
  gravity_vector = '0 0 0'
[]

[pump]
  type = Pump1Phase
  connections = 'pipe1:out pipe2:in'
  position = '1.02 0 0'
  head = 95
  A_ref = 1
  volume = 1
  initial_vel_x = 0
  initial_vel_y = 0
  initial_vel_z = 0
[]

[pipe2]
  type = FlowChannel1Phase
  position = '1.04 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 2
[]
