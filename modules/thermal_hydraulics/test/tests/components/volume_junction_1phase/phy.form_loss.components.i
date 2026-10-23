[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  A = ${A}
[]

[pipe2]
  type = FlowChannel1Phase
  position = '1 0 0'
  orientation = '1 0 0'
  A = ${A}
  initial_p = 1e5
[]

[junction]
  type = VolumeJunction1Phase
  connections = 'pipe1:out pipe2:in'

  position = '1 0 0'
  volume = 0.005
  initial_p = 1e5
  initial_vel_x = 1
  initial_vel_y = 0
  initial_vel_z = 0
[]

[pipe1_in]
  type = InletVelocityTemperature1Phase
  input = 'pipe1:in'
  vel = 1
  T = 300
[]

[pipe2_out]
  type = Outlet1Phase
  input = 'pipe2:out'
  p = 1e5
[]
