[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  A = ${A_big}
  # This pipe is pressurized higher than the others.
  initial_p = 1.05e5
[]

[pipe2]
  type = FlowChannel1Phase
  position = '1 0 0'
  orientation = '1 0 0'
  A = ${A_big}
  initial_p = 1e5
[]

[pipe3]
  type = FlowChannel1Phase
  position = '1 0 0'
  orientation = '0 1 0'
  # This pipe is smaller than the others.
  A = ${A_small}
  initial_p = 1e5
[]

[junction]
  type = VolumeJunction1Phase
  connections = 'pipe1:out pipe2:in pipe3:in'

  position = '1 0 0'
  volume = 0.37

  initial_p = 1e5
  initial_vel_x = 0
  initial_vel_y = 0
  initial_vel_z = 0
[]

[pipe1_wall]
  type = SolidWall1Phase
  input = 'pipe1:in'
[]
[pipe2_wall]
  type = SolidWall1Phase
  input = 'pipe2:out'
[]
[pipe3_wall]
  type = SolidWall1Phase
  input = 'pipe3:out'
[]
