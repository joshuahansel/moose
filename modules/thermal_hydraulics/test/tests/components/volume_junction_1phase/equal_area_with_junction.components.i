[inlet]
  type = InletStagnationPressureTemperature1Phase
  input = 'pipe1:in'
  # Stagnation property with p = 1e5 Pa, T = 250 K, vel = 1 m/s
  p0 = 100000.68965687
  T0 = 250.00049261084
[]

[pipe1]
  type = FlowChannel1Phase

  position = '0 0 0'
  orientation = '1 0 0'
  length = 1

  initial_T = T0

  n_elems = 25
[]

[junction]
  type = VolumeJunction1Phase
  connections = 'pipe1:out pipe2:in'

  # NOTE: volume parameters are added via command-line arguments by tests file.
  position = '1.02 0 0'

  initial_T = T0
  initial_vel_x = 1
  initial_vel_y = 0
  initial_vel_z = 0

  scaling_factor_rhoV  = 1
  scaling_factor_rhouV = 1
  scaling_factor_rhovV = 1
  scaling_factor_rhowV = 1
  scaling_factor_rhoEV = 1e-5
[]

[pipe2]
  type = FlowChannel1Phase

  position = '1.04 0 0'
  orientation = '1 0 0'
  length = 0.96

  initial_T = T0

  n_elems = 24
[]

[outlet]
  type = Outlet1Phase
  input = 'pipe2:out'
  p = 1e5
[]
