[inlet]
  type = InletDensityVelocity1Phase
  input = 'inlet_pipe:in'
  rho = 1.37931034483
  vel = 1
[]

[inlet_pipe]
  type = FlowChannel1Phase
  fp = eos

  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  A = 1

  f = 0

  initial_T = T0
  initial_p = 1e5
  initial_vel = 1

  n_elems = 20
[]

[junction1]
  type = VolumeJunction1Phase
  connections = 'inlet_pipe:out deadend_pipe:in outlet_pipe:in'
  position = '1 0 0'
  volume = 1e-8
[]

[outlet_pipe]
  type = FlowChannel1Phase
  fp = eos

  position = '1 0 0'
  orientation = '1 0 0'
  length = 1
  A = 1

  f = 0

  initial_T = 250
  initial_p = 1e5
  initial_vel = 1

  n_elems = 20
[]

[outlet]
  type = Outlet1Phase
  input = 'outlet_pipe:out'
  p = 1e5
[]

[deadend_pipe]
  type = FlowChannel1Phase
  fp = eos

  position = '1 0 0'
  orientation = '0 1 0'
  length = 1
  A = 1

  f = 0

  initial_T = 250
  initial_p = 1e5
  initial_vel = 0

  n_elems = 20
[]

[deadend]
  type = SolidWall1Phase
  input = 'deadend_pipe:out'
[]
