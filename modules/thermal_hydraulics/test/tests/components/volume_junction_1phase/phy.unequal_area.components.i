[inlet]
  type = InletDensityVelocity1Phase
  input = 'pipe1:in'
  rho = 1.37931034483
  vel = 1
[]

[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  A = 1
  n_elems = 20
[]

[junction]
  type = VolumeJunction1Phase
  connections = 'pipe1:out pipe2:in'
  position = '1 0 0'
  volume = 1e-8
[]

[pipe2]
  type = FlowChannel1Phase
  position = '1 0 0'
  orientation = '1 0 0'
  length = 1
  A = 0.5
  n_elems = 20
[]

[outlet]
  type = Outlet1Phase
  input = 'pipe2:out'
  p = 1e5
[]
