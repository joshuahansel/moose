[left_boundary]
  type = InletDensityVelocity1Phase
  input = 'left_channel:in'
  rho = 466.6666667
  vel = 1
[]

[left_channel]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 0.5
  n_elems = 50
  A = 1.0
[]

[junction]
  type = JunctionOneToOne1Phase
  connections = 'left_channel:out right_channel:in'
[]

[right_channel]
  type = FlowChannel1Phase
  position = '0.5 0 0'
  orientation = '1 0 0'
  length = 0.5
  n_elems = 50
  A = 0.5
[]

[right_boundary]
  type = Outlet1Phase
  input = 'right_channel:out'
  p = 1e5
[]
