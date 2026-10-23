[left_boundary]
  type = FreeBoundary1Phase
  input = 'left_channel:in'
[]

[left_channel]
  type = FlowChannel1Phase

  fp = fp

  position = '0 0 0'
  orientation = '1 0 0'
  length = 0.5
  n_elems = 50
  A = 1.0

  initial_T = T_ic_fn
  initial_p = p_ic_fn
  initial_vel = 0

  f = 0
[]

[junction]
  type = JunctionOneToOne1Phase
  connections = 'left_channel:out right_channel:in'
[]

[right_channel]
  type = FlowChannel1Phase

  fp = fp

  position = '0.5 0 0'
  orientation = '1 0 0'
  length = 0.5
  n_elems = 50
  A = 1.0

  initial_T = T_ic_fn
  initial_p = p_ic_fn
  initial_vel = 0

  f = 0
[]

[right_boundary]
  type = FreeBoundary1Phase
  input = 'right_channel:out'
[]
