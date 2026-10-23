[left_boundary]
  type = FreeBoundary1Phase
  input = 'channel:in'
[]

[channel]
  type = FlowChannel1Phase

  fp = fp

  position = '0 0 0'
  orientation = '1 0 0'
  length = 1.0
  n_elems = 100
  A = 1.0

  initial_T = T_ic_fn
  initial_p = p_ic_fn
  initial_vel = 0

  f = 0
[]

[right_boundary]
  type = FreeBoundary1Phase
  input = 'channel:out'
[]
