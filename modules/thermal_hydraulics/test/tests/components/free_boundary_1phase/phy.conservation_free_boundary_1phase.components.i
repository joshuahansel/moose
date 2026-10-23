[inlet]
  type = FreeBoundary1Phase
  input = pipe:in
[]
[pipe]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1.0
  n_elems = 10
  A = 1.0

  initial_T = T_fn
  initial_p = 1e5
  initial_vel = 1

  f = 0

  fp = fp
[]
[outlet]
  type = FreeBoundary1Phase
  input = pipe:out
[]
