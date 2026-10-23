[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 2
  A = 1
[]

[turbine]
  type = SimpleTurbine1Phase
  connections = 'pipe1:out pipe2:in'
  position = '1 0 0'
  volume = 1
  A_ref = 1.0
  K = 0
  on = false
  power = 1000
[]

[pipe2]
  type = FlowChannel1Phase
  position = '1. 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 2
  A = 1
[]
